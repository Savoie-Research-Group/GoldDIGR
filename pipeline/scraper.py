"""
scraper.py – Selenium-based article scraping and file downloading.

Supports Firefox (default, for containers) and Chrome (for local use
with user profile to bypass Cloudflare).

Preserves the battle-tested logic from PID-test_requests_package.py:
  - Selenium-only downloads (handles Cloudflare, cookie redirections)
  - Watchdog thread + psutil process-tree killing for hung browsers
  - Sequential cookie-strategy retries with file-existence checks
  - Text-format files saved via driver.page_source
"""

from __future__ import annotations
import logging
import pathlib
import platform
import shutil
import subprocess
import sys
import webbrowser
import tempfile
import threading
import time
from pathlib import Path
from typing import Optional

import psutil
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.common.exceptions import TimeoutException
from urllib3.exceptions import ReadTimeoutError

from .link_extractor import deduce_filename

logger = logging.getLogger(__name__)

# ── Text formats that browsers open in-page (need page_source save) ────
TXT_FORMATS = [".txt", ".csv", ".xml", ".html", ".mol", ".mol2"]

# ── Default download timeout (seconds before watchdog force-kills) ─────
DEFAULT_TIMEOUT = 90

# ── Default cookie strategies to try sequentially ──────────────────────
DEFAULT_COOKIE_STRATEGIES = ["", "?cookieSet=0", "?cookieSet=1", "?cookieSet=2"]


# ═════════════════════════════════════════════════════════════════════════
# Browser profile detection
# ═════════════════════════════════════════════════════════════════════════

def _find_chrome_profile() -> Optional[str]:
    """
    Locate the user's Chrome/Chromium user-data-dir.
    Returns the parent directory of Default/, or None.
    """
    if platform.system() == "Darwin":
        candidates = [
            Path.home() / "Library" / "Application Support" / "Google" / "Chrome",
        ]
    elif platform.system() == "Windows":
        local = Path.home() / "AppData" / "Local"
        candidates = [local / "Google" / "Chrome" / "User Data"]
    else:
        candidates = [
            Path.home() / ".config" / "google-chrome",
            Path.home() / ".config" / "chromium",
            Path.home() / "snap" / "chromium" / "common" / "chromium",
            Path.home() / ".var" / "app" / "com.google.Chrome" / "config" / "google-chrome",
        ]

    for base in candidates:
        default_profile = base / "Default"
        if default_profile.is_dir():
            return str(base)
    return None


def _find_firefox_profile() -> Optional[str]:
    """
    Locate the user's default Firefox profile directory.
    Checks standard, snap, and flatpak installations.
    """
    import configparser

    if platform.system() == "Darwin":
        candidates = [Path.home() / "Library" / "Application Support" / "Firefox"]
    else:
        candidates = [
            Path.home() / ".mozilla" / "firefox",
            Path.home() / "snap" / "firefox" / "common" / ".mozilla" / "firefox",
            Path.home() / ".var" / "app" / "org.mozilla.firefox" / ".mozilla" / "firefox",
        ]

    for base in candidates:
        ini = base / "profiles.ini"
        if not ini.exists():
            continue

        config = configparser.ConfigParser()
        config.read(ini)

        for section in config.sections():
            if section.startswith("Install") and config.has_option(section, "Default"):
                rel = config.get(section, "Default")
                candidate = base / rel
                if candidate.is_dir():
                    return str(candidate)

        for section in config.sections():
            if config.has_option(section, "Default") and config.get(section, "Default") == "1":
                if config.has_option(section, "Path"):
                    rel = config.get(section, "Path")
                    is_relative = config.getboolean(section, "IsRelative", fallback=True)
                    candidate = base / rel if is_relative else Path(rel)
                    if candidate.is_dir():
                        return str(candidate)

        for section in config.sections():
            if config.has_option(section, "Path"):
                rel = config.get(section, "Path")
                is_relative = config.getboolean(section, "IsRelative", fallback=True)
                candidate = base / rel if is_relative else Path(rel)
                if candidate.is_dir():
                    return str(candidate)

    return None


def _resolve_profile(browser: str, profile_setting: Optional[str]) -> Optional[str]:
    """
    Resolve browser_profile config to a usable profile path.

    - None / "" / "none"  -> None (clean session)
    - "auto"              -> auto-detect profile for the chosen browser
    - "/path/to/profile"  -> use as-is

    For Chrome: returns a *copied* user-data-dir (safe to use while Chrome runs).
    For Firefox: returns the profile path (used for cookie injection only).
    """
    if not profile_setting or profile_setting.lower() == "none":
        return None

    if profile_setting.lower() == "auto":
        if browser == "chrome":
            found = _find_chrome_profile()
        else:
            found = _find_firefox_profile()
        if found:
            logger.info("Auto-detected %s profile: %s", browser, found)
        else:
            logger.warning("Could not auto-detect %s profile, using clean session", browser)
            return None
        profile_setting = found

    if not Path(profile_setting).is_dir():
        logger.warning("Profile not found: %s, using clean session", profile_setting)
        return None

    if browser == "chrome":
        return _copy_chrome_profile(profile_setting)
    else:
        cookies_db = Path(profile_setting) / "cookies.sqlite"
        if not cookies_db.is_file():
            logger.warning("No cookies.sqlite in %s, using clean session", profile_setting)
            return None
        return profile_setting

_cached_chrome_profile: Optional[str] = None

def _copy_chrome_profile(source: str) -> Optional[str]:
    """
    Copy a Chrome user-data-dir to a temp directory so Selenium can
    use it while the real Chrome is still running.

    Copies Default/ profile essentials, skips caches for speed.
    The copy is cached — subsequent calls reuse the same temp dir.
    """
    global _cached_chrome_profile
    if _cached_chrome_profile and Path(_cached_chrome_profile).is_dir():
        return _cached_chrome_profile

    try:
        tmp_dir = tempfile.mkdtemp(prefix="golddigr_chrome_")
        src = Path(source)

        # Local State has encryption keys for cookies
        local_state = src / "Local State"
        if local_state.is_file():
            shutil.copy2(str(local_state), tmp_dir)

        # Copy Default profile, skipping large cache dirs
        src_default = src / "Default"
        dst_default = Path(tmp_dir) / "Default"

        skip_dirs = {"Cache", "Code Cache", "GPUCache", "Service Worker",
                     "DawnCache", "ShaderCache", "GrShaderCache"}

        def _ignore(directory, contents):
            return [c for c in contents if c in skip_dirs]

        if src_default.is_dir():
            shutil.copytree(str(src_default), str(dst_default),
                            ignore=_ignore, dirs_exist_ok=True)

        # Remove lock/singleton files that prevent a second instance
        for lock_name in ["SingletonLock", "SingletonCookie", "SingletonSocket",
                          "lockfile", "Lock"]:
            lock = Path(tmp_dir) / lock_name
            if lock.exists() or lock.is_symlink():
                lock.unlink(missing_ok=True)
            # Also check inside Default/
            lock_inner = dst_default / lock_name
            if lock_inner.exists() or lock_inner.is_symlink():
                lock_inner.unlink(missing_ok=True)

        logger.info("Copied Chrome profile to temp dir: %s", tmp_dir)
        _cached_chrome_profile = tmp_dir
        return tmp_dir

    except Exception as exc:
        logger.warning("Failed to copy Chrome profile: %s, using clean session", exc)
        return None


def _inject_firefox_cookies(driver, profile_path: str, url: str) -> None:
    """
    Read cookies from a Firefox cookies.sqlite and inject them into
    the running Selenium session.
    """
    import sqlite3
    from urllib.parse import urlparse

    cookies_db = Path(profile_path) / "cookies.sqlite"
    if not cookies_db.is_file():
        return

    domain = urlparse(url).hostname
    tmp = tempfile.mktemp(suffix=".sqlite")
    shutil.copy2(str(cookies_db), tmp)

    try:
        conn = sqlite3.connect(tmp)
        conn.row_factory = sqlite3.Row
        parts = domain.split(".")
        domain_patterns = []
        for i in range(len(parts) - 1):
            d = ".".join(parts[i:])
            domain_patterns.append(d)
            domain_patterns.append("." + d)

        placeholders = ",".join("?" for _ in domain_patterns)
        rows = conn.execute(
            f"SELECT name, value, host, path, isSecure, expiry "
            f"FROM moz_cookies WHERE host IN ({placeholders})",
            domain_patterns,
        ).fetchall()
        conn.close()

        injected = 0
        for row in rows:
            cookie = {
                "name": row["name"],
                "value": row["value"],
                "domain": row["host"],
                "path": row["path"],
                "secure": bool(row["isSecure"]),
            }
            if row["expiry"]:
                cookie["expiry"] = row["expiry"]
            try:
                driver.add_cookie(cookie)
                injected += 1
            except Exception:
                pass
        if injected:
            logger.info("Injected %d cookies from Firefox profile for %s", injected, domain)
    except Exception as exc:
        logger.warning("Failed to read cookies from profile: %s", exc)
    finally:
        Path(tmp).unlink(missing_ok=True)


# ═════════════════════════════════════════════════════════════════════════
# Driver factory
# ═════════════════════════════════════════════════════════════════════════

def _get_firefox_service():
    """Return a geckodriver Service, trying common paths."""
    from selenium.webdriver.firefox.service import Service
    for path in ["/usr/local/bin/geckodriver", "/snap/bin/geckodriver"]:
        if Path(path).exists():
            return Service(path)
    return Service()


def _find_chrome_binary() -> Optional[str]:
    """Find the Chrome/Chromium binary on the system."""
    candidates = []
    if platform.system() == "Darwin":
        candidates = [
            "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
            "/Applications/Chromium.app/Contents/MacOS/Chromium",
        ]
    elif platform.system() == "Windows":
        for base in [Path.home() / "AppData" / "Local", Path("C:/Program Files"), Path("C:/Program Files (x86)")]:
            candidates.append(str(base / "Google" / "Chrome" / "Application" / "chrome.exe"))
    else:
        candidates = [
            "/usr/bin/google-chrome",
            "/usr/bin/google-chrome-stable",
            "/usr/bin/chromium",
            "/usr/bin/chromium-browser",
            "/snap/bin/chromium",
        ]
    for path in candidates:
        if Path(path).exists():
            return path
    # Try PATH
    result = shutil.which("google-chrome") or shutil.which("chromium") or shutil.which("chromium-browser")
    return result


# Track Chrome subprocesses so we can kill them on driver.quit()
_chrome_processes: dict[int, subprocess.Popen] = {}

_CHROME_DEBUG_PORT_START = 19222  # avoid conflicts with common ports


def _next_debug_port() -> int:
    """Find an available debug port."""
    import socket
    port = _CHROME_DEBUG_PORT_START
    while port < _CHROME_DEBUG_PORT_START + 100:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(("127.0.0.1", port)) != 0:
                return port
        port += 1
    return port


def _create_driver(
    browser: str = "firefox",
    *,
    headless: bool = True,
    download_dir: Optional[Path] = None,
    profile_path: Optional[str] = None,
    profile_directory: Optional[str] = None,
):
    """
    Create a Selenium WebDriver for the specified browser.

    For Chrome: launches the real Chrome binary via subprocess with
    remote debugging, then attaches Selenium via CDP. This avoids
    Selenium's automation flags that trigger bot detection.

    For Firefox: standard Selenium launch (for container use).

    Parameters
    ----------
    browser : "firefox" or "chrome"
    headless : bool
    download_dir : Path, optional
        Configure browser for silent file downloads to this directory.
    profile_path : str, optional
        Chrome: a user-data-dir path (original or copied).
        Firefox: ignored here (cookie injection happens separately).
    profile_directory : str, optional
        Chrome only: subdirectory name within user-data-dir
        (e.g. "Default", "Profile 1"). Skips the profile picker.
    """
    if browser == "chrome":
        from selenium.webdriver.chrome.options import Options as ChromeOptions

        chrome_bin = _find_chrome_binary()
        if not chrome_bin:
            raise RuntimeError(
                "Chrome/Chromium not found. Install google-chrome or chromium, "
                "or set browser: firefox in config.yaml"
            )

        debug_port = _next_debug_port()

        # Build the command to launch Chrome directly
        cmd = [
            chrome_bin,
            f"--remote-debugging-port={debug_port}",
            "--no-first-run",
            "--no-default-browser-check",
            "--disable-popup-blocking",
            "--disable-features=ChromePDF",  # force PDF download instead of viewer
            "--disable-extensions",          # prevent adblock welcome tabs
        ]
        if headless:
            cmd.append("--headless=new")
        if profile_path:
            cmd.append(f"--user-data-dir={profile_path}")
        if profile_directory:
            cmd.append(f"--profile-directory={profile_directory}")

        # Launch Chrome as a subprocess
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        # Wait for the debug port to be ready
        import socket
        for _ in range(30):
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                if s.connect_ex(("127.0.0.1", debug_port)) == 0:
                    break
            time.sleep(0.5)
        else:
            proc.kill()
            raise RuntimeError(f"Chrome did not start on debug port {debug_port}")

        # Give Chrome extra time to fully initialize the profile
        time.sleep(2)

        # Attach Selenium to the running Chrome
        opts = ChromeOptions()
        opts.add_experimental_option("debuggerAddress", f"127.0.0.1:{debug_port}")

        driver = webdriver.Chrome(options=opts)

        # Track the subprocess so we can kill it later
        _chrome_processes[id(driver)] = proc

        # Set download directory via CDP if needed
        if download_dir:
            # Force all downloads to save to disk with real filenames
            try:
                driver.execute_cdp_cmd("Browser.setDownloadBehavior", {
                    "behavior": "allow",
                    "downloadPath": str(download_dir),
                    "eventsEnabled": True,
                })
            except Exception:
                # Fallback for older Chrome versions
                driver.execute_cdp_cmd("Page.setDownloadBehavior", {
                    "behavior": "allow",
                    "downloadPath": str(download_dir),
                })

        logger.info("Chrome attached via CDP on port %d (pid %d)", debug_port, proc.pid)
        return driver

    else:  # firefox
        from selenium.webdriver.firefox.options import Options as FirefoxOptions

        opts = FirefoxOptions()
        if headless:
            opts.add_argument("--headless")

        if download_dir:
            opts.set_preference("browser.download.folderList", 2)
            opts.set_preference("browser.download.dir", str(download_dir))
            opts.set_preference("browser.download.useDownloadDir", True)
            opts.set_preference("browser.download.manager.showWhenStarting", False)
            mime_types = ",".join([
                "text/plain", "text/xml", "text/csv", "text/html",
                "application/xhtml+xml", "application/octet-stream",
            ])
            opts.set_preference("browser.helperApps.neverAsk.saveToDisk", mime_types)
            opts.set_preference("browser.helperApps.neverAsk.openFile", mime_types)
            opts.set_preference("browser.download.viewableInternally.enabledTypes", "")
            opts.set_preference("browser.download.manager.showAlertOnComplete", False)
            opts.set_preference("browser.download.manager.closeWhenDone", True)
            opts.set_preference("browser.download.alwaysOpenPanel", False)
            opts.set_preference("browser.download.manager.useWindow", False)
            opts.set_preference("browser.download.manager.focusWhenStarting", False)
            opts.set_preference("dom.disable_beforeunload", True)
            opts.set_preference("browser.download.manager.quitBehavior", 2)
            opts.set_preference("pdfjs.disabled", True)
            opts.page_load_strategy = "eager"

        driver = webdriver.Firefox(options=opts, service=_get_firefox_service())
        return driver


def _quit_driver(driver) -> None:
    """
    Quit the Selenium driver and kill any associated Chrome subprocess.
    Safe to call even if the driver is already dead.
    """
    try:
        proc = _chrome_processes.pop(id(driver), None)
        try:
            driver.quit()
        except Exception:
            pass
        if proc and proc.poll() is None:
            proc.terminate()
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                proc.kill()
    except Exception:
        pass


def _safe_navigate(driver, url: str, timeout: int = 20) -> None:
    """
    Navigate to *url* without hanging forever.

    For Chrome (CDP): uses Page.navigate (returns immediately) then
    waits up to *timeout* seconds, then calls Page.stopLoading —
    the programmatic equivalent of clicking the X button.

    For Firefox: uses driver.get() in a background thread with a
    hard timeout fallback.
    """
    # Try CDP non-blocking navigation first (Chrome)
    try:
        driver.execute_cdp_cmd("Page.navigate", {"url": url})
        # Wait for the page to become usable, then force-stop
        time.sleep(min(timeout, 15))
        try:
            driver.execute_cdp_cmd("Page.stopLoading", {})
        except Exception:
            pass
        return
    except Exception:
        pass  # Not Chrome/CDP — fall back to threaded approach

    # Firefox fallback: background thread with hard timeout
    nav_done = threading.Event()

    def _nav():
        try:
            driver.get(url)
        except Exception:
            pass
        finally:
            nav_done.set()

    t = threading.Thread(target=_nav, daemon=True)
    t.start()
    if not nav_done.wait(timeout):
        logger.debug("Navigation timeout after %ds (proceeding anyway): %s", timeout, url)
        try:
            driver.execute_script("window.stop();")
        except Exception:
            pass


# ═════════════════════════════════════════════════════════════════════════
# Process-tree killing
# ═════════════════════════════════════════════════════════════════════════
def kill_process_tree(root_pid: int, grace: float = 2.0):
    """Terminate *root_pid* and all its children."""
    try:
        parent = psutil.Process(root_pid)
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        return

    children = parent.children(recursive=True)
    for p in children:
        try:
            p.terminate()
        except psutil.NoSuchProcess:
            pass
    parent.terminate()

    gone, alive = psutil.wait_procs([parent, *children], timeout=grace)
    for p in alive:
        try:
            p.kill()
        except psutil.NoSuchProcess:
            pass

    if alive and platform.system() == "Windows":
        subprocess.run(
            ["taskkill", "/PID", str(root_pid), "/T", "/F"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )


# ═════════════════════════════════════════════════════════════════════════
# Watchdog + download worker
# ═════════════════════════════════════════════════════════════════════════
def _watchdog(done_event: threading.Event, pid_holder: dict, timeout: int):
    """Kill the browser if the download doesn't finish in time."""
    if done_event.wait(timeout):
        return
    pid = pid_holder.get("pid")
    if pid:
        logger.warning("Download exceeded %ds → killing process tree rooted at %d", timeout, pid)
        kill_process_tree(pid)


def _download_worker(
    browser: str,
    download_dir: Path,
    link: str,
    pid_holder: dict,
    done_event: threading.Event,
    headless: bool = True,
    page_load_timeout: int = 15,
    profile_path: Optional[str] = None,
    chrome_profile_directory: Optional[str] = None,
    agent_config: Optional[dict] = None,
):
    """
    Download worker — runs in its own thread.

    Launches a browser, navigates to *link*, and lets the download
    manager save the file. For text-format files, saves page_source.

    If a Cloudflare challenge is detected, waits for it to clear,
    then tries agent click if enabled, before re-navigating.

    Signals *done_event* when finished (or on error).
    """
    start = time.time()
    new_driver = None
    try:
        new_driver = _create_driver(
            browser,
            headless=headless,
            download_dir=download_dir,
            profile_path=profile_path if browser == "chrome" else None,
            profile_directory=chrome_profile_directory,
        )
        # Get the PID for the watchdog to kill if needed
        proc = _chrome_processes.get(id(new_driver))
        if proc:
            pid_holder["pid"] = proc.pid
        else:
            try:
                pid_holder["pid"] = new_driver.service.process.pid
            except Exception:
                pass

        # Firefox cookie injection (Chrome uses profile directly)
        if browser == "firefox" and profile_path:
            from urllib.parse import urlparse
            parsed = urlparse(link)
            base_url = f"{parsed.scheme}://{parsed.hostname}/"
            _safe_navigate(new_driver, base_url, timeout=page_load_timeout)
            _inject_firefox_cookies(new_driver, profile_path, link)

        _safe_navigate(new_driver, link, timeout=page_load_timeout)

        # Check for Cloudflare challenge and wait for it to clear
        try:
            html = new_driver.page_source
            agent_enabled = agent_config and agent_config.get("enabled")
            challenge_waits = [5] if agent_enabled else [5, 10, 15]
            for attempt, wait in enumerate(challenge_waits, 1):
                if not _is_bot_challenge(html):
                    break
                logger.info(
                    "Download: bot challenge (attempt %d/%d), waiting %ds… %s",
                    attempt, len(challenge_waits), wait, link,
                )
                time.sleep(wait)
                html = new_driver.page_source

            # Try agent click if challenge persists
            if _is_bot_challenge(html) and agent_enabled:
                try:
                    from .agent import attempt_click
                    if attempt_click(new_driver, agent_config):
                        html = new_driver.page_source
                except Exception as exc:
                    logger.warning("Download: agent click failed: %s", exc)

            # If challenge cleared, re-navigate to trigger the actual download
            if not _is_bot_challenge(html):
                _safe_navigate(new_driver, link, timeout=page_load_timeout)
                # Give Chrome time to start and complete the download
                time.sleep(10)
        except Exception:
            pass

    except TimeoutException:
        logger.debug("Page load timeout (expected for downloads): %s", link)
    except Exception:
        pass
    finally:
        try:
            with open(str(download_dir / "links.log"), "a", encoding="utf-8") as log:
                log.write(f"{link}\n")

            # For text-format files that the browser opens in-page, save page source
            fname = deduce_filename(link)
            file_ext = Path(fname).suffix.lower()
            if file_ext in TXT_FORMATS:
                output = download_dir / fname
                html = new_driver.page_source
                if html and not output.is_file():
                    pathlib.Path(output).write_text(html, encoding="utf-8")
                    logger.info("Saved browser-opened file → %s", output.name)

            _quit_driver(new_driver)
        except Exception:
            try:
                if new_driver:
                    _quit_driver(new_driver)
            except Exception:
                pass
        finally:
            done_event.set()

    elapsed = time.time() - start
    logger.debug("Download worker took %.1fs for %s", elapsed, link)


# ═════════════════════════════════════════════════════════════════════════
# File existence check
# ═════════════════════════════════════════════════════════════════════════
def check_file_and_size(actual_file: Path) -> bool:
    """Check if *actual_file* exists and is non-empty. Deletes 0-byte files."""
    if actual_file.is_file() and actual_file.stat().st_size == 0:
        actual_file.unlink()
        logger.debug("Deleted empty file: %s", actual_file)
    return actual_file.is_file()


# ═════════════════════════════════════════════════════════════════════════
# Bot / Cloudflare challenge detection
# ═════════════════════════════════════════════════════════════════════════
_BOT_CHALLENGE_MARKERS = [
    "you are human",
    "verify you are human",
    "verifies you are not a bot",
    "protect against malicious bots",
    "performing security verification",
    "security verification",
    "security service to protect",
    "checking your browser",
    "checking if the site connection is secure",
    "enable javascript and cookies to continue",
    "cf-challenge",
    "managed-challenge",
    "just a moment",
    "attention required",
    "please wait while we verify",
    "ray id",
]


def _is_bot_challenge(html: str) -> bool:
    """Return True if *html* looks like a bot-protection challenge page."""
    lower = html.lower()
    return any(marker in lower for marker in _BOT_CHALLENGE_MARKERS)


# ═════════════════════════════════════════════════════════════════════════
# Public API: scrape_html
# ═════════════════════════════════════════════════════════════════════════
def scrape_html(
    url: str,
    output_path: Path,
    *,
    headless: bool = True,
    timeout: int = 15,
    interactive: bool = False,
    browser: str = "firefox",
    browser_profile: Optional[str] = None,
    chrome_profile_directory: Optional[str] = None,
    agent_config: Optional[dict] = None,
    # Legacy alias (ignored if browser_profile is set)
    firefox_profile: Optional[str] = None,
) -> Path:
    """
    Navigate to *url* and save the rendered HTML.

    Parameters
    ----------
    browser : str
        "firefox" (default, for containers) or "chrome" (local use).
    browser_profile : str
        "none" (clean session), "auto" (detect user profile), or path.
    chrome_profile_directory : str
        Chrome profile subdirectory (e.g. "Default", "Profile 1").
        Skips the profile picker when multiple profiles exist.
    interactive : bool
        If True and challenges persist, open URL in user's default
        browser for manual solving with retry/skip prompt.
    """
    # Legacy compat: firefox_profile -> browser_profile
    if browser_profile is None and firefox_profile is not None:
        browser_profile = firefox_profile

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # If HTML already exists, check for bot challenge
    if output_path.is_file():
        existing = output_path.read_text(encoding="utf-8", errors="ignore")
        if _is_bot_challenge(existing):
            logger.info("Bot challenge detected in cached file, re-downloading: %s", output_path)
            output_path.unlink(missing_ok=True)
        else:
            return output_path

    profile_path = _resolve_profile(browser, browser_profile)

    driver = _create_driver(
        browser, headless=headless,
        profile_path=profile_path if browser == "chrome" else None,
        profile_directory=chrome_profile_directory,
    )

    try:
        _safe_navigate(driver, url, timeout=timeout)

        # Firefox cookie injection (Chrome uses profile directly)
        if browser == "firefox" and profile_path:
            _inject_firefox_cookies(driver, profile_path, url)
            _safe_navigate(driver, url, timeout=timeout)

        # Brief wait for dynamic content
        try:
            WebDriverWait(driver, min(timeout, 10)).until(
                lambda d: d.execute_script("return document.readyState") in ("complete", "interactive")
            )
        except (TimeoutException, Exception):
            pass

        html = driver.page_source

        # IP block (fatal)
        if "IP Address Blocked" in html:
            raise RuntimeError(f"IP Address Blocked detected for {url}")

        # Bot challenge — wait for it to clear
        # If agent is enabled, just one quick wait then hand off to agent.
        # If no agent, try the full 5s → 10s → 15s progression.
        agent_enabled = agent_config and agent_config.get("enabled")
        challenge_waits = [5] if agent_enabled else [5, 10, 15]
        for attempt, wait in enumerate(challenge_waits, 1):
            if not _is_bot_challenge(html):
                break
            logger.info(
                "Bot challenge detected (attempt %d/%d), waiting %ds… %s",
                attempt, len(challenge_waits), wait, url,
            )
            time.sleep(wait)
            html = driver.page_source

        # Automatic waits didn't clear it — try agent click
        if _is_bot_challenge(html) and agent_config and agent_config.get("enabled"):
            try:
                from .agent import attempt_click
                if attempt_click(driver, agent_config):
                    html = driver.page_source
            except Exception as exc:
                # Distinguish real missing module from model loading errors
                exc_str = str(exc)
                if "No module named" in exc_str and "agent" in exc_str:
                    logger.warning(
                        "Agent module not available. Install with: "
                        "pip install -r requirements-agent.txt"
                    )
                else:
                    logger.warning("Agent click failed: %s", exc)

        # Agent didn't clear it — fall back to interactive
        if _is_bot_challenge(html):
            if interactive:
                logger.info("Opening in default browser for manual challenge solving: %s", url)
                webbrowser.open(url)
                print(
                    "\n╔══════════════════════════════════════════════════════════╗\n"
                    "║  Bot challenge detected!                                ║\n"
                    "║  Opened the URL in your browser.                        ║\n"
                    "║  Solve the CAPTCHA, then come back here.                ║\n"
                    "║                                                         ║\n"
                    "║  Press Enter to retry  |  Type 's' to skip this article ║\n"
                    "╚══════════════════════════════════════════════════════════╝"
                )
                while _is_bot_challenge(html):
                    choice = input("  ▸ ").strip().lower()
                    if choice == "s":
                        raise RuntimeError(f"Bot challenge skipped by user for {url}")
                    # Retry with a fresh session
                    _quit_driver(driver)
                    profile_path = _resolve_profile(browser, browser_profile)
                    driver = _create_driver(
                        browser, headless=headless,
                        profile_path=profile_path if browser == "chrome" else None,
                        profile_directory=chrome_profile_directory,
                    )
                    _safe_navigate(driver, url, timeout=timeout)
                    if browser == "firefox" and profile_path:
                        _inject_firefox_cookies(driver, profile_path, url)
                        _safe_navigate(driver, url, timeout=timeout)
                    time.sleep(3)  # give page time to settle
                    html = driver.page_source
                    if _is_bot_challenge(html):
                        print("  ✗ Still blocked. Press Enter to retry or 's' to skip.")
                    else:
                        print("  ✓ Challenge cleared!")
            else:
                raise RuntimeError(f"Bot-protection challenge persisted for {url}")

        pathlib.Path(output_path).write_text(html, encoding="utf-8")
        logger.info("Saved HTML → %s", output_path)

        return output_path

    finally:
        _quit_driver(driver)


# ═════════════════════════════════════════════════════════════════════════
# Public API: download_file (with watchdog + cookie strategies)
# ═════════════════════════════════════════════════════════════════════════
def download_file(
    url: str,
    dest_dir: Path,
    *,
    headless: bool = True,
    timeout: int = DEFAULT_TIMEOUT,
    strategies: Optional[list[str]] = None,
    delay_between_strategies: float = 5.0,
    browser: str = "firefox",
    browser_profile: Optional[str] = None,
    chrome_profile_directory: Optional[str] = None,
    agent_config: Optional[dict] = None,
    # Legacy alias
    firefox_profile: Optional[str] = None,
) -> Optional[Path]:
    """
    Download a file from *url* into *dest_dir* using Selenium.

    Uses the watchdog + worker thread pattern:
      1. Spawn download_worker thread (launches browser, navigates to URL)
      2. Spawn watchdog thread (kills browser if timeout exceeded)
      3. Join both threads
      4. Check if file appeared
      5. If not, try next cookie strategy
    """
    # Legacy compat
    if browser_profile is None and firefox_profile is not None:
        browser_profile = firefox_profile

    strategies = strategies or DEFAULT_COOKIE_STRATEGIES
    dest_dir.mkdir(parents=True, exist_ok=True)
    profile_path = _resolve_profile(browser, browser_profile)

    # Determine expected filename from URL (handles query-param filenames)
    expected_name = deduce_filename(url)
    actual_file = dest_dir / expected_name

    # Already downloaded?
    if check_file_and_size(actual_file):
        return actual_file

    for i, strategy in enumerate(strategies):
        full_url = url + strategy
        logger.info("Trying download: %s (strategy %d/%d)", full_url, i + 1, len(strategies))

        # Snapshot existing files before download attempt
        existing_files = set(dest_dir.iterdir()) if dest_dir.is_dir() else set()

        done_event = threading.Event()
        pid_holder = {"pid": None}

        t_worker = threading.Thread(
            target=_download_worker,
            args=(browser, dest_dir, full_url, pid_holder, done_event, headless),
            kwargs={
                "profile_path": profile_path,
                "chrome_profile_directory": chrome_profile_directory,
                "agent_config": agent_config,
            },
            daemon=True,
        )
        t_watchdog = threading.Thread(
            target=_watchdog,
            args=(done_event, pid_holder, timeout),
            daemon=True,
        )

        t_worker.start()
        t_watchdog.start()

        t_worker.join()
        t_watchdog.join()

        # Wait a moment for Chrome to finish writing the file
        time.sleep(2)

        # Check 1: expected filename exists
        if check_file_and_size(actual_file):
            return actual_file

        # Check 2: scan for any new files that appeared (Chrome may
        # use Content-Disposition filename or append numbers)
        current_files = set(dest_dir.iterdir()) if dest_dir.is_dir() else set()
        new_files = current_files - existing_files
        # Ignore .crdownload (partial), links.log, and 0-byte files
        new_files = {
            f for f in new_files
            if f.is_file()
            and not f.name.endswith(".crdownload")
            and f.name != "links.log"
            and f.stat().st_size > 0
        }
        if new_files:
            found = sorted(new_files)[0]  # pick the first real file
            logger.info("Downloaded file detected: %s", found.name)
            return found

        if i < len(strategies) - 1:
            time.sleep(delay_between_strategies)

    logger.warning("All download strategies failed for %s", url)
    return None
