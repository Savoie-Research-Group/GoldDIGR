# GoldDIGR: Search, Download, Analyze Open-Source Chemistry Data

<div align="center">

![Logo](Logo.png)

</div>

Systematic, resumable mining of supplementary information (SI) files from
chemistry journal articles. Extracts XYZ coordinates, computational-chemistry
details, and BibLaTeX metadata — fully containerized, no manual setup.

## What It Does

Given a CSV of article URLs, golddigr:

1. Scrapes article HTML (handles Cloudflare challenges)
2. Extracts SI download links (direct files + publisher download endpoints)
3. Downloads SI files (PDF, XYZ, CIF, MOL, DOCX, XLSX, ZIP, …)
4. Converts PDFs to text via StirlingPDF (sidecar container)
5. Cleans footers/headers, extracts XYZ coordinate blocks
6. Detects transition-state geometries
7. Extracts BibLaTeX metadata from HTML
8. Uses an LLM to pull structured computational-chemistry details

Every stage is tracked in SQLite. Re-running picks up exactly where it left off.

## Architecture

```
┌──────────────────────────────────────────────────────────────────────┐
│                        docker compose up                            │
│                                                                     │
│  ┌──────────────┐         ┌───────────────────────────────────┐     │
│  │ StirlingPDF  │◄────────│       Python Pipeline             │     │
│  │  (sidecar)   │  HTTP   │                                   │     │
│  │              │         │  CSV → Scrape HTML                │     │
│  │  split-pages │         │        → Extract SI links         │     │
│  │  pdf-to-text │         │          → Download files         │     │
│  └──────────────┘         │            → Process PDFs ────────┤     │
│                           │              → Clean text         │     │
│                           │                → LLM extract      │     │
│                           │                  → Find TS XYZ    │     │
│                           └───────────────────────────────────┘     │
│                                        │                            │
│                              ┌─────────┴──────────┐                 │
│                              │   SQLite Ledger     │                 │
│                              │  (resumable state)  │                 │
│                              └────────────────────┘                 │
└──────────────────────────────────────────────────────────────────────┘
```

### Pipeline Stages

```
PENDING → HTML_SCRAPED → LINKS_EXTRACTED → FILES_DOWNLOADED
        → PDF_PROCESSED → TEXT_EXTRACTED → DONE
```

Any stage can fail independently without affecting other articles.

## Quick Start

### Option A: Docker Compose (simplest)

Best for batch processing when Cloudflare isn't an issue.

```bash
git clone https://github.com/your-username/golddigr.git && cd golddigr

# Configure API keys
cp .env.example .env
# Edit .env with your LLM API key

# Place your article CSV
mkdir -p data/input
cp your-articles.csv data/input/articles.csv

# Run
docker compose up --build

# Or process a specific range
docker compose run pipeline python run.py --start 0 --end 100

# Check progress
docker compose run pipeline python run.py --status

# Retry failed jobs
docker compose run pipeline python run.py --retry-failed
```

### Option B: Local Chrome (recommended for Cloudflare-heavy publishers)

Best for Wiley, Elsevier, and other publishers with aggressive bot detection.
Uses your real Chrome browser with your existing cookies and login sessions.

```bash
# 1. Start StirlingPDF
docker run -p 8080:8080 stirlingtools/stirling-pdf:latest

# 2. Create a Python environment
conda create -n golddigr python=3.11 && conda activate golddigr
# or: python -m venv .venv && source .venv/bin/activate

# 3. Install Python dependencies
pip install -r requirements.txt

# 4. Set environment variables
export STIRLING_URL=http://localhost:8080
export OPENAI_API_KEY=sk-proj-...

# 5. Use the Chrome config
cp config.chrome.yaml config.yaml

# 6. Run (Chrome will open with your profile)
python run.py --start 0 --end 100
```

Chrome launches via CDP (Chrome DevTools Protocol) using your real profile —
Cloudflare sees a legitimate browser with real cookies. If a CAPTCHA still
appears, the pipeline will prompt you to solve it manually and retry.

> **Note:** If you have multiple Chrome profiles, set `chrome_profile_directory`
> in config.yaml to match yours (e.g., `Default`, `Profile 1`). Check with:
> `ls ~/.config/google-chrome/`

> **Important:** Chrome opens PDFs in its built-in viewer by default instead of
> downloading them. You need to change this once:
> Chrome → Settings → Privacy and security → Site settings → Additional content
> settings → PDF documents → select **"Download PDFs"**

### Option C: Local Firefox

For websites without aggressive bot detection.

```bash
# 1. Start StirlingPDF
docker run -p 8080:8080 stirlingtools/stirling-pdf:latest

# 2. Create a Python environment
conda create -n golddigr python=3.11 && conda activate golddigr
# or: python -m venv .venv && source .venv/bin/activate

# 3. Install dependencies
pip install -r requirements.txt

# 4. Set environment variables
export STIRLING_URL=http://localhost:8080
export OPENAI_API_KEY=sk-proj-...

# 5. Use the Firefox config
cp config.firefox-local.yaml config.yaml

# 6. Run
python run.py --start 0 --end 100
```

### Option D: Singularity/Apptainer (HPC)

```bash
# Build the SIF (one-time)
docker build -t golddigr .
apptainer build golddigr.sif docker-daemon://golddigr:latest

# Run
export OPENAI_API_KEY=sk-proj-...
./golddigr --start 0 --end 100
./golddigr --status
```

## Output

```
data/
├── input/
│   └── articles.csv              # your input
├── output/
│   ├── text/                     # extracted text from PDFs
│   │   └── {doi_path}/{stem}_xyz_clean.txt
│   ├── xyz/                      # transition-state structures
│   ├── biblatex/                 # BibLaTeX entries per article
│   └── comp_details/             # LLM-extracted comp-chem JSON
├── downloads/
│   ├── html/                     # scraped article pages
│   └── files/{doi_path}/         # raw SI files (PDF, XYZ, CIF, …)
└── db/
    ├── pipeline.db               # SQLite job ledger
    └── pipeline.log              # full run log
```

## Sample Data

A `sample-data/` directory is included with a complete input/output example
from one article ([Guérard et al., *Chem. Eur. J.* 2016, 22, 12332](https://doi.org/10.1002/chem.201600922)):

```
sample-data/
├── input/articles.csv                          # one-article CSV
├── output/
│   ├── xyz/10.1002/chem.201600922/…/01.xyz     # 12 extracted XYZ structures
│   ├── biblatex/0.bib                          # BibLaTeX entry
│   ├── comp_details/…/…_comp.json              # LLM-extracted comp-chem details
│   ├── text/…/…_full.txt                       # raw extracted text
│   └── figures/…/manifest.json                 # figure metadata (images removed)
├── db/pipeline.db                              # SQLite ledger showing DONE status
└── sample.log                                  # pipeline log
```

Use this to understand the expected output format, or as a test:

```bash
# Copy sample input and run
cp sample-data/input/articles.csv data/input/articles.csv
python run.py --start 0 --end 1
```

## Configuration

Three example configs are provided:

| File | Browser | Profile | Headless | Use case |
|------|---------|---------|----------|----------|
| `config.yaml` | Firefox | Clean | Yes | Docker / Singularity |
| `config.chrome.yaml` | Chrome | Your profile | No | Cloudflare-heavy |
| `config.firefox-local.yaml` | Firefox | Clean | Yes | Local, no Cloudflare issues |

Copy whichever fits your use case to `config.yaml`.

### Key Settings

```yaml
scraper:
  browser: chrome          # "firefox" or "chrome"
  browser_profile: auto    # "none", "auto", or explicit path
  chrome_profile_directory: Default  # Chrome profile subdirectory
  headless: false          # true for containers, false to see the browser
  interactive: true        # prompt on Cloudflare challenges
```

## Cloudflare Handling

Many websites use Cloudflare bot detection. golddigr
handles this with a multi-tier approach:

1. **Automatic wait** — If a challenge page is detected, waits 5s → 10s → 15s
   for it to auto-clear (often works for IP-based challenges)

2. **Chrome + user profile** — Using `browser: chrome` with `browser_profile: auto`
   launches Chrome with your real cookies and browser fingerprint. Cloudflare
   sees a legitimate browser and usually skips the challenge entirely.

3. **Interactive fallback** — If challenges persist and `interactive: true`,
   the pipeline opens the URL in your default browser and prompts you to solve
   the CAPTCHA. Type Enter to retry or `s` to skip the article.

4. **Cached challenge detection** — If a previously saved HTML file contains a
   Cloudflare challenge page, it's automatically deleted and re-downloaded.

## Project Structure

```
golddigr/
├── docker-compose.yml          # Two-service setup (pipeline + StirlingPDF)
├── Dockerfile                  # Python + Firefox + geckodriver
├── config.yaml                 # Default config (Docker)
├── config.chrome.yaml          # Local Chrome config
├── config.firefox-local.yaml   # Local Firefox config
├── .env.example                # API key template
├── run.py                      # CLI entry point
├── golddigr                    # Singularity/Apptainer wrapper
├── pipeline/
│   ├── orchestrator.py         # State-machine driver
│   ├── job_db.py               # SQLite job ledger
│   ├── scraper.py              # Browser automation (Firefox + Chrome CDP)
│   ├── link_extractor.py       # SI link detection (direct + query-param URLs)
│   ├── stirling_client.py      # StirlingPDF Python wrapper
│   ├── pdf_processor.py        # PDF → text → cleanup → XYZ
│   ├── pdf_txt_processing.py   # Footer detection, header removal
│   ├── separate_xyz.py         # XYZ line detection, Gaussian conversion
│   ├── file_processors.py      # Route files by type (docx/xlsx/zip → text)
│   ├── metadata.py             # HTML meta → BibLaTeX
│   ├── figure_extractor.py     # PDF figure extraction
│   └── cc_detector.py          # Comp-chem keyword detection + LLM
└── data/                       # Mounted volumes (gitignored)
```

## Key Design Decisions

| Problem | Approach |
|---|---|
| Tracking progress | SQLite ledger with state machine |
| PDF processing | StirlingPDF sidecar, zero config |
| Cloudflare bypass | Chrome CDP with real user profile |
| Download retries | Watchdog thread + cookie strategy loop |
| Footer removal | NLTK-based detection in `pdf_txt_processing.py` |
| File format routing | Extension + query-param detection in `link_extractor.py` |
| Configuration | Single `config.yaml` + env var overrides |
| Resumability | Automatic — SQLite tracks every stage |

## License

MIT — see [LICENSE](LICENSE).
