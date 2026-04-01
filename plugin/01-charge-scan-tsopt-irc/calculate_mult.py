import sys

def get_multiplicity(xyz_file, charge):
    """
    Calculates the lowest spin multiplicity for a molecule.
    Multiplicity = 2S + 1.
    If total electrons is even, lowest energy state is usually a singlet (S=0, mult=1).
    If total electrons is odd, it must be at least a doublet (S=1/2, mult=2).
    """
    atomic_numbers = {
        'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
        'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19,
        'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28,
        'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
        'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46,
        'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55,
        'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64,
        'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73,
        'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82,
        'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86
    }

    total_electrons = 0
    try:
        with open(xyz_file, 'r') as f:
            lines = f.readlines()
            # The first two lines are count and comment
            for line in lines[2:]:
                parts = line.split()
                if len(parts) > 0:
                    symbol = line.split()[0].capitalize()
                    if symbol in atomic_numbers:
                        total_electrons += atomic_numbers[symbol]
                    else:
                        print(f"Error: Atom symbol '{symbol}' not found in dictionary.", file=sys.stderr)
                        sys.exit(1)

        # Adjust for charge
        total_electrons -= charge

        if total_electrons % 2 == 0:
            return 1  # Singlet
        else:
            return 2  # Doublet

    except FileNotFoundError:
        print(f"Error: XYZ file not found at '{xyz_file}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python calculate_mult.py <path_to_xyz_file> <charge>", file=sys.stderr)
        sys.exit(1)

    xyz_path = sys.argv[1]
    try:
        charge_val = int(sys.argv[2])
    except ValueError:
        print("Error: Charge must be an integer.", file=sys.stderr)
        sys.exit(1)

    multiplicity = get_multiplicity(xyz_path, charge_val)
    print(multiplicity)

