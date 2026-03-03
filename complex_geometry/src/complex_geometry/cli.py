"""
Command-line interface for complex-geometry package.
"""

import argparse
import json
import sys
from pathlib import Path

from complex_geometry import ProteinLigandGeometry


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Calculate local geometrical features for protein-ligand systems"
    )
    parser.add_argument(
        "pdb_file", type=str, help="Path to PDB file or PDB ID (will download)"
    )
    parser.add_argument(
        "-c", "--chain", type=str, default="A", help="Chain ID for protein (default: A)"
    )
    parser.add_argument(
        "-l",
        "--ligand",
        type=str,
        default=None,
        help="Ligand residue name (auto-detect if not specified)",
    )
    parser.add_argument(
        "-b",
        "--binding-cutoff",
        type=float,
        default=10.0,
        help="Distance cutoff for binding site (Angstrom, default: 10.0)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Output JSON file (print to stdout if not specified)",
    )
    parser.add_argument(
        "-f",
        "--features",
        type=str,
        nargs="+",
        default=None,
        help="Specific features to calculate (default: all)",
    )

    args = parser.parse_args()

    # Handle PDB ID or file path
    pdb_file = args.pdb_file
    if len(pdb_file) == 4 and pdb_file.isalnum():
        # Assume it's a PDB ID
        import urllib.request

        url = f"https://files.rcsb.org/download/{pdb_file}.pdb"
        pdb_file = f"{pdb_file}.pdb"
        if not Path(pdb_file).exists():
            print(f"Downloading PDB {args.pdb_file}...")
            urllib.request.urlretrieve(url, pdb_file)

    # Calculate features
    geom = ProteinLigandGeometry(
        pdb_file=pdb_file,
        chain_id=args.chain,
        ligand_resname=args.ligand,
        binding_cutoff=args.binding_cutoff,
    )

    features = geom.calculate_all_features()

    # Filter features if requested
    if args.features:
        features = {k: v for k, v in features.items() if k in args.features}

    # Output
    if args.output:
        with open(args.output, "w") as f:
            json.dump(features, f, indent=2)
        print(f"Results saved to {args.output}")
    else:
        print(json.dumps(features, indent=2))

    return 0


if __name__ == "__main__":
    sys.exit(main())
