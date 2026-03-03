"""
Complex Geometry: A package for calculating local geometrical features.

This package calculates local geometrical features for protein-ligand systems,
combining the best of MDAnalysis and Biotite libraries.

Usage:
    from complex_geometry import ProteinLigandGeometry

    geom = ProteinLigandGeometry('1A9M.pdb', chain_id='A')
    features = geom.calculate_all_features()
"""

__version__ = "0.1.0"
__author__ = "Your Name"
__license__ = "MIT"

from complex_geometry.core import ProteinLigandGeometry, analyze_trajectory

__all__ = ["ProteinLigandGeometry", "analyze_trajectory"]
