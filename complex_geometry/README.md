# Complex Geometry

Calculate local geometrical features for protein-ligand systems using MDAnalysis and Biotite.

## Features

- **Protein-ligand closeness**: Distance-based binding affinity measure
- **SASA**: Solvent accessible surface area calculation
- **Orientation**: Protein-ligand spatial orientation
- **Contacts**: Atomic contact analysis (hydrophobic/polar)
- **Hydrogen bonds**: Interfacial hydrogen bond detection
- **Distance features**: Center of mass, geometric center, minimum distance
- **Trajectory analysis**: Support for MD trajectories

## Installation

```bash
pip install complex-geometry
```

## Usage

### As a Python library

```python
from complex_geometry import ProteinLigandGeometry

geom = ProteinLigandGeometry('1ZZT.pdb', ligand_resname='DP9')
features = geom.calculate_all_features()

print(features)
```

### Command-line interface

```bash
# Basic usage
cpx-calc 1ZZT.pdb -l DP9

# Specify features
cpx-calc 1ZZT.pdb -l DP9 --features closeness sasa contacts

# Output to file
cpx-calc 1ZZT.pdb -l DP9 -o results.json
```

### Trajectory Analysis

```python
from complex_geometry import analyze_trajectory

results = analyze_trajectory(
    topology='system.pdb',
    trajectory='traj.xtc',
    ligand_resname='DP9'
)
```

## Example

```bash
cd example
python -c "
from complex_geometry import ProteinLigandGeometry
geom = ProteinLigandGeometry('1ZZT.pdb', ligand_resname='DP9')
features = geom.calculate_all_features()
print(features)
"
```

Output:
```
closeness: -7.78
sasa: 19381.16
sasa_binding_site: 1190.0
orientation: 0.31
contacts_total: 5
hydrogen_bonds: 0
...
```

## Development

```bash
git clone https://github.com/yourusername/complex-geometry.git
cd complex-geometry
pip install -e ".[all]"
pytest
python -m build
```

## License

MIT License
