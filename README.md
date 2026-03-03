# Complex Geometry

Unified implementation for calculating local geometrical features of protein-ligand systems, combining MDAnalysis and Biotite libraries.

## 1. Overview

This package calculates five local geometrical features that characterize protein-ligand binding affinity.

### Features Implemented

| Feature | Symbol | Description |
|---------|--------|-------------|
| Protein-ligand closeness | cn | Distance-based binding measure |
| Solvent accessible surface area | SASA | Surface area of binding site |
| Protein-ligand orientation | ot | Spatial orientation of ligand |
| Protein-ligand contacts | ct | Atomic contact count |
| Interfacial hydrogen bonds | hb | Hydrogen bond count |

---

## 2. Design Philosophy

Different libraries excel at different tasks:

| Function | Library | Reason |
|----------|---------|--------|
| Structure loading | Biotite | Faster parsing |
| Centroid calculation | Biotite | Optimized C implementation |
| SASA | Biotite | Shrake-Rupley algorithm |
| CellList neighbor search | Biotite | Efficient spatial indexing |
| Atom selection | MDAnalysis | Flexible selection syntax |
| Contacts | MDAnalysis | Mature selection engine |
| Hydrogen bonds | MDAnalysis | Better atom type handling |
| Trajectory support | MDAnalysis | Native trajectory support |

---

## 3. Feature Formulas

### 3.1 Protein-ligand Closeness (cn)

$$c_n = -\frac{\sum_{R_k \in BSR} D(\bar{X}_{R_k}, \bar{X}_{LIG})}{|BSR|}$$

- $BSR$ = binding site residues
- $\bar{X}_{R_k}$ = geometric center of residue $R_k$
- $\bar{X}_{LIG}$ = geometric center of ligand
- Measures negative average distance between binding site and ligand centers
- Higher values indicate tighter binding

### 3.2 Solvent Accessible Surface Area (SASA)

$$SASA = \sum_{i \in BS} SASA(i)$$

- Sum of SASA for all binding site atoms
- Uses Shrake-Rupley algorithm
- Lower SASA typically indicates tighter binding

### 3.3 Protein-ligand Orientation (ot)

$$o_t = \frac{\sum_{R_k \in BSR} |\frac{\pi}{2} - \angle LBB_k|}{|BSR|}$$

- Each angle $\angle LBB_k$ is formed by: binding site center → residue $R_k$ and binding site center → ligand
- Measures deviation from perpendicular orientation
- Lower values indicate more optimal orientation

### 3.4 Protein-ligand Contacts (ct)

$$c_t = \sum_{R_k \in BSR} \sum_{i \in R_k} \mathbb{I}[\min_{j \in LIG} d(i,j) < 3\mathring{A}]$$

- Counts atomic pairs within 3Å distance
- Uses indicator function $\mathbb{I}$

### 3.5 Interfacial Hydrogen Bonds (hb)

$$h_b = \sum_{R_k \in BSR} n_{D \in R_k, A \in LIG} + \sum_{R_k \in BSR} n_{D \in LIG, A \in R_k}$$

- Hydrogen bond criteria:
  - Donor-acceptor distance < 3.0 Å
  - Donor-H...Acceptor angle > 135°

---

## 4. Usage

### 4.1 Installation

```bash
# Clone repository
git clone git@github.com:MilesYyh/complex_geometry.git
cd complex-geometry

# Install dependencies
pip install numpy biotite mdanalysis

# Install in development mode
pip install -e .
```

### 4.2 Python API

```python
from complex_geometry import ProteinLigandGeometry

# Auto-detect chain
geom = ProteinLigandGeometry('protein.pdb')

# Specify chain
geom = ProteinLigandGeometry('protein.pdb', chain_id='A')

# Specify ligand name (if auto-detection fails)
geom = ProteinLigandGeometry('protein.pdb', ligand_resname='LIG')

# Adjust binding site cutoff (default 10.0 Å)
geom = ProteinLigandGeometry('protein.pdb', binding_cutoff=12.0)

# Calculate all features
features = geom.calculate_all_features()

print(features)
# {'closeness': -8.72, 'sasa': 6781.39, 'orientation': 0.50, ...}
```

### 4.3 Command Line

```bash
# Auto-detect chain
plg-calc structure.pdb

# Specify chain
plg-calc structure.pdb -c A

# Select specific features
plg-calc structure.pdb -c A --features closeness sasa contacts

# Output to JSON
plg-calc structure.pdb -c A -o results.json
```

---

## 5. Example

### 5.1 Provided Example

The `example/` folder contains PDB structure `1ZZT.pdp` (dihydrofolate reductase with substrate DP9) and sample output.

```bash
# Run example
cd example
python -c "
from complex_geometry import ProteinLigandGeometry
geom = ProteinLigandGeometry('1ZZT.pdb', ligand_resname='DP9')
features = geom.calculate_all_features()
print(features)
"

# Or use CLI
cpx-calc 1ZZT.pdb -l DP9
```

### 5.2 Example Output

**PDB: 1ZZT (Dihydrofolate reductase) + DP9 (substrate)**

```
closeness                     : -7.7819
sasa                          : 19381.1641
sasa_binding_site             : 1190.0000
orientation                   : 0.3126

contacts_total                : 5
contacts_hydrophobic          : 0
contacts_polar                : 5

hydrogen_bonds                : 0

com_distance                  : 18.6272
cog_distance                 : 18.6731
min_distance                  : 0.0000

n_protein_atoms               : 3220
n_ligand_atoms                : 46
n_binding_residues           : 15
```

---

## 6. Extended Applications

### 6.1 Mutation Delta Features

```python
def calculate_mutation_delta(features_wt, features_mut):
    """Calculate feature differences between wild-type and mutant"""
    delta = {}
    for key in features_wt:
        delta[f'delta_{key}'] = features_mut[key] - features_wt[key]
    return delta
```

### 6.2 Trajectory Analysis

```python
import MDAnalysis as mda

u = mda.Universe('structure.pdb', 'trajectory.xtc')

for ts in u.trajectory:
    features = geom.calculate_all_features()
    time_series.append(features)
```

For publishing instructions, see [release.md](release.md).

---

## 7. Appendix: Feature List

| Category | Feature | Type | Description |
|----------|---------|------|-------------|
| Geometry | closeness | float | Protein-ligand binding affinity |
| Geometry | orientation | float | Orientation deviation |
| Geometry | cog_distance | float | Geometric center distance |
| Geometry | com_distance | float | Center of mass distance |
| Geometry | min_distance | float | Minimum atom distance |
| Surface | sasa | float | Total protein SASA |
| Surface | sasa_binding_site | float | Binding site SASA |
| Contacts | contacts_total | int | Total contact count |
| Contacts | contacts_hydrophobic | int | Hydrophobic contacts |
| Contacts | contacts_polar | int | Polar contacts |
| H-bonds | hydrogen_bonds | int | H-bond count |
| Info | n_protein_atoms | int | Protein atom count |
| Info | n_ligand_atoms | int | Ligand atom count |
| Info | n_binding_residues | int | Binding site residues |
