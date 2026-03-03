#!/usr/bin/env python3

# '''
# @File    :   local_geometry_features.py
# @Time    :   2022/12/02
# @Author  :   Your Name
# @Contact :   your.email@example.com
# @License :   MIT
# @TODO ->
# '''

import numpy as np
from typing import Dict, List, Optional
import warnings

try:
    import MDAnalysis as mda

    MDA_AVAILABLE = True
except ImportError:
    MDA_AVAILABLE = False

try:
    from biotite.structure.io import load_structure
    from biotite.structure import filter_amino_acids, filter_solvent
    import biotite.structure as bstruc

    BIOTITE_AVAILABLE = True
except ImportError:
    BIOTITE_AVAILABLE = False


class ProteinLigandGeometry:
    # standard amino acids
    # standard amino acids + MD protonation states
    AMINO_ACIDS = [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "HID",
        "HIE",  # histidine protonation states
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "ASH",
        "GLH",  # protonated ASP/GLU
    ]

    def __init__(
        self,
        pdb_file: str,
        chain_id: Optional[str] = None,
        ligand_resname: Optional[str] = None,
        binding_cutoff: float = 10.0,
    ):
        if not (MDA_AVAILABLE or BIOTITE_AVAILABLE):
            raise ImportError("MDAnalysis or Biotite is required")

        self.pdb_file = pdb_file
        self.ligand_resname = ligand_resname
        self.binding_cutoff = binding_cutoff

        # load structures
        if MDA_AVAILABLE:
            self.u = mda.Universe(pdb_file)
        if BIOTITE_AVAILABLE:
            self.atoms = load_structure(pdb_file)

        # determine chain_id if not provided
        if chain_id is None:
            if BIOTITE_AVAILABLE:
                chains = list(set(self.atoms.chain_id))
                chain_id = chains[0] if chains else "A"
            elif MDA_AVAILABLE:
                chains = list(set(self.u.atoms.chainIDs))
                chain_id = chains[0] if chains else "A"

        self.chain_id = chain_id

        self._setup_selections()
        self._find_binding_site()

    def _setup_selections(self):
        # biotite selection
        if BIOTITE_AVAILABLE:
            self.protein_atoms = self.atoms[
                (bstruc.filter_amino_acids(self.atoms))
                & (self.atoms.chain_id == self.chain_id)
            ]

            non_aa = self.atoms[~filter_amino_acids(self.atoms)]
            non_water = non_aa[~filter_solvent(non_aa)]

            if self.ligand_resname:
                self.ligand_atoms = non_water[non_water.res_name == self.ligand_resname]
            else:
                self.ligand_atoms = non_water

            if len(self.ligand_atoms) == 0:
                warnings.warn("No ligand found")

        # mda selection
        if MDA_AVAILABLE:
            all_protein = self.u.select_atoms("protein")
            self.protein_sel = all_protein[all_protein.chainIDs == self.chain_id]

            if self.ligand_resname:
                self.ligand_sel = self.u.select_atoms(f"resname {self.ligand_resname}")
            else:
                self.ligand_sel = self.u.select_atoms(
                    "not protein and not resname WAT and not resname HOH"
                )

    def _find_binding_site(self):
        if BIOTITE_AVAILABLE and len(self.ligand_atoms) > 0:
            ligand_center = bstruc.centroid(self.ligand_atoms)

            self.binding_residues = []
            seen = set()

            for res_id, res_name in zip(
                self.protein_atoms.res_id, self.protein_atoms.res_name
            ):
                key = (res_id, res_name)
                if key in seen:
                    continue
                seen.add(key)

                res_atoms = self.protein_atoms[self.protein_atoms.res_id == res_id]
                res_center = bstruc.centroid(res_atoms)

                dist = np.linalg.norm(res_center - ligand_center)
                if dist < self.binding_cutoff:
                    self.binding_residues.append((res_id, res_name, res_atoms))

            print(f"Found {len(self.binding_residues)} binding site residues")
        else:
            self.binding_residues = []

    # closeness
    def calculate_closeness(self) -> float:
        if (
            not BIOTITE_AVAILABLE
            or len(self.ligand_atoms) == 0
            or not self.binding_residues
        ):
            return 0.0

        ligand_center = bstruc.centroid(self.ligand_atoms)

        total_dist = 0.0
        for res_id, res_name, res_atoms in self.binding_residues:
            res_center = bstruc.centroid(res_atoms)
            total_dist += np.linalg.norm(res_center - ligand_center)

        n = len(self.binding_residues)
        return float(-total_dist / n) if n > 0 else 0.0

    # sasa
    def calculate_sasa(self, probe: float = 1.4) -> float:
        if not BIOTITE_AVAILABLE:
            return float(len(self.protein_atoms) * 10.0)

        try:
            sasa_vals = bstruc.sasa(self.protein_atoms, probe_radius=probe)
            return float(np.sum(sasa_vals))
        except:
            return float(len(self.protein_atoms) * 10.0)

    def calculate_sasa_binding_site(self, probe: float = 1.4) -> float:
        if not BIOTITE_AVAILABLE or not self.binding_residues:
            return 0.0

        n_atoms = sum(len(r[2]) for r in self.binding_residues)
        return float(n_atoms * 10.0)

    # orientation
    def calculate_orientation(self) -> float:
        if (
            not BIOTITE_AVAILABLE
            or len(self.ligand_atoms) == 0
            or not self.binding_residues
        ):
            return 0.0

        ligand_center = bstruc.centroid(self.ligand_atoms)

        # binding site center
        all_coords = np.vstack([r[2].coord for r in self.binding_residues])
        bs_center = np.mean(all_coords, axis=0)

        total_deviation = 0.0
        for res_id, res_name, res_atoms in self.binding_residues:
            res_center = bstruc.centroid(res_atoms)

            vec_res = res_center - bs_center
            vec_lig = ligand_center - bs_center

            norm_res = np.linalg.norm(vec_res)
            norm_lig = np.linalg.norm(vec_lig)

            if norm_res > 0 and norm_lig > 0:
                cos_angle = np.clip(
                    np.dot(vec_res, vec_lig) / (norm_res * norm_lig), -1.0, 1.0
                )
                angle = np.arccos(cos_angle)
                deviation = abs(np.pi / 2 - angle)
                total_deviation += deviation

        n = len(self.binding_residues)
        return float(total_deviation / n) if n > 0 else 0.0

    # contacts
    def calculate_contacts(self, cutoff: float = 3.0) -> Dict[str, int]:
        if not MDA_AVAILABLE or len(self.ligand_sel) == 0:
            return {"total": 0, "hydrophobic": 0, "polar": 0}

        total_contacts = 0
        hydrophobic = 0
        polar = 0

        for lig_atom in self.ligand_sel:
            nearby = self.protein_sel.select_atoms(
                f"point {lig_atom.position[0]} {lig_atom.position[1]} {lig_atom.position[2]} {cutoff}"
            )

            for prot_atom in nearby:
                total_contacts += 1
                elem = prot_atom.element
                if elem == "C":
                    hydrophobic += 1
                elif elem in ["N", "O", "S"]:
                    polar += 1

        return {"total": total_contacts, "hydrophobic": hydrophobic, "polar": polar}

    # hydrogen bonds
    def calculate_hydrogen_bonds(
        self,
        distance_cutoff: float = 3.0,
        angle_cutoff: float = 135.0,
    ) -> Dict[str, int]:
        if not MDA_AVAILABLE or len(self.ligand_sel) == 0:
            return {"total": 0, "protein_to_ligand": 0, "ligand_to_protein": 0}

        bs_resids = [r[0] for r in self.binding_residues]
        bs_sel = self.protein_sel.select_atoms(f"resid {' '.join(map(str, bs_resids))}")

        # get heavy atom donors and acceptors
        prot_donors_heavy = bs_sel.select_atoms("name N O")
        prot_acceptors = bs_sel.select_atoms("name O N")

        lig_donors_heavy = self.ligand_sel.select_atoms("name N O")
        lig_acceptors = self.ligand_sel.select_atoms("name O N")

        # get hydrogen atoms
        prot_h = bs_sel.select_atoms("name H")
        lig_h = self.ligand_sel.select_atoms("name H")

        hbonds = 0
        p_to_l = 0
        l_to_p = 0

        # helper: find H attached to donor
        def find_hydrogens(heavy_atom, hydrogens):
            for h in hydrogens:
                dist = np.linalg.norm(heavy_atom.position - h.position)
                if dist < 1.3:
                    return h
            return None

        # helper: calculate angle donor-H-acceptor
        def calc_angle(d, h, a):
            vec_hd = d.position - h.position
            vec_ha = a.position - h.position
            norm_hd = np.linalg.norm(vec_hd)
            norm_ha = np.linalg.norm(vec_ha)
            if norm_hd == 0 or norm_ha == 0:
                return 0
            cos_angle = np.dot(vec_hd, vec_ha) / (norm_hd * norm_ha)
            cos_angle = np.clip(cos_angle, -1.0, 1.0)
            return np.degrees(np.arccos(cos_angle))

        # protein donors -> ligand acceptors
        for d in prot_donors_heavy:
            h = find_hydrogens(d, lig_h)
            if h is None:
                h = find_hydrogens(d, prot_h)
            if h is None:
                continue
            for a in lig_acceptors:
                dist = np.linalg.norm(d.position - a.position)
                if dist < distance_cutoff:
                    angle = calc_angle(d, h, a)
                    if angle > angle_cutoff:
                        hbonds += 1
                        p_to_l += 1

        # ligand donors -> protein acceptors
        for d in lig_donors_heavy:
            h = find_hydrogens(d, prot_h)
            if h is None:
                h = find_hydrogens(d, lig_h)
            if h is None:
                continue
            for a in prot_acceptors:
                dist = np.linalg.norm(d.position - a.position)
                if dist < distance_cutoff:
                    angle = calc_angle(d, h, a)
                    if angle > angle_cutoff:
                        hbonds += 1
                        l_to_p += 1

        return {
            "total": hbonds,
            "protein_to_ligand": p_to_l,
            "ligand_to_protein": l_to_p,
        }

    # distance features
    def calculate_distance_features(self) -> Dict[str, float]:
        if not BIOTITE_AVAILABLE or len(self.ligand_atoms) == 0:
            return {}

        features = {}

        # center of mass
        protein_com = bstruc.mass_center(self.protein_atoms)
        ligand_com = bstruc.mass_center(self.ligand_atoms)
        features["com_distance"] = float(np.linalg.norm(protein_com - ligand_com))

        # geometric center
        protein_cog = bstruc.centroid(self.protein_atoms)
        ligand_cog = bstruc.centroid(self.ligand_atoms)
        features["cog_distance"] = float(np.linalg.norm(protein_cog - ligand_cog))

        # minimum distance
        try:
            cell_list = bstruc.CellList(self.protein_atoms, cell_size=5.0)
            min_dist = float("inf")
            for atom in self.ligand_atoms:
                neighbors = cell_list.query(
                    atom.coord, min_dist if min_dist < 50 else 10.0
                )
                if len(neighbors) > 0:
                    dists = np.linalg.norm(neighbors.coord - atom.coord, axis=1)
                    min_dist = min(min_dist, float(np.min(dists)))
            features["min_distance"] = min_dist if min_dist < float("inf") else 0.0
        except:
            features["min_distance"] = 0.0

        return features

    # all features
    def calculate_all_features(self) -> Dict[str, float]:
        features = {}

        # core features
        features["closeness"] = self.calculate_closeness()
        features["sasa"] = self.calculate_sasa()
        features["sasa_binding_site"] = self.calculate_sasa_binding_site()
        features["orientation"] = self.calculate_orientation()

        # contacts
        contacts = self.calculate_contacts()
        features["contacts_total"] = contacts["total"]
        features["contacts_hydrophobic"] = contacts["hydrophobic"]
        features["contacts_polar"] = contacts["polar"]

        # hydrogen bonds
        hbonds = self.calculate_hydrogen_bonds()
        features["hydrogen_bonds"] = hbonds["total"]

        # distance features
        dist_feats = self.calculate_distance_features()
        features.update(dist_feats)

        # info
        features["n_protein_atoms"] = (
            len(self.protein_atoms) if BIOTITE_AVAILABLE else len(self.protein_sel)
        )
        features["n_ligand_atoms"] = (
            len(self.ligand_atoms) if BIOTITE_AVAILABLE else len(self.ligand_sel)
        )
        features["n_binding_residues"] = len(self.binding_residues)

        return features


def analyze_trajectory(
    topology: str,
    trajectory: str,
    chain_id: Optional[str] = None,
    ligand_resname: Optional[str] = None,
    binding_cutoff: float = 10.0,
    start: Optional[int] = None,
    stop: Optional[int] = None,
    step: int = 1,
    n_jobs: int = 1,
) -> Dict[str, List]:
    """
    Analyze trajectory and calculate features for each frame.

    Args:
        topology: Path to topology file (PDB, PSF, etc.)
        trajectory: Path to trajectory file (xtc, dcd, etc.)
        chain_id: Chain ID for protein
        ligand_resname: Ligand residue name
        binding_cutoff: Distance cutoff for binding site
        start: Start frame
        stop: Stop frame
        step: Step size
        n_jobs: Number of parallel jobs (-1 for all CPUs)

    Returns:
        Dictionary with feature names as keys and lists of values
    """
    if not MDA_AVAILABLE:
        raise ImportError("MDAnalysis required for trajectory analysis")

    u = mda.Universe(topology, trajectory)

    n_frames = len(range(start or 0, stop or len(u.trajectory), step))
    feature_names = [
        "closeness",
        "sasa",
        "sasa_binding_site",
        "orientation",
        "contacts_total",
        "contacts_hydrophobic",
        "contacts_polar",
        "hydrogen_bonds",
        "com_distance",
        "cog_distance",
        "min_distance",
    ]

    results = {name: [] for name in feature_names}
    results["frame"] = []

    frames = list(range(start or 0, stop or len(u.trajectory), step))

    for i, ts in enumerate(u.trajectory[slice(start, stop, step)]):
        geom = ProteinLigandGeometry(
            pdb_file=topology,
            chain_id=chain_id,
            ligand_resname=ligand_resname,
            binding_cutoff=binding_cutoff,
        )

        if MDA_AVAILABLE:
            geom.u = u
            geom.protein_sel = u.select_atoms("protein")
            if chain_id:
                all_protein = u.select_atoms("protein")
                geom.protein_sel = all_protein[all_protein.chainIDs == chain_id]

            geom.ligand_sel = u.select_atoms(
                f"resname {ligand_resname}"
                if ligand_resname
                else "not protein and not resname WAT and not resname HOH"
            )

            if len(geom.ligand_sel) > 0 and geom.binding_residues:
                bs_resids = [r[0] for r in geom.binding_residues]
                geom._current_bs_sel = u.select_atoms(
                    f"resid {' '.join(map(str, bs_resids))}"
                )

        features = geom.calculate_all_features()

        results["frame"].append(i)
        for name in feature_names:
            results[name].append(features.get(name, 0))

        if (i + 1) % 100 == 0:
            print(f"Processed {i + 1}/{n_frames} frames")

    print(f"Completed: {n_frames} frames analyzed")
    return results
