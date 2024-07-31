import glob
import os
import random
import tempfile
import traceback

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol

RDLogger.DisableLog("rdApp.*")

import pickle

from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.PDBIO import Select
from scipy.spatial import distance_matrix

from constants import ATOM_TYPES

def write_xyz(types, coords, msg="", fn=None, is_onehot=True):
    xyz = ""
    xyz += f"{coords.shape[0]}\n"
    xyz += msg + "\n"
    for i in range(coords.shape[0]):
        if is_onehot:
            atom_type = ATOM_TYPES[np.argmax(types[i])]
        else:
            atom_type = types[i]
        xyz += f"{atom_type}\t{coords[i][0]}\t{coords[i][1]}\t{coords[i][2]}\n"
    if fn is not None:
        with open(fn, "w") as w:
            w.writelines(xyz[:-1])
    return xyz[:-1]


def xyz_to_sdf(xyz_fn, sdf_fn):
    os.system(f"obabel {xyz_fn} -O {sdf_fn}")


def read_file(filename):
    extension = filename.split(".")[-1]
    if extension == "sdf":
        mol = Chem.SDMolSupplier(filename)[0]
    elif extension == "mol2":
        mol = Chem.MolFromMol2File(filename)
    elif extension == "pdb":
        mol = Chem.MolFromPDBFile(filename)
    elif extension == "xyz":
        filename2 = filename[:-4] + ".sdf"
        xyz_to_sdf(filename, filename2)
        mol = Chem.SDMolSupplier(filename2)[0]
    else:
        # print("Wrong file format...")
        return
    if mol is None:
        # print("No mol from file...")
        return
    return mol


def get_properties(mol, key="ligand", pdb_id=None):
    key = key
    pdb_id = pdb_id

    properties = {}

    # 1. Molecular weight
    mw = Descriptors.MolWt(mol)

    # 2. TPSA
    tpsa = Descriptors.TPSA(mol)

    # 3. LogP
    logp = Descriptors.MolLogP(mol)

    # 6. Free SASA
    sasa = calc_free_sasa(mol)

    properties.update(
        {
            "mw": np.array([mw]),
            "tpsa": np.array([tpsa]),
            "logp": np.array([logp]),
            "sasa": np.array([sasa]),
        }
    )

    return properties


def calc_free_sasa(mol):
    """
    Code from
    https://sunhwan.github.io/blog/2021/02/04/RDKit-Protein-Ligand-SASA.html
    """
    from rdkit.Chem import rdFreeSASA

    # compute ligand SASA
    mol_h = Chem.AddHs(mol, addCoords=True)

    # Get Van der Waals radii (angstrom)
    ptable = Chem.GetPeriodicTable()
    radii = [ptable.GetRvdw(atom.GetAtomicNum()) for atom in mol_h.GetAtoms()]

    # Compute solvent accessible surface area
    sasa = rdFreeSASA.CalcSASA(mol_h, radii)
    return sasa
