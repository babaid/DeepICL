from rdkit import Chem

DATA_DIR = "" ## directory to PDBbind v2020 general set 
SAVE_DIR = "" ## directory where processed data will be saved

ATOM_TYPES = ["C", "N", "O", "F", "P", "S", "Cl", "Br"]
AA_TYPES = [
    "GLY", "ALA", "VAL", "LEU", "ILE",
    "PHE", "PRO", "MET", "TRP", "SER",
    "THR", "TYR", "CYS", "ARG", "HIS",
    "LYS", "ASN", "ASP", "GLN", "GLU"
]
INTERACTION_TYPES = ["pipi", "anion", "cation", "hbd", "hba", "hydro"]
HYDROPHOBICS = ["F", "CL", "BR", "I"]
HBOND_DONOR_SMARTS = ["[!#6;!H0]"]
HBOND_ACCEPTOR_SMARTS = [
    "[$([!#6;+0]);!$([F,Cl,Br,I]);!$([o,s,nX3]);!$([Nv5,Pv5,Sv4,Sv6])]"
]
SALT_ANION_SMARTS = ["[O;$([OH0-,OH][CX3](=[OX1])),$([OX1]=[CX3]([OH0-,OH]))]"]
SALT_CATION_SMARTS = [
    "[N;$([NX3H2,NX4H3+;!$(NC=[!#6]);!$(NC#[!#6])][#6])]",
    "[#7;$([NH2X3][CH0X3](=[NH2X3+,NHX2+0])[NHX3]),$([NH2X3+,NHX2+0]=[CH0X3]([NH2X3])[NHX3])]",
    "[#7;$([$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]1:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[#6X3]1),$([$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]1:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3]:[#6X3H]1)]",
]
DEGREES = [0, 1, 2, 3, 4]
HYBRIDIZATIONS = [
    Chem.rdchem.HybridizationType.S,
    Chem.rdchem.HybridizationType.SP,
    Chem.rdchem.HybridizationType.SP2,
    Chem.rdchem.HybridizationType.SP3,
    Chem.rdchem.HybridizationType.SP3D,
    Chem.rdchem.HybridizationType.SP3D2,
]
FORMALCHARGES = [-2, -1, 0, 1, 2, 3]
SYMBOL_TO_MASS = {
    "C": 12.0,
    "N": 14.0,
    "O": 16.0,
    "F": 19.0,
    "P": 31.0,
    "S": 32.1,
    "Cl": 35.5,
    "Br": 80.0,
}
MAX_ATOM_NUM = 50
MAX_ADD_ATOM_NUM = 30
SEED = 0


ATOM_NUM = [6, 7, 8, 9, 15, 16, 17, 35]
ATOM_MASS = [12.0, 14.0, 16.0, 19.0, 31.0, 32.1, 35.5, 79.9]

INTERACTION_TYPES = ["pipi", "anion", "cation", "hbd", "hba", "hydro"]

NUM_LIGAND_ATOM_TYPES = len(ATOM_TYPES) + 1
NUM_POCKET_ATOM_TYPES = 51
NUM_INTERACTION_TYPES = len(INTERACTION_TYPES) + 1

PERIODIC_TABLE = """                                                       
H,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,HE
LI,BE,1,1,1,1,1,1,1,1,1,1,B,C,N,O,F,NE
NA,MG,1,1,1,1,1,1,1,1,1,1,AL,SI,P,S,CL,AR
K,CA,SC,TI,V,CR,MN,FE,CO,NI,CU,ZN,GA,GE,AS,SE,BR,KR
RB,SR,Y,ZR,NB,MO,TC,RU,RH,PD,AG,CD,IN,SN,SB,TE,I,XE
CS,BA,LU,HF,TA,W,RE,OS,IR,PT,AU,HG,TL,PB,BI,PO,AT,RN
"""