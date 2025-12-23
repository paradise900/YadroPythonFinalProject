from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from typing import Iterable

def smiles_valid(smiles: str) -> bool:
    return Chem.MolFromSmiles(smiles) is not None

def substructure_match(smiles: str, substructure: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    sub = Chem.MolFromSmiles(substructure)
    if not mol or not sub:
        return False
    return mol.HasSubstructMatch(sub)

def find_matches(items: Iterable[dict], substructure: str) -> list[dict]:
    sub = Chem.MolFromSmiles(substructure)
    if not sub:
        return []
    out = []
    for m in items:
        mol = Chem.MolFromSmiles(m["smiles"])
        if mol and mol.HasSubstructMatch(sub):
            out.append(m)
    return out

def tanimoto(smiles1: str, smiles2: str) -> float:
    m1 = Chem.MolFromSmiles(smiles1)
    m2 = Chem.MolFromSmiles(smiles2)
    if not m1 or not m2:
        return 0.0
    fp1 = AllChem.GetMorganFingerprintAsBitVect(m1, radius=2, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(m2, radius=2, nBits=2048)
    return float(DataStructs.TanimotoSimilarity(fp1, fp2))
