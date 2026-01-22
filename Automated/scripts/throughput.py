### IMPORTS
from rdkit import Chem
import pandas as pd
import argparse
from rdkit.Chem.MolStandardize import rdMolStandardize
import re
### END
def identify_bidentate(mol, TRANSITION_METALS_NUM, mdis):
    """
    Identifies bidentate ligands by fragmenting the complex, and counting number of bonds to the central metal ion -> If a ligand has two bonds: bidentate
    """
    map = 0
    neighbours = 0
    
    # Find the central metal ion; and enumerate the atoms using rdkit atom props -> Consistent through splitting.
    for atom in mol.GetAtoms():
        atom.SetIntProp("map", map)
        map += 1
        if atom.GetAtomicNum() in TRANSITION_METALS_NUM:
            neighbours = atom.GetNeighbors()
        else:
            continue

    if neighbours == 0:
        return -1 #TODO: Make this into a valueerror handling
    
    neighbours = [n.GetIdx() for n in neighbours]
    
    # Begin fragmenting
    out = mdis.Disconnect(mol)
    frags = []
    for f in Chem.GetMolFrags(out, asMols=True, fragsMolAtomMapping=[]):
        frags += [f]

    bidentate_idx = []
    for i, m in enumerate(frags):
        # Count how many bonds each fragment has to the central metal ion
        count = 0
        for atom in m.GetAtoms():
            if int(atom.GetProp("map")) in neighbours:
                count += 1
                atom.SetIntProp("binding_site", 1)
        if count == 2:
            bidentate_idx.append(i)
    if len(bidentate_idx) == 1:
        return frags[bidentate_idx[0]]
    elif len(bidentate_idx) > 1:
        return [frags[b_idx] for b_idx in bidentate_idx]
    else:
        return -1

def highlight_bindsite(mol):
    """
    Finds the binding sites and returns the atom indexes -> Used to rejoin the molecules after!        
    """
    highlight_atoms = []
    for atom in mol.GetAtoms():
        if atom.HasProp("binding_site"):
            if int(atom.GetProp("binding_site")) != 0:  
                highlight_atoms.append(atom.GetIdx())
    return highlight_atoms

def remove_duplicate_mols(ligands, ligand_highlight):
    """
    Finds and removes duplicate mols if they exist. Uses canonical smiles and CanonicalRankAtoms to create identifier hash (smi_bind).
    """
    seen = set()
    bindsites = []
    unique_mols = []
    
    for high, mol in zip(ligand_highlight, ligands):
        ranks = Chem.CanonicalRankAtoms(mol)
        site_ranks = sorted(ranks[idx] for idx in high)
        try:
            smi = Chem.MolToSmiles(mol, canonical=True)
        except:
            continue
        smi_bind = f"{smi}{tuple(site_ranks)}"
      #  print(smi_bind, high)
        if smi_bind not in seen:
            seen.add(smi_bind)
            unique_mols.append(mol)
            bindsites.append(high)
    return unique_mols, bindsites

def rejoin_mols(ligands, right_side):
    """
    Joins the ligand and the right side of the molecule
    """

    output = []
    seen = set()
    right_side = Chem.MolFromSmiles(right_side)

    for ligand in ligands:
        # print(Chem.MolToSmiles(ligand))
        Chem.Kekulize(ligand, clearAromaticFlags=True)
        comb = Chem.CombineMols(right_side, ligand)
        edcombo = Chem.EditableMol(comb)
        b_idx = []
        m_idx = -1

        for atom in comb.GetAtoms():
            if atom.GetAtomicNum() == 28: ### Assumes use of Ni metal 
                m_idx = atom.GetIdx()
                continue
            try:   
                atom.GetProp("binding_site")
                b_idx.append(atom.GetIdx())
            except:
                continue

        if m_idx == -1:
            print("Couldn't find central metal")
            break

        for val in b_idx:
            try:
                edcombo.AddBond(m_idx, val, order=Chem.BondType.SINGLE)
            except:
                continue

        f = edcombo.GetMol()
        s = Chem.MolToSmiles(f, canonical=True)

        if s in seen:
            continue
        else:
            seen.add(s)
            output.append(f)
    return output, list(seen)

if __name__ == "__main__":
# BEGIN
    parser = argparse.ArgumentParser()
    parser.add_argument("INPUT", type=str, help="INPUT SMILES CSV FILE")
    parser.add_argument("--smiles_column", type=str, default = "smiles_huckel_DFT_xyz")
    parser.add_argument("--right_side", type=str, help="Right hand side of the molecule - amino acid and metal ion (Input string)", default = "O=C1O[Ni]N(C1C1=CC=CC=C1)C1=CC=CC=C1") # TODO: Make this accept a csv of right hand side structures.
    args = parser.parse_args()
### END

    combined_frame = pd.read_csv(args.INPUT)

    sms = combined_frame[args.smiles_column].to_list()
    mols = [Chem.MolFromSmiles(s) for s in sms]
    print(f"Parsed {len(mols)} molecules")

    params = rdMolStandardize.MetalDisconnectorOptions()
    params.splitAromaticC = True
    params.splitGrignards = True
    params.adjustCharges = False
# Load MetalDisconnector
    mdis = rdMolStandardize.MetalDisconnector(params)

    TRANSITION_METALS_NUM = [21,22,23,24,25,26,27,57,28,29,30,39,40,41,42,43,44,45,46,47,48,71,72,73,74,75,76,77,78,79,80]
    
# Extract bidentate ligands
    ligands = []
    mol_id = []
    for i, mol in enumerate(mols):
        frag = identify_bidentate(mol, TRANSITION_METALS_NUM, mdis)
        if frag != -1:
            ligands.append(frag)
            mol_id.append(i)

    print("Found this many bidentate ligands!!", len(ligands))
    
    ligand_highlight = []
    ligands_flattened = []

    for n, mol in enumerate(ligands):
        if isinstance(mol, list):
            for m in mol:
                highlight_atoms = highlight_bindsite(m)
                ligand_highlight.append(highlight_atoms)
                ligands_flattened.append(m)
        else:
            highlight_atoms = highlight_bindsite(mol)
            ligand_highlight.append(highlight_atoms)
            ligands_flattened.append(mol)
    ligands = ligands_flattened

    ligands, bindsites = remove_duplicate_mols(ligands, ligand_highlight)

    print("Number of unique mols", len(ligands))
    #Chem.Draw.MolsToGridImage(unique_mols, highlightAtomLists=bindsites, molsPerRow=5, subImgSize=(500,500), maxMols=50)
    

    new_mols, new_smiles = rejoin_mols(ligands, args.right_side)

    smiles = []
    
    for s in new_smiles:
        try:
            m = Chem.AddHs(Chem.MolFromSmiles(s))
        except:
            continue
        for i, a in enumerate(m.GetAtoms()):
            a.SetAtomMapNum(i+1)
        s = Chem.MolToSmiles(m)
        target = f"Ni" #TODO: ADD multiple metal recognition here dependent on the right side of the complex
        pattern = r"(Ni)@.*?(:)"
        s = re.sub(pattern, r"\1+2\2", s)
        smiles.append(s)
    print(len(smiles))
    df = pd.DataFrame()
    df['Complex'] = smiles
    df.to_csv("smiles_out.csv")
    print("Saved smiles_out.csv with the complexes for screening")
