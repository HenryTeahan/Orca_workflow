### IMPORTS
from rdkit import Chem
import pandas as pd
import argparse
from rdkit.Chem.MolStandardize import rdMolStandardize
import re

pt = Chem.GetPeriodicTable()

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
    Joins the ligand and the right side of the molecule.
    Preserves ligand order during failure.
    """
    mols = []
    smiles = []

   # right_side = Chem.MolFromSmiles(right_side)
    print("Joining right side", right_side, "with ligands")
    for ligand in ligands:
        try:
            # Default failure
            mols.append(-1)
            smiles.append(-1)
            #
            Chem.Kekulize(ligand)
            # Combine ligand and isomer
            comb = Chem.CombineMols(right_side, ligand)
            edcombo = Chem.EditableMol(comb)
            b_idx = []
            rem_idx = []
            m_idx = -1
            for atom in comb.GetAtoms():
                if atom.GetAtomicNum() == 28: #TODO: Assumes use of Ni metal 
                    m_idx = atom.GetIdx()
                    continue # Skip central metal ion

                if atom.HasProp("binding_site"):
                    b_idx.append(atom.GetIdx()) # Append if binding site    
                    valence = atom.GetTotalValence()
                    norm_valence = pt.GetDefaultValence(atom.GetAtomicNum())
                    remaining = norm_valence - valence
                    rem_idx.append(remaining)
                    
            if m_idx == -1 or not b_idx:
                continue
                
            for rem, val in zip(rem_idx, b_idx):
                if rem >= 1:
                    edcombo.AddBond(val, m_idx, order=Chem.BondType.SINGLE)
                else:
                    edcombo.AddBond(val, m_idx, order=Chem.BondType.DATIVE)
            
            
            f = edcombo.GetMol()

            for atom in f.GetAtoms():
                if atom.HasProp("binding_site"):
                    atom.SetFormalCharge(0)
            s = Chem.MolToSmiles(f, canonical=True)
            smiles[-1] = s
            mols[-1] = f
        except:
            #print("Error")
            continue
    return mols, smiles
def parse_rejoined_mols(ligands, right_side):
    if isinstance(right_side, list):
        right_side = [Chem.MolFromSmiles(r) for r in right_side]
    else:
        right_side = Chem.MolFromSmiles(right_side)
    print("isomer in mol", right_side)
    output, smiles_output = rejoin_mols(ligands, right_side)
    #print(output, smiles_output)
    #output = [t for t in output if type(t) != int]
    smiles_output = []
    clean_mols = []
    for i, m in enumerate(output):
        try:
            for atom in m.GetAtoms():
                num_radicals = atom.GetNumRadicalElectrons() # REMOVE RANDOM RADICALS LEFT OVER
                if num_radicals:
                    atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
                    atom.SetNumRadicalElectrons(0)
            m.UpdatePropertyCache()
            
            Chem.SanitizeMol(m)
            m = Chem.AddHs(m)
            for i,a in enumerate(m.GetAtoms()):
                a.SetAtomMapNum(i+1)

            smiles_output.append(Chem.MolToSmiles(m))
            clean_mols.append(m)
        except:
            smiles_output.append(-1)
            clean_mols.append(-1)
            #print("FAILED", i)
            continue
    successful_smiles = [s for s in smiles_output if type(s) != int] 
    print(f"Succesfully parsed {len(successful_smiles)} mols ({len(successful_smiles)/len(ligands)*100:.2f}%)")
    return smiles_output, clean_mols
#smiles, mols = parse_rejoined_mols(unique_mols, central_metal)
if __name__ == "__main__":
    # BEGIN INPUT PARSING
    parser = argparse.ArgumentParser()
    parser.add_argument("INPUT", type=str, help="INPUT SMILES CSV FILE")
    parser.add_argument("--smiles_column", type=str, default = "smiles_huckel_DFT_xyz")
    parser.add_argument("--isomer_pair", nargs="+", help="Right hand side of the molecule - isomer pairs and metal ion (Input string)",
                        default = ["O=C1O[Ni]N(C1C1=CC=CC=C1)C1=CC=CC=C1","O=C1O[Ni]C(N1C1=CC=CC=C1)C1=CC=CC=C1"])
    parser.add_argument("--outname", type=str, default="smiles_out.csv", help="Name of output file")
    # By default, uses isomer1 and isomer2s
    args = parser.parse_args()
    ### END
    
    # READ INPUT SMILES:
    combined_frame = pd.read_csv(args.INPUT)
    sms = combined_frame[args.smiles_column].to_list() #We use the DFT_xyz smiles
    mols = [Chem.MolFromSmiles(s) for s in sms]
    print(f"Parsed {len(mols)} molecules")
    # END
    
    # LOAD MetalDisconnectorOptions - Set TRANSITION_METALS_NUM
    params = rdMolStandardize.MetalDisconnectorOptions()
    params.splitAromaticC = True
    params.splitGrignards = True
    params.adjustCharges = False
    # Load MetalDisconnector
    mdis = rdMolStandardize.MetalDisconnector(params)
    TRANSITION_METALS_NUM = [21,22,23,24,25,26,27,57,28,29,30,39,40,41,42,43,44,45,46,47,48,71,72,73,74,75,76,77,78,79,80]
    # END

    # Extract bidentate ligands
    ligands = []
    mol_id = []
    for i, mol in enumerate(mols):
        frag = identify_bidentate(mol, TRANSITION_METALS_NUM, mdis)
        if frag != -1: # identify_bidentate returns -1 if it cannot find a bidentate ligand or if it errors.
            if isinstance(frag,list):
                for f in frag:
                    f.SetIntProp("ID", i)
                    ligands.append(f)
                    mol_id.append(i)
            else:
                frag.SetIntProp("ID", i) # Ligand identifier
                ligands.append(frag)
                mol_id.append(i)
    print("Found this many bidentate ligands!!", len(ligands))
    # END
    
    # Parse extracted bidentate ligands
    ligand_highlight = []
    ligands_flattened = []

    for n, mol in enumerate(ligands):
        #if isinstance(mol, list): # Capture cases where complex has >1 bidentate ligands
        #    for m in mol:
        #        highlight_atoms = highlight_bindsite(m) # Extracts binding site set by identify_bidentate
        #        ligand_highlight.append(highlight_atoms)
        #        ligands_flattened.append(m)
        #else:
        highlight_atoms = highlight_bindsite(mol)
        ligand_highlight.append(highlight_atoms)
        ligands_flattened.append(mol)
    # END
    
    # Remove duplicate ligands
    ligands, bindsites = remove_duplicate_mols(ligands_flattened, ligand_highlight)
    print("Number of unique mols", len(ligands))
    #print(ligands)
    # END

    # TODO: Initiate the df here!
    df = pd.DataFrame()
    df['ligand_ID'] = [l.GetProp("ID") for l in ligands]
    df['ligands'] = [Chem.MolToSmiles(l) for l in ligands]
    df['bind_site'] = bindsites
     
    for id, isomer in enumerate(args.isomer_pair):
        # Bind unique ligands to the input isomer.
        print("Isomer", isomer)
        new_smiles, new_mols = parse_rejoined_mols(ligands, isomer)
        # END
        df[f'smiles_{id}'] = new_smiles
    # END
    # NOTE: smiles and ligands_smiles now contain ligand, smile pairs.
    df = df[df['smiles_0'] != "-1"]
    df = df[df['smiles_0'] != -1]
    df.to_csv(args.outname)
    print("Saved smiles_out.csv with the complexes for screening")
    print(df.head(n=10))
