### IMPORTS
import os
import rdkit
import argparse
from pathlib import Path
from rdkit import Chem
import pandas as pd
import subprocess, re
import shutil
from TMC_embed.xyz2mol.xyz2mol_local import *
from TMC_embed.xyz2mol.xyz2mol_local_tmc import *
from TMC_embed.utils import *
from TMC_embed.tmc_embed import *
from contextlib import contextmanager
from extract_energies import extract_energy
### END

### Functions
@contextmanager
def pushd(new_dir: Path):
    old_dir = Path.cwd()
    os.chdir(new_dir)
    try:
        yield
    finally:
        os.chdir(old_dir)


def mkxyzblock(xyzblock: str): 
    xyzblock = Chem.MolToXYZBlock(xyzblock)
    return xyzblock

def make_inp_file_xTB_opt(xyz: str):
    '''
    Takes input xyz file 
    '''
    # spec_str = f"! Native-GFN-xTB\n \n %method \n     RunTyp opt\n end\n \n * XYZ 0 1 \n" #NOTE: IF ORCA is used
    xyz_str = f"{xyz[2:]}\n*"
    inp_file = xyz_str
    return inp_file
    
def run_xTB(tmp_dir, xyz, ncpu):
    ''' Runs xTB SP calcuation on output files from embedding, chooses the lowest energy embeddding '''

    print(f"Running orca xtb sp on: {xyz}")
    
    out = tmp_dir / Path(xyz).with_suffix(".out").name
    
    print(out)
    try:
        with open(out, "w") as f: #TODO: Make this in a try-except loop. xTB can & will error silently - not producing xtbopt.xyz
            subprocess.run(
                [xtb_path, str(xyz), "--SP", "--parallel", str(ncpu), "--gfn2"],
                stdout= f,
                stderr= subprocess.STDOUT,
                check= True,
                cwd = tmp_dir
            )
    except subprocess.CalledProcessError as e:
        print(f"xtb failed for {xyz}")

    # NOTE: Need to rewrite below to only return energy and move files around dependent on this...!

    if (tmp_dir / out).exists():
        shutil.copy(out, Path.cwd() / out.name)
        return out

parser = argparse.ArgumentParser()
print("Starting parsing")
parser.add_argument("--xtb_path", default="/home/henryteahan/opt/xtb-6.7.0/xtb-dist/bin/xtb", type=str, help="Full path to xTB binaries")
parser.add_argument("--charge", default=0, type=int, help="Central metal ion charge in complex")
parser.add_argument("xyz_files", nargs="+", help="input XYZ files - one or more")
parser.add_argument("--ncpu", default=4, help="Number of OpenMP processes (cpus)")
args = parser.parse_args()

xtb_path = args.xtb_path

print(f"Input arguments {args} \n \n \n ------------------ \n")
for xyz in args.xyz_files:
    path = Path(xyz)
    print(f"Embedding {path.name}")

    # The xyz file must be in the working directory.
    if not (Path.cwd() / path.name).exists():
        shutil.copy(path, Path.cwd() / path.name)

    # Move into tmp directory and run calculations.
    tmp_dir = Path.cwd() / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(path, tmp_dir)
    # Use TMC_embed tools and extract smiles and embed molecules.
    xyz = tmp_dir /  path.name
    
    # Extract the transition metal files from the .xyz file.
    # TODO: Add the option to not use .xyz files and rather use smiles.
    with pushd(tmp_dir):
        smiles = extract_TMC_smiles(str(xyz), charge = 0)
    target = f"Ni"
    print(smiles[smiles.find(target)-10:smiles.find(target)+10])
    pattern = r"(Ni)@.*?(:)"
    smiles = re.sub(pattern, r"\1+2\2", smiles)
    print(smiles[smiles.find(target)-10:smiles.find(target)+10])

    with pushd(tmp_dir):
        possible_coord_orders = get_possible_coord_permutations(smiles)
        stereo_confs = []
        for coord_order in possible_coord_orders:
            mol = get_tmc_mol(smiles, xtb_path, coord_order, N_tries=1)
            stereo_confs.append(mol)
     
    name = f"{path.stem}"
    for i in range(len(stereo_confs)):
        try:
            xyz_data = mkxyzblock(stereo_confs[i]) # Making original xyz-files
        except:
            print("Error", stereo_confs)
        filename = f"{name}_embed_{i}.xyz"
        out_path = Path.cwd() / filename
        out_path.write_text(xyz_data)
    
    
        # RUN xTB SP 
        out = run_xTB(tmp_dir, out_path, ncpu=args.ncpu) ### RUNS FROM THE EXECUTED FOLDER LOCATION
        # Get xTB SP energy and return the best geometry in the working directory
        df = extract_energy(out)    
        print(df)
    # TODO: Make energy handling here -> run PRISM Pruner and keep unique confs within 3 kcal/mol of each other -> Use these in the next process (R2Scan)
# Clean up tmp folder
shutil.rmtree(tmp_dir)
