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
### END

### Functions

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
    
def run_xTB(workdir, filename):
    ''' Attempts running unconstrained xTB calculations '''

    workdir = Path(workdir)
    inp = Path(filename)


    print(f"Running ORCA xTB Opt: {inp.name}")
    print(f"Current directory: {os.getcwd()}")
    print(f"Input file: {inp}")
    
    xTB_opt = workdir / "OUTPUT"
    xTB_opt = Path(xTB_opt)
    xTB_opt.mkdir(exist_ok=True)

    out = inp.with_suffix(".out").name
    print(out)
    try:
        with open(out, "w") as f: #TODO: Make this in a try-except loop. xTB can & will error silently - not producing xtbopt.xyz
            subprocess.run(
                [xtb_path, str(inp), "--ohess", "verytight", "--parallel", "4", "--gfn2"],
                stdout= f,
                stderr= subprocess.STDOUT,
                check= True,
                cwd = workdir
            )
    except subprocess.CalledProcessError as e:
        print(f"xTB failed for {inp.name}")


    opt_file = workdir / f"xtbopt.xyz"
    log_file = workdir / f"xtbopt.log"
    out_file = workdir / out
    if opt_file.exists() and log_file.exists() and out_file.exists():
        new_path = xTB_opt / f"opt-{inp.stem}.xyz"
        opt_file.rename(new_path)
        trj_path = xTB_opt / f"opt-{inp.stem}.log"
        log_file.rename(trj_path)
        out_path = xTB_opt / out
        out_file.rename(out_path)
    else:
        print("COULDNT FIND OPTIMIZED FILE OR LOG FILE")

    tempnames = [
        "xtbopt.log",
        "xtbopt.xyz",
        "xtbrestart",
        "xtbtopo.mol",
        "xtbopt",
        "constrain.inp",
        "constrain2.inp",
        "charges",
        "out.out",
        "wbo",
        "sep_embed_final.xyz",
        ".xtboptok",
        "hessian",
        "g98.out",
        "vibspectrum"
    ]

    for name in tempnames:
        try:
            os.remove(workdir / name)
        except FileNotFoundError:
            pass
        try:
            os.remove(xTB_opt / name)
        except FileNotFoundError:
            pass

def get_energies(workdir):
    subprocess.run(["/home/henryteahan/bin/get_energies"],
                    check=True,
                    cwd=workdir) #TODO: incorporate this script in python            
### END



parser = argparse.ArgumentParser()

parser.add_argument("--xtb_path", default="/home/henryteahan/opt/xtb-6.7.0/xtb-dist/bin/xtb", type=str, help="Full path to xTB binaries")
parser.add_argument("--charge", default=0, type=int, help="Central metal ion charge in complex")
parser.add_argument("xyz_files", nargs="+", help="input XYZ files - one or more")
parser.add_argument("--unconstr", default="False", help="Attempt to run unconstrained optimization. May break desired geometries.")
args = parser.parse_args()

xtb_path = args.xtb_path

print("Input arguments", args)
for xyz in args.xyz_files:
    path = Path(xyz)
    print(f"Embedding {path.name}")

    smiles = extract_TMC_smiles(os.path.join(os.getcwd(), str(path)), charge = 0)
    target = f"Ni"
    print(smiles[smiles.find(target)-10:smiles.find(target)+10])
    pattern = r"(Ni)@.*?(:)"
    smiles = re.sub(pattern, r"\1+2\2", smiles)
    print(smiles[smiles.find(target)-10:smiles.find(target)+10])

    possible_coord_orders = get_possible_coord_permutations(smiles)
    stereo_confs = []
    for coord_order in possible_coord_orders:
        mol = get_tmc_mol(smiles, xtb_path, coord_order) # NOTE: We are now returning all mols and min_Es for R2scan SP
        stereo_confs.append(mol)
    # Cleanup
    dir = Path(os.getcwd())
    tempnames = [
        "xtbopt.log",
        "xtbopt.xyz",
        "xtbrestart",
        "xtbtopo.mol",
        "xtbopt",
        "constrain.inp",
        "constrain2.inp",
        "charges",
        "out.out",
        "wbo",
        "sep_embed_final.xyz",
        ".xtboptok"
    ]

    for name in tempnames:
        try:
            os.remove(dir / name)
        except FileNotFoundError:
            pass





    # Making folder for the input files!

    name = f"{path.stem}"

    # TODO: Fix the embedding workflow and the handling of input and output files.
    
    base = Path("INPUT")
    base.mkdir(parents=True, exist_ok=True)  
    for i in range(len(stereo_confs)):
        #input_files = base / "INP_files"
        xyz_files = base / "XYZ_files"
        #base.mkdir(input_files, exist_ok=True)
        xyz_files.mkdir(exist_ok=True)
        #min_E_index = np.argmin(energies)
        print(stereo_confs[i])
        print(stereo_confs)        
        try:
            xyz = mkxyzblock(stereo_confs[i]) # Making original xyz-files
        except:
            print("Error", stereo_confs)
        filename = f"{name}_embed_{i}.xyz"
        path_xyz = xyz_files / filename
        path_xyz.write_text(xyz)
    
    #inp = make_inp_file_xTB_opt(xyz)
    #
    #filename = f"{name}.inp"
    #path_inp = input_files / filename
    #path_inp.write_text(inp)

    if args.unconstr == "True":
        # Make output directory
        output = Path("OUTPUT")
        output.mkdir(exist_ok=True)
        run_xTB(os.getcwd(), path_xyz) ### RUNS FROM THE EXECUTED FOLDER LOCATION
        get_energies(base)

