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
#from extract_energies import extract_energy
from prism_pruner.conformer_ensemble import ConformerEnsemble
from prism_pruner.pruner import prune
from pebble import ProcessExpired, ProcessPool
import numpy as np
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
EH_TO_KCAL = 627.5096080305927
### Parallel Processing TODO: Implement this on the loop
def process_parallel(function, arguments, num_workers=6, timeout=30):
    res = []
    with ProcessPool(max_workers=num_workers) as pool:
        future = pool.map(
            function,
            [id for id in arguments],
            timeout=timeout,
        )
        iterator = future.result()

        i = 0
        while True:
            if i % 100 == 0:
                print("Finished {i} iterations")
            try:
                result = next(iterator)
                # print(result)
                res.append(result)
            except StopIteration:
                print("Stop iteration")
                break
            except TimeoutError:
                print("Timeout error")
                res.append(None)
            except ProcessExpired as error:
                print("%s. Exit code: %d" % (error, error.exitcode))
                res.append(None)
            except Exception as error:
                print(arguments[i])
                print("function raised %s" % error)
                res.append(None)
                # print(error.traceback)  # Python's traceback of remote process
            i += 1
    return res
### BELOW STOLEN FROM TMC_EMBED!!
def embed_tmc_complex(smiles, xtb_path, coord_order=None, geometry=None, cpus=4):
    m = Chem.MolFromSmiles(smiles, sanitize=False)
    m = reorder_atoms_to_map(m)
    mhap = rdmolops.DativeBondsToHaptic(m)
    coord_dir = get_hapt_coord_dir(mhap, coord_order)
    if not coord_order: #choose random coord_order
        coord_order = list(coord_dir.keys())
        #print(coord_order)
    charge = Chem.GetFormalCharge(m)
    constraints, angle_constraints, new_mol = setup_constrain_dict_pos(mhap, coord_order, coord_dir)
    for i,atom in enumerate(new_mol.GetAtoms()):
        atom.SetAtomMapNum(i+1)
    sep_embed = embed_smiles_far(Chem.MolToSmiles(new_mol))
    Chem.MolToXYZFile(sep_embed, filename="sep_embed.xyz")
    os.mkdir("ff_embed")
    shutil.move("sep_embed.xyz", "ff_embed") # Move instead of copy
    os.chdir("ff_embed")
    write_constrain_file(constraints, angle_constraints, 0.0025)
    shutil.copy("constrain.inp", "../")


    try:
        output = run_cmd(f"{xtb_path} sep_embed.xyz --opt crude --input constrain.inp --gfnff --parallel {cpus} --chrg {charge}".format())
        #print(output)
        shutil.copy("xtbopt.xyz", "../sep_embed_ff.xyz")
        os.chdir("../")
        shutil.rmtree("ff_embed")
    except:
        print("stopped at ff embed part")
        os.chdir("../")
        shutil.rmtree("ff_embed")
        return -1, np.nan
    os.mkdir("gfn2_constrain")
    shutil.move("sep_embed_ff.xyz", "gfn2_constrain")
    os.chdir("gfn2_constrain")
    write_constrain_file(constraints, angle_constraints, 0.2)
    shutil.copy("constrain.inp", "../constrain2.inp")
    try:
        output = run_cmd(f"{xtb_path} sep_embed_ff.xyz --opt crude --input constrain.inp --parallel {cpus} --gfnff --chrg {charge}".format())
        #print(output)
        shutil.copy("xtbopt.xyz", "../sep_embed_final.xyz")
        os.chdir("../")
        shutil.rmtree("gfn2_constrain")
        mol = assign_coordinates2mol(sep_embed, "sep_embed_final.xyz", constraints)
        output = run_cmd(f"{xtb_path} sep_embed_final.xyz --opt verytight --iterations 1000 --input constrain2.inp --gfn2 --parallel {cpus} --chrg {charge}".format())

        with open("out.out", 'w') as _file:
            _file.writelines(output)
            #lines = _file.readlines()
        with open("out.out", 'r') as _file:
            lines = _file.readlines()
        for line in lines:
            #print(line)
            if "TOTAL ENERGY" in line:
                e = float(line.split()[3])
                #print(line)
        mol = assign_coordinates2mol(sep_embed, "xtbopt.xyz", constraints)
    except:
        print("failed at gfn2 constrain part")
        os.chdir("../")
        shutil.rmtree("gfn2_constrain")
        return -1, np.nan

    return mol, e

def get_tmc_mol(atom_mapped_smiles, xtb_path, coord_order=None, geometry=None, N_tries=30, cpus=4):
    """
    Do N_tries embeddings with final gfn2 optimizations and keep lowest energy
    """
    pos_mols = []
    pos_es = []
    for _ in range(N_tries):
        mol, e = embed_tmc_complex(atom_mapped_smiles, xtb_path, coord_order=coord_order, geometry=geometry, cpus=cpus)
        pos_mols.append(mol)
        pos_es.append(e)
       # print(f"energy = {e}")

    if np.all(pos_mols) == -1:
        print("embedding failed")
    else:
        #print(coord_order)
        print("")

    return pos_mols, pos_es  # NOTE: returning all the conformations, so I can run r2scan SP on each



def mkxyzblock(mol, energy): 
    xyzblock = Chem.MolToXYZBlock(mol)
    parts = xyzblock.split("\n", 1)  # Split only at the first \n
    xyzblock = parts[0] + "\n" + energy + "\n" + parts[1]
    return xyzblock

def make_inp_file_xTB_opt(xyz: str):
    '''
    Takes input xyz file 
    '''
    # spec_str = f"! Native-GFN-xTB\n \n %method \n     RunTyp opt\n end\n \n * XYZ 0 1 \n" #NOTE: IF ORCA is used
    xyz_str = f"{xyz[2:]}\n*"
    inp_file = xyz_str
    return inp_file
    
def run_xTB(tmp_dir, xyz, ncpu, xtb_path):
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
parser.add_argument("--smiles", default=None, type=str, help="Path to csv file containing smiles of complexes for calculation")
parser.add_argument("--smiles_col", default = "Complex", type=str, help="Column name containing the smiles of complexes for calculation")
parser.add_argument("--xtb_path", default="/home/henryteahan/opt/xtb-6.7.0/xtb-dist/bin/xtb", type=str, help="Full path to xTB binaries")
parser.add_argument("--charge", default=0, type=int, help="Central metal ion charge in complex")
parser.add_argument("--xyz_files", type=str, default=None, help="input XYZ files - one or more")
parser.add_argument("--ncpu", default=4, help="Number of processes")
parser.add_argument("--N_tries", type=int, default = 20, help = "Number of embedding attempts")
args = parser.parse_args()
xtb_path = args.xtb_path

print(f"Input arguments {args} \n \n \n ------------------ \n")
if args.xyz_files != None:
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
            possible_coord_orders = get_possible_coord_permutations(smile)
            stereo_confs = []
            energies = []
            for coord_order in possible_coord_orders:
                mol, es = get_tmc_mol(smile, xtb_path, coord_order, N_tries=args.N_tries, cpus = args.ncpu)
                stereo_confs.extend(mol)
                energies.extend(es)
        
        name = f"{path.stem}"
        for i in range(len(stereo_confs)):
            try:
                xyz_data = mkxyzblock(stereo_confs[i], str(energies[i])) # Making original xyz-files
            except:
                print("Error", stereo_confs)
            filename = f"{name}_embed_{i}.xyz"
            out_path = Path.cwd() / filename
            out_path.write_text(xyz_data)
        
            # TODO: Update xyz file handling in this 
            # RUN xTB SP 
            #out = run_xTB(tmp_dir, out_path, ncpu=args.ncpu) ### RUNS FROM THE EXECUTED FOLDER LOCATION
            # Get xTB SP energy and return the best geometry in the working directory
            #df = extract_energy(out)    
            #print(df)       
elif args.smiles != None:
    smiles = pd.read_csv(args.smiles)
    smiles = smiles[args.smiles_col]
    path = Path(args.smiles)
    tmp_dir = Path.cwd() / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(path, tmp_dir)
    smile_index = 0
    for row, smile in enumerate(smiles): #TODO: Make this into nested parallel loop.
        with pushd(tmp_dir):
            possible_coord_orders = get_possible_coord_permutations(smile)
            stereo_confs = []
            energies = []
            for coord_order in possible_coord_orders:
                mol, es = get_tmc_mol(smile, xtb_path, coord_order, N_tries=args.N_tries, cpus = args.ncpu)
                stereo_confs.extend(mol)
                energies.extend(es)
        name = f"{path.stem}"
        all_xyz_confs = ""
        for i in range(len(stereo_confs)):

            xyz_data = mkxyzblock(stereo_confs[i], str(energies[i])) # Making original xyz-files
            print(xyz_data)
            all_xyz_confs += xyz_data + "\n"
        
# Remove empty lines
        all_xyz_confs = "\n".join(line for line in all_xyz_confs.splitlines() if line.strip())
        
        with pushd(tmp_dir):
            filename = f"{name}_embed_{row}.xyz"
            out_path = Path.cwd() / filename
            out_path.write_text(all_xyz_confs)
        #TODO: rewrite to not have to write and then reread. 
        ensemble = ConformerEnsemble.from_xyz(out_path, read_energies=True)
        pruned, mask = prune(
        ensemble.coords,
        ensemble.atoms,

        rot_corr_rmsd_pruning=True,
        debugfunction=print,
        )
        atoms = ensemble.atoms
        coords = ensemble.coords[mask]
        energies = ensemble.energies[mask]
        unique_E = energies*EH_TO_KCAL
        mask_more = unique_E - unique_E.min() < 3 ### UNIQUE CONFORMERS WITHN 3 kcal / mol of the lowest!!!
        filename = f"{name}_id_{smile_index}.xyzs"
        out_path = Path.cwd() / filename
        ConformerEnsemble(coords[mask_more],atoms, energies[mask_more]).to_xyz(out_path)
        smile_index += 1
# Clean up tmp folder
shutil.rmtree(tmp_dir)
