import sqlite3
from pathlib import Path
import json
from embed import embed
from xyz2inp import xyzs_to_inp
from extract_energies import extract_energy
import argparse
import pandas as pd
import textwrap
import subprocess
import numpy as np
### INITIATE DB
#DB_PATH = Path.home() / "Projects/Orca_workflow/Automated/db/jobs.db"



def submit_orca(inp_file, job_dir, ncpus = 8, mem = 16000): #TODO: FIX THE ROUTING HERE - write in scratch but copy out of!
    inp_stem = Path(inp_file).stem
    slurm_filename = f"submit_{inp_stem}.slurm"
    sbatch_script = f"""#!/bin/bash
#SBATCH --job-name={inp_file}
#SBATCH --cpus_per_task={ncpus}
#SBATCH --mem {mem}
#SBATCH --time=10-00:00:00
#SBATCH --output={job_dir}/{inp_file}.out
#SBATCH --error={job_dir}/{inp_file}.err
#SBATCH --partition=kemi1
#Move to scratch dir for calculations

cd /scratch/$SLURM_JOB_ID/

/groups/kemi/hteahan/opt/orca_6_1_0_linux_x86-64_shared_openmpi418/orca {inp_file}

mv *.out {job_dir}
mv *.err {job_dir}
""".format(mem=mem, ncpus=ncpus, job_dir=job_dir, partition=partition, inp_file=inp_file)
    
    slurm_filename.write_text(sbatch_script)
    result = subprocess.run(["sbatch", str(slurm_filename)],
                            capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Failed to submit orca job: {result.stderr}")
    job_id = result.stdout.strip().split()[-1]
    return job_id # Track slurm task by id


def main(args):
    DB_PATH = Path(args.DB_PATH)
    DB_PATH.parent.mkdir(exist_ok=True, parents=True)
    

    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()

    cur.execute("""
CREATE TABLE IF NOT EXISTS jobs (
    job_id INTEGER PRIMARY KEY,
    ligand_id INTEGER,
    ligand_smiles TEXT,
    mol_ID INTEGER,
    complex_smiles TEXT,
    xyz_file TEXT,
    xtb_energy TEXT,
    inp_file TEXT,
    status TEXT,
    job_dir TEXT,
    error TEXT,
    best_energy REAL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP)""")
    conn.commit()
    ### END

    print("MADE DATABASE") 
    outdir = Path(Path.cwd() / "Output")

    scratch = Path(args.SCRATCH_DIR)
    scratch.mkdir(parents=True, exist_ok=True)
    
    embed_dir = scratch / "embed"
    embed_dir.mkdir(exist_ok=True)
    orca_dir = scratch / "Orca"
    orca_dir.mkdir(exist_ok=True)
    
    smiles_df = pd.read_csv(Path(args.SMILES), index_col=0)
    smiles = smiles_df[args.SMILES_COL]
    charges = smiles_df[args.CHARGE_COL]
    mol_ID = smiles_df[args.ID_COL] # THis is index column
    ligand_ID = smiles_df.index.to_list()
    xtb_path = args.xTB_path
    
    for smile, charge, mol_ID, lig_ID in zip(smiles, charges, mol_ID, ligand_ID): #TODO: Parallelize
        try:
            conformers = embed(smile, mol_ID, lig_ID, xtb_path, embed_dir)
            print(conformers) 
            cur.execute("""
                        INSERT OR IGNORE INTO jobs
                        (ligand_id, mol_id, ligand_smiles, xyz_file, xtb_energy, status, job_dir)
                        VALUES (?, ?, ?, ?, ?, ?, ?)
                        """, (conformers["ligand_ID"], conformers['mol_ID'], smile, conformers["XYZ"], json.dumps(conformers["xTB_energies"]), "complete", str(embed_dir)))
            conn.commit()
        except Exception as e:
            print("Error in embedding:", e)
            cur.execute("""
            INSERT OR IGNORE INTO jobs
            (ligand_id, mol_id, ligand_smiles, status)
            VALUES (?, ?, ?, ?)
            """, (lig_ID, mol_ID, smile, "ERROR in embedding"))
            conn.commit()
            break
        try:    
            xyz = Path(embed_dir)/conformers["XYZ"]
            multiplicity = np.abs(charge) + 1 # Use simplest states -> If S = 0; multiplicity=1, If S = 1/2; multiplicity 2... 
            inp_files = xyzs_to_inp(xyz, charge, multiplicity, ncpus = 8, mem = 16000)        
            cur.execute("""
                UPDATE jobs
                SET status="inputs_created", inp_file=?
                WHERE xyz_file=?
                """, (json.dumps([str(f) for f in inp_files]), conformers["XYZ"]))
            conn.commit()
        except Exception as e:
            print("Error in INP File generation:", e)
            cur.execute("""
            INSERT OR IGNORE INTO jobs
            (ligand_id, mol_id, ligand_smiles, status)
            VALUES (?, ?, ?, ?)
            """, (lig_ID, mol_ID, smile, "ERROR in INP FILE Generation"))
        for inp_file in inp_files:
            try:
                job_id = submit_orca(inp_file, job_dir=outdir, ncpus=8, mem=16000)

                print(f"Submitted ORCA job {inp_file} as Slurm ID {job_id}")
                cur.execute("""
                            UPDATE jobs
                            SET status='orca_submitted', slurm_job_id=?
                            WHERE inp_file LIKE ?
                            """, (job_id, f"%{inp_file.name}%"))
                conn.commit()
            except Exception as e:
                print(f"Error submitting {inp_file}: {e}")
                cur.execute("""
                            UPDATE jobs
                            SET status='error_in_orca_submit', slurm_job_id=?
                            WHERE inp_file LIKE ?
                            """, (job_id, f"%{inp_file.name}%"))

            
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--SMILES", type=str, help="File containing smiles", required = True)
    parser.add_argument("--SMILES_COL", type=str, help="Column in SMILES containing the smiles", required = True)
    parser.add_argument("--DB_PATH", type=str, help="Full path to DB", required = True)
    parser.add_argument("--CHARGE_COL", type=str, help="Column in SMILES containing the charges", default = "charge_0")
    parser.add_argument("--ID_COL", type=str, default="ID", help="ID column in SMILES making tracking of ligands easier")
    #parser.add_argument("ISOMER", type) #TODO: Make throughput.py run in here
    parser.add_argument("--JOB_ID", type=str, help="NAME OF JOB", required = True)
    parser.add_argument("--SCRATCH_DIR", type=str, help="HPC Scratch Dir", required = True)
    parser.add_argument("--xTB_path", type=str, default="~/opt/xtb-6.7.1/xtb-dist/bin/xtb", help="Path to xTB binary")
    args = parser.parse_args()    
    
    xtb_path = Path(args.xTB_path).expanduser()

    main(args)

