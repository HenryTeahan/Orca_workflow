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
import threading
import time
### INITIATE DB
#DB_PATH = Path.home() / "Projects/Orca_workflow/Automated/db/jobs.db"

def is_job_running(slurm_job_id):
    result = subprocess.run(
            ["squeue", "-j", str(slurm_job_id), "-h"],
            capture_output=True, text=True)
    print("Job check", result)
    return bool(result.stdout.strip())
def monitor_orca_jobs(conn, cur): #TODO: MAKE THIS ERROR HANDLING WORK!!!!!
    cur.execute("""
                SELECT job_id, inp_file, slurm_job_id
                FROM jobs
                WHERE status="orca_submitted"
                """)
    rows = cur.fetchall()
    if not rows:
        print("No submitted orca tasks yet, take a nap!!")
        return
    for job_id, inp_file, slurm_job_id in rows:
        inp_path = Path(inp_file)
        out_file = Path("Orca_run") / inp_path.name.replace(".inp", ".out")
        err_file = Path("Orca_run") / inp_path.name.replace(".inp", ".err")
        if is_job_running(slurm_job_id):
            continue
    
        if err_file.exists() and err_file.stat().st_size > 0:
            print(f"Orca failed :-( {inp_file}")
            cur.execute("""
                        UPDATE jobs
                        SET status="error", error="Orca failed"
                        WHERE job_id=?
                        """, (job_id,))
            conn.commit()
            continue
        try: #If no error and job complete
            df = extract_energy(out_file)
            gibbs_free = df['Final Gibbs free energy']
            cur.execute("""
                        UPDATE jobs
                        SET best_energy = ?
                        WHERE inp_file = ?
                        """, (gibbs_free, inp_file))
        except Exception as e:
            print("Error in parsing", e)

def submit_orca(inp_file, job_dir, ncpus = 8, mem = 16000):
    inp_stem = Path(inp_file).stem
    slurm_filename = job_dir / f"submit_{inp_stem}.slurm"
    sbatch_script = f"""#!/bin/bash
#SBATCH --job-name={inp_file}
#SBATCH --nodes=1
#SBATCH --ntasks={ncpus}
#SBATCH --mem {mem}
#SBATCH --time=02:00:00
#SBATCH --output={job_dir}/{inp_file}.out
#SBATCH --error={job_dir}/{inp_file}.err
#SBATCH --partition=kemi1
#Move to scratch dir for calculations
module load mpi/openmpi-x86_64
cp {inp_file} /scratch/$SLURM_JOB_ID/
cd /scratch/$SLURM_JOB_ID/

/groups/kemi/hteahan/opt/orca_6_1_0_linux_x86-64_shared_openmpi418/orca {inp_file}

mv *.out {job_dir}
mv *.err {job_dir}
""".format(mem=mem, ncpus=ncpus, job_dir=job_dir, inp_file=inp_file)
    
    slurm_filename.write_text(sbatch_script)
    result = subprocess.run(["sbatch", str(slurm_filename)],
                            capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Failed to submit orca job: {result.stderr}")
    job_id = result.stdout.strip().split()[-1]
    return job_id # Track slurm task by id

def monitoring_thread(DB_PATH, poll_interval=30):
    thread_conn = sqlite3.connect(DB_PATH)
    thread_cur = thread_conn.cursor()
    while True:
        monitor_orca_jobs(thread_conn, thread_cur)  # the function from before
        time.sleep(poll_interval)



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
    slurm_job_id TEXT,
    error TEXT,
    best_energy REAL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP)""")
    conn.commit()
    ### END

    print("MADE DATABASE") 
    t = threading.Thread(target=monitoring_thread, args=(DB_PATH,))
    t.daemon = False
    t.start()

    outdir = Path(Path.cwd() / "Orca_run")
    outdir.mkdir(parents=True, exist_ok=True)

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
                        INSERT INTO jobs
                        (ligand_id, mol_id, ligand_smiles, xyz_file, xtb_energy, status, job_dir)
                        VALUES (?, ?, ?, ?, ?, ?, ?)
                        """, (conformers["ligand_ID"], conformers['mol_ID'], smile, conformers["XYZ"], json.dumps(conformers["xTB_energies"]), "complete", str(embed_dir)))
            conn.commit()
        except Exception as e:
            print("Error in embedding:", e)
            cur.execute("""
            INSERT INTO jobs
            (ligand_id, mol_id, ligand_smiles, status)
            VALUES (?, ?, ?, ?)
            """, (lig_ID, mol_ID, smile, "ERROR in embedding"))
            conn.commit()
            break
        try:    
            xyz = Path(embed_dir)/conformers["XYZ"]
            multiplicity = np.abs(charge) + 1 # Use simplest states -> If S = 0; multiplicity=1, If S = 1/2; multiplicity 2... 
            inp_files = xyzs_to_inp(xyz, charge, multiplicity, ncpus = 8, mem = 16000)        
            for inp_file in inp_files:
                cur.execute("""
                    INSERT INTO jobs (ligand_id, mol_id, ligand_smiles, xyz_file, inp_file, status, job_dir)
                    VALUES (?,?,?,?,?,?,?)
                    """, (lig_ID, mol_ID, smile, str(xyz), str(inp_file), "inputs_created", str(outdir)))
                conn.commit()
        except Exception as e:
            print("Error in INP File generation:", e)
            cur.execute("""
            UPDATE jobs
            SET status='input_gen_failed'
            WHERE xyz_file=? and ligand_id=?
            """, (str(xyz), lig_ID))

        for inp_file in inp_files:
            inp_path = Path(inp_file)
            print(inp_path)
            try:
                job_id = submit_orca(inp_path, job_dir=outdir, ncpus=8, mem=16000)
                print(f"Submitted ORCA job {inp_file} as Slurm ID {job_id}")

                # Update DB row for this specific inp_file
                cur.execute("""
                    UPDATE jobs
                    SET status='orca_submitted', slurm_job_id=?
                    WHERE inp_file=?
                """, (job_id, str(inp_file)))
                conn.commit()
            except Exception as e:
                print(f"Error submitting {inp_file}: {e}")
                cur.execute("""
                    UPDATE jobs
                    SET status='error_in_orca_submit', error=?
                    WHERE inp_file=?
                """, (str(e), str(inp_file)))
                conn.commit()
    t.join()



            
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

