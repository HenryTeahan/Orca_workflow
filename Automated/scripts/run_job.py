import sqlite3
from pathlib import Path
import json
from embed import embed
from xyz2inp import xyzs_to_inp
from submit import curate_inp
import argparse
import pandas as pd
import textwrap
import subprocess
### INITIATE DB
DB_PATH = Path.home() / "Projects/Orca_workflow/Automated/db/jobs.db"

def main(args):
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
    inp_file, TEXT,
    status TEXT,
    job_dir TEXT,
    error TEXT,
    best_energy REAL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP)""")
    conn.commit()
    ### END

    print("MADE DATABASE") 
    scratch = Path(args.SCRATCH_DIR)
    scratch.mkdir(parents=True, exist_ok=True)
    
    embed_dir = scratch / "embed"
    embed_dir.mkdir(exist_ok=True)
    orca_dir = scratch / "Orca"
    orca_dir.mkdir(exist_ok=True)
    
    smiles_df = pd.read_csv(Path(args.SMILES), index_col=0)
    smiles = smiles_df[args.SMILES_COL]
    mol_ID = smiles_df[args.ID_COL] # THis is index column
    ligand_ID = smiles_df.index.to_list()
    xtb_path = args.xTB_path
    for smile, mol_ID, lig_ID in zip(smiles, mol_ID, ligand_ID):
        conformers = embed(smile, mol_ID, lig_ID, xtb_path, embed_dir)
        
        for conf in conformers:
            cur.execute("""
                        INSERT OR IGNORE INTO jobs
                        (ligand_id, mol_id, ligand_smiles, xyz_file, xtb_energy, status, job_dir)
                        VALUES (?, ?, ?, ?, ?, ?, ?)
                        """, (conf["ligand_ID"], conf['mol_ID'], smile, conf["XYZ"], json.dumps(conf["xTB_energies"]), "complete", str(embed_dir)))
        conn.commit()

        xyzs = Path(embed_dir)/conf["XYZ"]
        inp_files = xyzs_to_inp(xyzs) ###TODO: make this parse args such as charge, method, etc. Charge and multiplicity needs to be dynamically generated.
       ########################################################################################################### 
        cur.execute("""
            UPDATE jobs
            SET status="inputs_created", inp_file=?
            WHERE xyz_file=?
                    """, (json.dumps([str(f) for f in inp_files]), xyzs))
        conn.commit()

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--SMILES", type=str, help="File containing smiles", required = True)
    parser.add_argument("--SMILES_COL", type=str, help="Column in SMILES containing the smiles", required = True)
    parser.add_argument("--ID_COL", type=str, default="ID", help="ID column in SMILES making tracking of ligands easier")
    #parser.add_argument("ISOMER", type) #TODO: Make throughput.py run in here
    parser.add_argument("--JOB_ID", type=str, help="NAME OF JOB", required = True)
    parser.add_argument("--SCRATCH_DIR", type=str, help="HPC Scratch Dir", required = True)
    parser.add_argument("--xTB_path", type=str, default="~/opt/xtb-6.7.1/xtb-dist/bin/xtb", help="Path to xTB binary")
    args = parser.parse_args()    
    
    xtb_path = Path(args.xTB_path).expanduser()

    main(args)

