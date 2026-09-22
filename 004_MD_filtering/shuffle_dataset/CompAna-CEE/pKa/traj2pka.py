import os
import csv
import random
import pandas as pd
import mdtraj as mdt
from concurrent.futures import ProcessPoolExecutor

def process_trajectory(traj_path, variant, residue, chain):
    max_retries=2
    for attempt in range(max_retries + 1):
        try:
            traj = mdt.load(traj_path)
            break  # Break the loop if loading is successful
        except Exception as e:
            print(f"Error loading trajectory {traj_path} (Attempt {attempt + 1}/{max_retries + 1}): {e}")
            if attempt == max_retries:
                return []  # Return an empty list if max retries reached

    data = []

    for i in range(traj.n_frames):
        content = [variant, traj_path, chain]
        frame = traj[i]
        unique_temp_key = ''.join((random.choice('abcdxyzpqr') for _ in range(10)))
        frame.save_pdb(f"{unique_temp_key}_tmp.pdb")
        os.system(f"propka3 {unique_temp_key}_tmp.pdb")

        with open(f"{unique_temp_key}_tmp.pka", "r") as fi:
            for ln in fi:
                if ln.startswith(f"   {residue} {chain}"):
                    content.append(ln.split("    ")[0])
                    content.append(ln.split("    ")[1])

        data.append(content)
        os.system(f"rm {unique_temp_key}_tmp.pdb")
        os.system(f"rm {unique_temp_key}_tmp.pka")

    return data

def traj2pka(path, pre, ext, variant, residue, chain, out, sep):
    max_workers = 30

    # Filter for .h5 files
    prefixed = [filename.path for filename in os.scandir(path) if filename.name.startswith(f"{pre}") and filename.name.endswith(f"{ext}")]
    data_to_write = []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_trajectory, traj_path, variant, residue, chain) for traj_path in prefixed]

        with open(f"{out}.csv", "a", newline="") as csv_file:
            csv_writer = csv.writer(csv_file, delimiter=sep)

            for future in futures:
                data_to_write.extend(future.result())

                # Write the new data directly to the CSV file row by row
                if not os.path.exists(f"{out}.csv"):
                    # Write header only if the file doesn't exist
                    csv_writer.writerow(future.result()[0].keys())

                for row in future.result():
                    csv_writer.writerow(row)

    print(f"Processing variant {variant} completed.")