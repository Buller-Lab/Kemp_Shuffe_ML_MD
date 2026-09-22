import os
import csv
import random
import pandas as pd
import mdtraj as mdt
from concurrent.futures import ProcessPoolExecutor

def process_trajectory(traj_path, variant, atom1, atom1_alt, atom_2, chain):
    max_retries=2
    for attempt in range(max_retries + 1):
        try:
            traj = mdt.load(traj_path)
            break  # Break the loop if loading is successful
        except Exception as e:
            print(f"Error loading trajectory {traj_path} (Attempt {attempt + 1}/{max_retries + 1}): {e}")
            if attempt == max_retries:
                return []  # Return an empty list if max retries reached
    topology = traj.topology
    selected_atom_1 = topology.select(atom1)
    selected_atom_2 = topology.select(atom1_alt)
    selected_atom_3 = topology.select(atom_2)
    atom_pairs = [[selected_atom_1[0],selected_atom_3[0]],[selected_atom_2[0],selected_atom_3[0]]]
    distances = mdt.compute_distances(traj,atom_pairs)
    df = pd.DataFrame(distances)
    df.columns  = ["atom1","atom1_alt"]
    df["Variant"]=variant
    df["Path"]=traj_path
    df["Chain"]=chain    
    df["Distance"] = df[["atom1","atom1_alt"]].min(axis=1)
    df=df.iloc[: , 2:]
    return df

def measure_distance(path, pre, ext, variant, atom1, atom1_alt, atom2, chain, out, sep):
    max_workers = 30

    # Filter for .h5 files
    prefixed = [filename.path for filename in os.scandir(path) if filename.name.startswith(f"{pre}") and filename.name.endswith(f"{ext}")]

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_trajectory, traj_path, variant, atom1, atom1_alt, atom2, chain) for traj_path in prefixed]

        with open(f"{out}.csv", "a", newline="") as csv_file:
            csv_writer = csv.writer(csv_file, delimiter=sep)

            for future in futures:
                new_df = future.result()
                
                # Write the new data directly to the CSV file row by row
                if not os.path.exists(f"{out}.csv"):
                    # Write header only if the file doesn't exist
                    csv_writer.writerow(new_df.columns)
                
                for _, row in new_df.iterrows():
                    csv_writer.writerow(row)

    print(f"Processing variant {variant} completed.")