import os
import subprocess
import csv
import pandas as pd

def repair_pdb(input_pdb, output_folder):
    repair_command = f"EvoEF2/EvoEF2 --command=RepairStructure --pdb={input_pdb}"
    subprocess.run(repair_command, shell=True)

    file_name = os.path.splitext(os.path.basename(input_pdb))[0]
    repair_file = f"{file_name}_Repair.pdb"

    repair_file_path = os.path.join(os.getcwd(), repair_file)
    output_path = os.path.join(output_folder, repair_file)
    os.rename(repair_file_path, output_path)
    return output_path


def compute_stability(input_pdb):
    stability_command = f"EvoEF2/EvoEF2 --command=ComputeStability --pdb={input_pdb}"
    result = subprocess.check_output(stability_command, shell=True, text=True)
    
    for line in result.split('\n'):
        if line.startswith("Total"):
            stability = line.split()[2]
            return stability
    return None  

def extract_id(filename):
    name = os.path.splitext(filename)[0]
    last_underscore_pos = name.rfind('_')
    
    if last_underscore_pos != -1:
        return name[last_underscore_pos + 1:]
    else:
        return name


def process_pdb_files(input_folder, output_folder, output_csv):
    os.makedirs(output_folder, exist_ok=True)
    with open(output_csv, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Variant", "Stability"])

        for file in os.listdir(input_folder):
            if file.endswith(".pdb"):
                pdb_file = os.path.join(input_folder, file)
                repair_file = repair_pdb(pdb_file, output_folder)
                stability = compute_stability(repair_file)
                file_id = extract_id(file)
                csv_writer.writerow([file_id, stability])


if __name__ == "__main__":
    input_folder = "minimized_shuffle_structures"
    output_folder = "processed_shuffle_structures"
    output_csv = "total_energies_shuffle_library.csv"
    process_pdb_files(input_folder, output_folder, output_csv)
    stability_df = pd.read_csv("total_energies_shuffle_library.csv")
    variants_df = pd.read_csv("kemp_shuffle_variants.csv",sep=";")
    stability_df["Variant"] = stability_df["Variant"].astype(str)
    variants_df["Variant"] = variants_df["Variant"].astype(str)
    merged_df = pd.merge(stability_df, variants_df, on="Variant", how="inner")
    merged_df.to_csv("kemp_shuffle_variants_stability.csv", sep=";",index=False)

