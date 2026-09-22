import pandas as pd
import subprocess
from pathlib import Path

CSV_INPUT = "Top400_combinatorial_variants.csv"
FASTA_OUTPUT = "Top400_sequences.fasta"
TOOL_SCRIPT = "multiple_prediction_wrapper_export_Top400.sh"
PREDICTION_OUTPUT = "seq_prediction_Top400.txt"
FINAL_CSV = "../Top400_variants_solubility.csv"

df = pd.read_csv(CSV_INPUT, sep=";")

with open(FASTA_OUTPUT, "w") as f:
    for _, row in df.iterrows():
        header = str(row["Variant"]).strip().replace(" ", "_").replace(";", "_")
        seq = str(row["Sequence"]).strip()
        f.write(f">{header}\n{seq}\n")

result = subprocess.run(
    ["bash", TOOL_SCRIPT, FASTA_OUTPUT],
    capture_output=True,
    text=True,
    timeout=600
)

if Path(PREDICTION_OUTPUT).exists():
    with open(PREDICTION_OUTPUT) as f:
        lines = f.readlines()
else:
    lines = result.stdout.split("\n")

solubilities = []
for line in lines:
    if line.startswith("SEQUENCE PREDICTIONS"):
        parts = [p.strip() for p in line.split(",")]
        if len(parts) >= 4:
            try:
                scaled_sol = float(parts[3])
                solubilities.append(scaled_sol)
            except ValueError:
                solubilities.append(None)
        else:
            solubilities.append(None)

while len(solubilities) < len(df):
    solubilities.append(None)

df["Solubility"] = solubilities[:len(df)]

df.to_csv(FINAL_CSV, sep=";", index=False)

print(f"Done. Processed {len(df)} sequences.")
print(f"Results saved to: {FINAL_CSV}")
print(f"FASTA file kept at: {FASTA_OUTPUT}")
print("\nPreview:")
print(df[["Variant", "Solubility"]].head(10))