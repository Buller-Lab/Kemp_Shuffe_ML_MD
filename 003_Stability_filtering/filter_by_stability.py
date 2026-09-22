import pandas as pd

df = pd.read_csv("total_energies_Top400_filtered_solubility.csv", sep=";")
lower_cutoff = -1687.8600
upper_cutoff = -1637.5800

filtered = df[(df["Stability"] >= lower_cutoff) & (df["Stability"] <= upper_cutoff)]

filtered.to_csv("Top400_filtered_solubility_stability.csv", sep=";", index=False)

print(f"Total variants in Top400: {len(df)}")
print(f"Variants after stability filter [{lower_cutoff:.4f}, {upper_cutoff:.4f}]: {len(filtered)}")
print(f"Filtered file saved as: Top400_filtered_solubility_stability.csv")
