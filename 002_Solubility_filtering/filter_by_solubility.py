import pandas as pd

df = pd.read_csv("Top400_variants_solubility.csv", sep=";")

lower_cutoff = 0.421
upper_cutoff = 0.443

filtered = df[(df["Solubility"] >= lower_cutoff) & (df["Solubility"] <= upper_cutoff)]

filtered.to_csv("Top400_filtered_solubility.csv", sep=";", index=False)

print(f"Total variants in Top400: {len(df)}")
print(f"Variants after solubility filter [{lower_cutoff:.4f}, {upper_cutoff:.4f}]: {len(filtered)}")
print(f"Filtered file saved as: Top400_filtered_solubility.csv")
