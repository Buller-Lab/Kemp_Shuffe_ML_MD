import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import product

df = pd.read_csv("kemp_shuffle_variants_solubility.csv", sep=";")
FIO_cutoff = 1.3
df["Positive"] = df["FIO_HG3.17"] >= FIO_cutoff
y_true = df["Positive"].astype(int).values
sol = df["Solubility"].values
fio = df["FIO_HG3.17"].values

best_lower = 0.0
best_upper = 1.0
best_mcc = -1.0
best_cm = None
min_fn = float('inf')

candidates = np.unique(np.round(sol, 4))
candidates = np.sort(candidates)

for i in range(len(candidates)):
    for j in range(i, len(candidates)):
        lower = candidates[i]
        upper = candidates[j]
        y_pred = ((sol >= lower) & (sol <= upper)).astype(int)
        
        tp = np.sum((y_true == 1) & (y_pred == 1))
        tn = np.sum((y_true == 0) & (y_pred == 0))
        fp = np.sum((y_true == 0) & (y_pred == 1))
        fn = np.sum((y_true == 1) & (y_pred == 0))
        
        if fn < min_fn:
            min_fn = fn
            best_mcc = -1.0
            best_lower = lower
            best_upper = upper
            best_cm = [[tn, fp], [fn, tp]]
        
        if fn == min_fn:
            mcc_num = (tp * tn - fp * fn)
            mcc_den = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) + 1e-10)
            mcc = mcc_num / mcc_den
            if mcc > best_mcc:
                best_mcc = mcc
                best_lower = lower
                best_upper = upper
                best_cm = [[tn, fp], [fn, tp]]

print(f"Optimal Solubility Range: [{best_lower:.4f}, {best_upper:.4f}]")
print(f"False Negatives: {best_cm[1][0]}")
print(f"MCC: {best_mcc:.4f}")

cm = np.array(best_cm)
labels = ["Negative", "Positive"]
fig, ax = plt.subplots(figsize=(6, 5))
im = ax.imshow(cm, interpolation='nearest', cmap='Blues')
ax.figure.colorbar(im, ax=ax)
ax.set(xticks=np.arange(2), yticks=np.arange(2),
       xticklabels=labels, yticklabels=labels,
       title=f'Confusion Matrix\nSolubility Range [{best_lower:.4f}, {best_upper:.4f}]\nMCC = {best_mcc:.4f}',
       ylabel='True Label', xlabel='Predicted Label')

for i, j in product(range(2), range(2)):
    ax.text(j, i, f"{cm[i,j]}", ha="center", va="center",
            color="white" if cm[i,j] > cm.max()/2 else "black")

plt.tight_layout()
plt.savefig("confusion_matrix_range.png", dpi=200)

plt.figure(figsize=(8, 6))
plt.scatter(fio, sol, c=df["Positive"].map({True: 'green', False: 'red'}), alpha=0.7, s=40)

plt.axhline(y=best_lower, color='darkgray', linestyle='--', linewidth=2, label=f'Lower cutoff = {best_lower:.4f}')
plt.axhline(y=best_upper, color='gray', linestyle='-', linewidth=2, label=f'Upper cutoff = {best_upper:.4f}')

plt.xlabel("FIO_HG3.17")
plt.ylabel("Solubility")
plt.title(f"Solubility vs FIO_HG3.17\nGreen = Positive (FIO ≥ {FIO_cutoff}) | Red = Negative")
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("scatter_solubility_vs_fio.png", dpi=200)
plt.show()

print("Plots saved:")
print("- confusion_matrix_range.png")
print("- scatter_solubility_vs_fio.png")