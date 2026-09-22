import csv
import itertools
from typing import List

def generate_variants(fixed: List[str], variable: List[str], min_mut: int, max_mut: int, output_csv: str = "combinatorial_variants.csv"):
    def get_pos(mut: str) -> int:
        for i, c in enumerate(mut):
            if c.isdigit():
                j = i
                while j < len(mut) and mut[j].isdigit():
                    j += 1
                return int(mut[i:j])
        raise ValueError(f"Invalid mutation: {mut}")
    
    fixed_set = set(fixed)
    fixed_pos = {get_pos(m) for m in fixed}
    var_list = [m for m in variable if get_pos(m) not in fixed_pos]
    
    variants = []
    variant_id = 1
    
    if min_mut == 0:
        mutations_str = ",".join(fixed) if fixed else "WT"
        variants.append((variant_id, mutations_str))
        variant_id += 1
    
    for k in range(max(min_mut, 1 if fixed else 0), max_mut + 1):
        for combo in itertools.combinations(var_list, k):
            combo_pos = {get_pos(m) for m in combo}
            if len(combo_pos) == len(combo):
                all_muts = sorted(fixed + list(combo))
                variants.append((variant_id, ",".join(all_muts)))
                variant_id += 1
    
    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Variant", "Mutations"])
        for vid, muts in variants:
            writer.writerow([vid, muts])

if __name__ == "__main__":
    fixed = ["M49L","K50Q","L69V","A125V","E131S","M172A","Y174S","Q207M","H209N","V266S"]
    variable = ["I10M", "Q37K", "N47E","T54V","G82A","M84C","S89N","Q90H","Q90F","N102E","T105I","T142N","P154K","R190K","T208M","F267M","W275A","R276F","T279S","D300N"]
    generate_variants(fixed, variable, min_mut=3, max_mut=7)
