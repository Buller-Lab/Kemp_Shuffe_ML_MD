import pandas as pd

def convert_to_three_letter_code(mutations):
    three_letter_code = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
        'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
        'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
        'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }

    point_mutations = mutations.split(',')

    result = []
    for mutation in point_mutations:
        residue_from, position, residue_to = mutation[0], mutation[1:-1], mutation[-1]
        three_letter_from = three_letter_code[residue_from]
        three_letter_to = three_letter_code[residue_to]
        result.append(("A", f"{three_letter_from}-{position}-{three_letter_to}"))

    return(result)

df = pd.read_csv('kemp_shuffle_variants.csv', sep=";")

df['PDB_mutations'] = df['Mutations'].apply(convert_to_three_letter_code)

df.to_csv('kemp_shuffle_variants_pdb_mut.csv', index=False, sep=";")
