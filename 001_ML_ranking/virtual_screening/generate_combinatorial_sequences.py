import pandas as pd

def count_digits_in_string(input_string):
    return sum(1 for char in input_string if char.isnumeric())

variants = pd.read_csv("combinatorial_variants.csv", index_col=None)

sequence = "EAAQSVDQLIKARGKVYFGVATDQNRLTTGKNAAIIQADFGMVWPENSMKWDATEPSQGNFNFAGADYLVNWAQQNGKLIGGGMLVWHSQLPSWVSSITDKNTLTNVMKNHITTLMTRYKGKIRAWDVVGEAFNEDGSLRQTVFLNVIGEDYIPIAFQTARAADPNAKLYIMDYNLDSASYPKTQAIVNRVKQWRAAGVPIDGIGSQTHLSAGQGAGVLQALPLLASAGTPEVSILMLDVAGASPTDYVNVVNACLNVQSCVGITVFGVADPDSWRASTTPLLFDGNFNPKPAYNAIVQDLQQ"

variants["Sequence"] = ""

for idx, row in variants.iterrows():
    mutation_list = row["Mutations"].split(",")
    sequence_new = sequence
    for mutation in mutation_list:
        len_number = count_digits_in_string(mutation)
        wt = mutation[0]
        res_numb = mutation[1:1 + len_number]
        mut = mutation[1 + len_number]
        pos = int(res_numb) - 1
        if wt == sequence_new[pos]:
            sequence_new = sequence_new[:pos] + mut + sequence_new[pos + 1:]
    variants.at[idx, "Sequence"] = sequence_new

variants.to_csv("sequences_combinatorial_variants.csv", index=False, sep=";")