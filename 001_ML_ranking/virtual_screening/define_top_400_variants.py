import pandas as pd

seq_file = "sequences_combinatorial_variants.csv"
pred_file = "y_preds_combinatorial.tsv"
output_joined = "combinatorial_variants_predicted_performance.csv"
shuffle_file = "/mnt/bkup/Kemp_Shuffe_ML_MD/Machine_Learning/training_shuffle_dataset/feature_generation/kemp_shuffle_variants.csv"
output_filtered = "combinatorial_variants_predicted_performance_filtered_out_shuffle.csv"
output_top400 = "Top400_combinatorial_variants.csv"

hg3_r5 = "EAAQSVDQLMKARGKVYFGVATDQNRLTTGKNAAIIQADFGMVWPENSLQWDAVEPSQGNFNFAGADYVVNWAQQNGKLIGGGMLVWHSHLPSWVSSITDKETLTNVMKNHITTLMTRYKGKIRVWDVVGSAFNEDGSLRQTVFLNVIGEDYIKIAFQTARAADPNAKLYIADSNLDSASYPKTQAIVNKVKQWRAAGVPIDGIGSMTNLSAGQGAGVLQALPLLASAGTPEVSILMLDVAGASPTDYVNVVNACLNVQSCVGITSFGVADPDSWRASTTPLLFDGNFNPKPAYNAIVQDLQQ"

def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

seq_df = pd.read_csv(seq_file, sep=";")
pred_df = pd.read_csv(pred_file, sep="\t")

seq_df = seq_df.reset_index(drop=True)
pred_df = pred_df.reset_index(drop=True)

joined_df = pd.concat([seq_df, pred_df], axis=1)

if "Unnamed: 0" in joined_df.columns:
    joined_df = joined_df.drop(columns=["Unnamed: 0"])

joined_df.to_csv(output_joined, sep=";", index=False)

shuffle_df = pd.read_csv(shuffle_file, sep=";")
existing_sequences = set(shuffle_df['Sequence'].astype(str).str.strip())

filtered_df = joined_df[~joined_df['Sequence'].astype(str).str.strip().isin(existing_sequences)].copy()

filtered_df = filtered_df[filtered_df['Sequence'].apply(lambda x: hamming_distance(x, hg3_r5) >= 4)].copy()

filtered_df.to_csv(output_filtered, sep=";", index=False)

top400 = filtered_df.nlargest(400, 'pred')
top400.to_csv(output_top400, sep=";", index=False)