# Kemp_Shuffe_ML_MD

git clone https://github.com/Buller-Lab/Kemp_Shuffe_ML_MD.git

cd ML_prescreening/

Note: variant99 in input kemp_shuffle_variants.csv was removed because of "Z" in sequence which could not be featurized (artefact from sequencing)

At first, features are derived from kemp sequences via protparam and protlearn python libraries

python 001_protlearn_features_from_fasta_training.py

