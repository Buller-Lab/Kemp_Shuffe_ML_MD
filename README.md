# Kemp_Shuffe_ML_MD

git clone https://github.com/Buller-Lab/Kemp_Shuffe_ML_MD.git

conda create environment from ml.yml

cd Kemp_Shuffe_ML_MD

## Machine Learning

### Training based on shuffle dataset

Note: variant99 in input kemp_shuffle_variants.csv was removed because of "Z" in sequence which could not be featurized (artefact from sequencing)

At first, features are derived from kemp sequences via protparam and protlearn python libraries

cd Machine_Learning/training_shuffle_dataset/feature_generation

python protlearn_features_from_fasta_training.py

Then we select a limited number of features for training to avoid overfitting:

cd ../feature_selection

python training_feature_selection

Up to a maximum of 100 features, model performance is estimated via 10-fold cross validation. These values are stored in the feature selection monitor csv file and the top feature set (Top 24 features) was selected based on R2 and RMSE opting to minimize the number of features.

Based on this number of features, the model was trained and stored as kemp_model_24.pickle

cd ..

python training.py

### Virtual screening of combinatorial library

Navigate to the virtual screening folder:

cd ../virtual_screening_combinatorial

The next step is to create a combinatorial library. For this, we investigated the variants from the shuffle dataset that displayed performances higher than HG3.17 and HG3.R5 to derive a promising set of mutations. Based on this analysis, we fixed 10 mutations that are introduced in all variants and calculated the full combinatorial library (>120 000 variants)

Here you can see the defined mutations:

    fixed = ["M49L","K50Q","L69V","A125V","E131S","M172A","Y174S","Q207M","H209N","V266S"]
    variable = ["I10M", "Q37K", "N47E","T54V","G82A","M84C","S89N","Q90H","Q90F","N102E","T105I","T142N","P154K","R190K","T208M","F267M","W275A","R276F","T279S","D300N"]

python generate_combinatorial_library.py

Afterward, full sequences are mutated based on HG3 wildtype sequence:

python generate_combinatorial_sequences.py

Based on these sequences, features are generated (this is expected to take 1-2 hours depending on the hardware):

cd feature_generation

python protlearn_features_from_fasta_prediction.py

The virtual screening can be performed based on these predictions:

cd ..
python prediction.py

Afterwards we filter out shuffle variants, ensure that the variants contain minimally 4 mutations compared to HG3.R5 to provide sufficient novelty/effect and rank the top 400.

## Solubility filtering 

To increase probability of successful expression, solubility scores of the Shuffle dataset are investigated. 

At first solubility scores are calculated using Protein-sol:

cd ..
cd ../Solubility_filtering/protein-sol-sequence-prediction-software
python predict_solubilities_shuffle_dataset.py


As the solubility scores did not show a linear relationship but seem to converge around the top variants, upper and lower cutoffs we are defined to allow data-driven filtering of the combinatorial variants. To avoid filtering out very good variants, no false negatives were tolerated hereby (even if this causes a high number of false positives):

cd ..
python define_solubility_cutoffs.py

Similarly, solublities are calculated for the top 400 combinatorial variants:

cd protein-sol-sequence-prediction-software
python predict_solubilities_Top400.py 
cd ..

## Stability filtering

To avoid severe energetic destabilization or unfavorable rigidification, we also investigate the optimal total energy range based on the Shuffle Dataset using EvoEF2 (implemented in a previously reported script).

cd ..
cd Stability_filtering

We first convert the mutations from single letter code to an openmm compatible format

python convert_single_letter_code_shuffle.py

Afterwards the structures are created via PDBFixer and minimized using OpenMM

python create_minimized_structures_shuffle.py

Those structures are used as input to calculate the total energies via EvoEF2:

python calculate_total_energies_shuffle.py

Analogously to the solubility cutoff definition, we define a range of tolerated total energies:

python define_stability_cutoff.py

This can be applied to the solubility filtered Top400 combinatorial variants after determining their structure and calculating their total energy:

python convert_single_letter_code_Top400_filtered_solubility.py
python create_minimized_structures_Top400_filtered_solubility.py
python calculate_total_energies_Top400_filtered_solubility.py
python filter_by_stability.py

## MD-based virtual screening

Applying the solubility and stability cutoffs increases the chance of soluble expression and proper folding, however this does not necessarily account for high catalytic efficiency. Based on a previous publication and we investigated relations between kcat, Km and kcat/Km values of HG3, HG3.17 and HG3.R1-R5 and several MD simulation derived parameters. 

For MD simulation, we utilized a published high-throughput MD protocol relying on 20 separate replicates of 1 ns simulation time. (cite paper)

### MD simulations of HG3, HG3.17, HG3.R1-R5

These simulations were conducted in the beginning for HG3, HG3.17 and HG3.R1-R5 (all scripts, input & parameter files are provided here on GitHub and all resulting trajectories are deposited on Zenodo (provide ZENODO DOI) ):

cd ../MD_simulations/HG3_17_R1-R5
python HT_MD_HG3_X.py # with X representing the respective variant


### CompAna-CEE library (separate Github Repository, like with Homologic)

To allow efficient comparative analysis of different parameters for multiple variants each simulated in several replicates, we developed CompAna-CEE (Comparative Analysis of Conformational Enzyme Ensembles).

Describe here what Compana does, which functions it has and which parameters can be analyzed.

### CompAna-CEE of HG3, HG3.17 and HG3.R1-R5 trajectories

demonstrate how to use with HG3, HG3.17, HG3.R1-R5 including tables and figures

# for distances
python dist

# for angles

# for pKa values

### MD simulations of Shuffle Library

trajectories were not uploaded to zenodo due to large file sizes but are available upon request, however Compana-CEE results are provided on Zenodo

### CompAna-CEE of Shuffle Library variant trajectories

### Cutoff definition for all different categories

### MD simulations of Top 400 ML predicted combinatorial variants (solubility & stability filtered --> X variants)

trajectories were not uploaded to zenodo due to large file sizes but are available upon request, however Compana-CEE results are provided on Zenodo
cd ../virtual_screening
python simulate_variants.py

### CompAna-CEE of Top 400 ML predicted combinatorial variants (solubility & stability filtered --> X variants)


## MD-based Filtering

Based on all defined cutoffs, the final library was successibely filtered

python dist filter 1
python dist filter 2
python angle filter 1
python angle filter 2
python pka filter

## Final variant selection

Either phylogeny or random resulting in Top 10 --> experimental characterization

TO DO
Define CompAna & EvoEF as submodules!!!