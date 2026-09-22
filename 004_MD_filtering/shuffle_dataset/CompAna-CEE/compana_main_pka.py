import pandas as pd
from pKa.traj2pka import traj2pka

import warnings
warnings.filterwarnings("ignore")

### USER INPUT
#########################################################################################

# definition of trajectories
pre ="production_"

# definition of trajectory file format (compatible with MDTraj) 
ext = "h5"

# definition of separator in output csv
sep=";"

# definition of output prefix
out = "HG3_pka"

# definition of path to variant trajectories
traj_path = "/home/stcg/MD_screening/MD_"

# list of variant names, should equal the folder in which all trajectories which should be analyzed are located
variant_list = ["HG3.R5","HG3.17","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150","151","152","153","154","155","156","157","158","159","160","161","162","163","164","165","166","167","168","169","170","171","172","173","174","175","176","177","178","179","180","181","182","183","184","185","186","187","188","189","190","191","192","193","194","195","196","197","198"]

# list of residue strings in the following format: "RES 123"
residue_list = ["ASP 127"]

# list of chain identifier strings
chain_list = ["A"]

### WORKFLOW DEFINITION
#########################################################################################

# initiate csv and generate header   
header = pd.DataFrame([["Variant","Trajectory","Chain","Residue","pKa"]])
header.to_csv(f"{out}.csv", sep=f"{sep}", header = None, index=None)

# definition of all loops to interate over variants, residues, chains and pH values
for variant in variant_list:
    path = traj_path + variant
    for residue in residue_list:
        for chain in chain_list:
            traj2pka(path, pre, ext, variant, residue, chain, out, sep)
