import pandas as pd
from scipy import stats

def statistics(out,sep,X,y,cutoff):
    df = pd.read_csv(f"{out}.csv", index_col=None, sep=sep,low_memory=False)
    if cutoff != None:
        curated_df = df[df.Angle <= cutoff]
    else:
        curated_df = df    
    statistics_df= curated_df.groupby([X])[y].median().to_frame('Median')
    statistics_df["Mean"] = curated_df.groupby([X])[y].mean().to_frame('Mean')
    statistics_df["Std"] = curated_df.groupby([X])[y].std().to_frame('Std')
    statistics_df["Var"] = curated_df.groupby([X])[y].var().to_frame('Var')
    statistics_df["Min"] = curated_df.groupby([X])[y].min().to_frame('Min')   
    statistics_df["Max"] = curated_df.groupby([X])[y].max().to_frame('Max')  
    groups = curated_df.groupby([X])[y]

    for i in groups.groups:
        statistics_df[f"levene-{i}"]= "n.d."
    statistics_df2 = statistics_df.copy()
 
    for i in groups.groups:
        for j in range(len(groups.groups)):
            if j != i:
                group1 = groups.get_group(i)
                group2 = groups.get_group(list(groups.groups)[j])
                res = stats.levene(group1, group2)
                statistics_df2.at[i, f"levene-{list(groups.groups)[j]}"] = res.pvalue
                if res.pvalue > 0.05:
                    res2 = stats.ttest_ind(group1, group2, equal_var=True)
                    statistics_df2.at[i, f"ttest-{list(groups.groups)[j]}"] = res2.pvalue
                else:
                    res2 = stats.ttest_ind(group1, group2, equal_var=False)
                    statistics_df2.at[i, f"ttest-{list(groups.groups)[j]}"] = res2.pvalue
    statistics_df2.to_csv(f"{out}_statistics.csv", sep=sep)