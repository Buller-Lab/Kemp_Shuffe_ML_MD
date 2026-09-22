import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import MaxNLocator

csv_file="feature_selection_monitor_kemp.csv"
png_file="feature_selection_monitor_kemp.png"
svg_file="feature_selection_monitor_kemp.svg"

df=pd.read_csv(csv_file,sep=';')
df['n_features']=range(1,len(df)+1)
df['rmse']=-df['rmse']

selected_features=24
selected_idx=df[df['n_features']==selected_features].index[0]

r2=df.loc[selected_idx,'r2']
r2_std=df.loc[selected_idx,'r2_stdv']
rmse=df.loc[selected_idx,'rmse']
rmse_std=df.loc[selected_idx,'rmse_stdv']

sns.set_style("whitegrid")

fig,axs=plt.subplots(1,2,figsize=(14,6))

r2_plot=axs[0].errorbar(df['n_features'],df['r2'],yerr=df['r2_stdv'],fmt='o-',capsize=4,linewidth=2)
for bar in r2_plot[2]:
    bar.set_alpha(0.2)

axs[0].axvline(selected_features,color='black',linestyle='--',linewidth=2)
axs[0].legend([f'Features: {selected_features}\nR²: {r2:.3f} ± {r2_std:.3f}'],loc='lower right',frameon=True)
axs[0].set_xlabel('Number of Features',fontsize=16)
axs[0].set_ylabel('R²',fontsize=16)
axs[0].tick_params(axis='both',labelsize=13)
axs[0].xaxis.set_major_locator(MaxNLocator(integer=True))

rmse_plot=axs[1].errorbar(df['n_features'],df['rmse'],yerr=df['rmse_stdv'],fmt='o-',capsize=4,linewidth=2,color='darkorange')
for bar in rmse_plot[2]:
    bar.set_alpha(0.2)

axs[1].axvline(selected_features,color='black',linestyle='--',linewidth=2)
axs[1].legend([f'Features: {selected_features}\nRMSE: {rmse:.3f} ± {rmse_std:.3f}'],loc='upper right',frameon=True)
axs[1].set_xlabel('Number of Features',fontsize=16)
axs[1].set_ylabel('RMSE',fontsize=16)
axs[1].tick_params(axis='both',labelsize=13)
axs[1].xaxis.set_major_locator(MaxNLocator(integer=True))

plt.tight_layout()
plt.savefig(png_file,dpi=300,bbox_inches='tight')
plt.savefig(svg_file,bbox_inches='tight')
plt.show()