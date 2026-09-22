import plotly.express as px
import pandas as pd

def violine_plot(out,sep,X,y,color_list,cutoff,x_range,y_range,width,height):
    df = pd.read_csv(f"{out}.csv", index_col=None, sep=sep,low_memory=False)
    if cutoff != None:
        curated_df = df[df.pKa <= cutoff]
    else:
        curated_df = df    
    fig = px.violin(curated_df, x=X,  y=y, color=X, color_discrete_sequence=color_list, box=True, points=False,width=width, height=height)
    if x_range != None:
        fig.update_xaxes(range=x_range)
    if y_range != None:   
        fig.update_yaxes(range=y_range)        
    fig.write_image(f"{out}.svg")