import pandas as pd
import plotly.graph_objects as go
from io import StringIO


df = pd.read_csv('LD_info_chr_all.txt' ,sep='\t')
pivot = df.pivot(index='SNP_A', columns='SNP_B', values='R2')

# Prepare labels and matrix
x_labels = pivot.columns.astype(str)
y_labels = pivot.index.astype(str)
z_matrix = pivot.values

# Available colorscales
colorscales = ['Viridis', 'Plasma', 'Cividis', 'Inferno', 'Magma']

# Create heatmap
heatmap = go.Heatmap(
    z=z_matrix,
    x=x_labels,
    y=y_labels,
    colorscale=colorscales[0],
    colorbar=dict(title='R2')
)

fig = go.Figure(data=[heatmap])

# Dropdown buttons
buttons = []
for scale in colorscales:
    buttons.append(dict(
        method='restyle',
        label=scale,
        args=[{'colorscale': scale}, [0]]
    ))

fig.update_layout(
    title='Example LD Heatmap (R2)',
    xaxis_title='SNP_B',
    yaxis_title='SNP_A',
    updatemenus=[dict(
        active=0,
        buttons=buttons,
        x=1.15,
        y=1.1,
        xanchor='left',
        yanchor='top',
        direction='down',
        pad={'r': 10, 't': 10}
    )],
    margin=dict(l=80, r=200, t=80, b=80)
)

fig.show()
