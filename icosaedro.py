import numpy as np
import icosphere
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def mesh_plot(vertices, faces):
    gm = go.Mesh3d(x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
                   i=faces[:, 0], j=faces[:, 1], k=faces[:, 2])
    return gm


def wireframe_plot(vertices, faces):
    Xe = np.concatenate((vertices[faces, 0], np.full((faces.shape[0], 1), None)),
                        axis=1).ravel()
    Ye = np.concatenate((vertices[faces, 1], np.full((faces.shape[0], 1), None)),
                        axis=1).ravel()
    Ze = np.concatenate((vertices[faces, 2], np.full((faces.shape[0], 1), None)),
                        axis=1).ravel()

    gm = go.Scatter3d(x=Xe, y=Ye, z=Ze, mode='lines', name='',
                      line=dict(color='rgb(40,40,40)', width=1))
    return gm


nu = 15
vertices, faces = icosphere.icosphere(nu)
scale = 6371
vertices = vertices*scale
print(vertices)
print(faces)
faces = faces
fig = go.Figure()

fig.add_trace(mesh_plot(vertices, faces))
fig.add_trace(wireframe_plot(vertices, faces));

fig.update_layout(title_text='Icosphere', height=600, width=600)
fig.show()