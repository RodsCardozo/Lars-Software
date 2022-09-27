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

def area(vertices):
    import numpy as np
    a1 = vertices[0]
    a2 = vertices[1]
    a = a1 - a2
    e = np.linalg.norm(a)
    A = ((3)**(1/2)/4)*e**2
    return A


nu = 50
vertices, faces = icosphere.icosphere(nu)
center = []
for i in range(0, len(faces), 1):
    A = vertices[faces[i][0]]
    B = vertices[faces[i][1]]
    C = vertices[faces[i][2]]
    x = A[0] + B[0] + C[0]
    y = A[1] + B[1] + C[1]
    z = A[2] + B[2] + C[2]
    center.append([x/3, y/3, z/3])
print(faces)
print(vertices)

scale = 1
vertices = vertices*scale

faces = faces
fig = go.Figure()

fig.add_trace(mesh_plot(vertices, faces))
fig.add_trace(wireframe_plot(vertices, faces));

fig.update_layout(title_text='Icosphere', height=600, width=600)
fig.show()

