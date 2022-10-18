import numpy as np
import icosphere
import numpy as np
import pandas as pd
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

def area(vertices, faces, Raio):
    import numpy as np
    a1 = vertices[faces[0][0]]*Raio
    a2 = vertices[faces[0][1]]*Raio
    a = np.array(a1) - np.array(a2)
    e = np.linalg.norm(a)
    A = ((3)**(1/2)/4)*e**2
    return A


nu = 10
vertices, faces = icosphere.icosphere(nu)
'''print((faces[0][0]))
print(vertices)
print((vertices[0][0]))
print((vertices[5][0]))
print((vertices[1][0]))'''
center = []
Raio_terra = 6371
x = []
y = []
z = []
for i in range(0, len(faces), 1):
    a = faces[i][0]
    b = faces[i][1]
    c = faces[i][2]
    A = vertices[a]*(Raio_terra)
    B = vertices[b]*(Raio_terra)
    C = vertices[c]*(Raio_terra)
    x.append((A[0] + B[0] + C[0])/3)
    y.append((A[1] + B[1] + C[1])/3)
    z.append((A[2] + B[2] + C[2])/3)
    center.append([x[i]/3, y[i]/3, z[i]/3])

import plotly.graph_objects as go
import numpy as np

fig = go.Figure(data=[go.Scatter3d(x=x, y=y, z=z,
                                   mode='markers')])
fig.show()

'''
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
x = []
y = []
z = []
# Creating dataset
for i in range (0, len(center), 1):
    x.append(center[i][0])
    y.append(center[i][1])
    z.append(center[i][2])

# Creating figure
fig = plt.figure(figsize=(10, 7))
ax = plt.axes(projection="3d")

# Creating plot
ax.scatter3D(x, y, z, color="green")
plt.title("simple 3D scatter plot")

# show plot
plt.show()



scale = 1
vertices = vertices*scale

faces = faces
fig = go.Figure()

fig.add_trace(mesh_plot(vertices, faces))
fig.add_trace(wireframe_plot(vertices, faces));

fig.update_layout(title_text='Icosphere', height=600, width=600)
fig.show()

'''