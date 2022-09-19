import numpy as np
import matplotlib.pyplot as plt
from orientacao_quat import orientacao_quat
from periodo_orbital import periodo_orbital
from propagador_orbital import propagador_orbital
from sup_terra import sup_terra
from terra import terra
from euler_angle import rotacao_euler_2


m = float(3)  # massa do cubesat
a = float(0.1)  # comprimento do sat
b = float(0.1)  # largura do sat
c = float(0.2)  # altura do sat

Ia = (m / 12) * (b ** 2 + c ** 2)  # momento de inercia na direcao x
Ib = (m / 12) * (a ** 2 + c ** 2)  # momento de inercia na direcao y
Ic = (m / 12) * (a ** 2 + b ** 2)  # momento de inercia na direcao z

rp = 7000  # semi eixo maior
ecc = float(0.0)  # ecentricidade da orbita
Raan = float(0.0)  # ascencao direita do nodo ascendente
arg_per = (float(0.0))  # argumento do perigeu
true_anomaly = (float(0.0))  # anomalia verdadeira
inc = (float(51.6))  # inclinacao
mu = 398600  # constante gravitacional da terra
J2 = 1.08263e-3  # zona harmonica 2
Raio_terra = float(6371)  # raio da terra
num_orbita = 1  # numero de obitas
Is = 1367.0
Ir = 267.0
e = 1.0
ai = 1.0

T_orbita = periodo_orbital(rp)
passo = 10000
gama = 0.3

PSIP = 0.0
TETAP = 0.0
PHIP = (2 * np.pi) / T_orbita
psi0 = Raan
teta0 = inc
phi0 = 0.0

psi = []
teta = []
phi = []
xyz = []

Posi_XYZ = propagador_orbital(rp, ecc, Raan, true_anomaly, inc, arg_per, num_orbita, 0)

orient_xyz = orientacao_quat(Ia, Ib, Ic, PSIP, TETAP, PHIP, psi0, teta0, phi0, T_orbita)

ori = np.zeros((len(Posi_XYZ), 3))
K = len(orient_xyz) / len(Posi_XYZ)
j = 0
for i in range(0, len(orient_xyz), int(K)):
    ori[j][0] = orient_xyz[i][0]
    ori[j][1] = orient_xyz[i][1]
    ori[j][2] = orient_xyz[i][2]
    j = j + 1

T = np.linspace(0, len(Posi_XYZ), len(Posi_XYZ))
R_terra = terra(Raio_terra, 10)
divisao = 10
Terra = np.array(terra(Raio_terra, divisao))
As = (sup_terra(Raio_terra, divisao))
x = []
y = []
z = []
for i in range(0, len(Terra), 1):
    x.append(Terra[i][0])
    y.append(Terra[i][1])
    z.append(Terra[i][2])
VETOR_TERRA = [Terra, Posi_XYZ]


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
m = ['green', 'red']
for i in range(0, len(VETOR_TERRA), 1):
    for j in range(0, len(VETOR_TERRA[i]), 1):
        xs = VETOR_TERRA[i][j][0]
        ys = VETOR_TERRA[i][j][1]
        zs = VETOR_TERRA[i][j][2]
        ax.scatter(xs, ys, zs, color=m[i])

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()

import plotly.io as pio

pio.renderers.default = 'browser'
import plotly.graph_objects as go

m = ['green', 'red']
xs = []
ys = []
zs = []
corzinha = ['green', 'red']
for i in range(0, len(VETOR_TERRA), 1):
    for j in range(0, len(VETOR_TERRA[i]), 1):
        xs.append(np.array(VETOR_TERRA[i][j][0]))
        ys.append(np.array(VETOR_TERRA[i][j][1]))
        zs.append(np.array(VETOR_TERRA[i][j][2]))
        # ax.scatter(xs, ys, zs, color=m[i])
    m = corzinha[i]
fig = go.Figure(data=[go.Scatter3d(x=xs, y=ys, z=zs,
                                       mode='markers',
                                       marker=dict(size=5, color=corzinha, colorscale='Viridis', opacity=0.8))])
fig.show()


Vs = np.array([1,0,0])

Ai = [a*c,
      b*c,
      a*c,
      b*c,
      a*b,
      a*b]

Ni = [[1,0,0],
      [0,1,0],
      [-1,0,0],
      [0,-1,0],
      [0,0,-1],
      [0,0, 1]]

rhok1 = []
rhok2 = []
rhok3 = []
rhok4 = []
rhok5 = []
rhok6 = []
for i in range(0, len(Posi_XYZ), 1):

    A1 = rotacao_euler_2(Ni[0], ori[i][0], ori[i][1], ori[i][2])
    A2 = rotacao_euler_2(Ni[1], ori[i][0], ori[i][1], ori[i][2])
    A3 = rotacao_euler_2(Ni[2], ori[i][0], ori[i][1], ori[i][2])
    A4 = rotacao_euler_2(Ni[3], ori[i][0], ori[i][1], ori[i][2])
    A5 = rotacao_euler_2(Ni[4], ori[i][0], ori[i][1], ori[i][2])
    A6 = rotacao_euler_2(Ni[5], ori[i][0], ori[i][1], ori[i][2])

    for k in range(0, len(Terra), 1):
        rhok1.append(np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A1))
        rhok2.append(np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A2))
        rhok3.append(np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A3))
        rhok4.append(np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A4))
        rhok5.append(np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A5))
        rhok6.append(np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A6))

