
import numpy as np
import matplotlib.pyplot as plt
import periodo_orbital
from propagador_orbital_mk4 import propagador_orbital
import pandas as pd
import icosphere
import numpy as np
import pandas as pd
from scipy.integrate import odeint

m = float(3)  # massa do cubesat
a = float(0.1)  # comprimento do sat
b = float(0.1)  # largura do sat
c = float(0.2)  # altura do sat

Ia = (m / 12) * (b ** 2 + c ** 2)  # momento de inercia na direcao x
Ib = (m / 12) * (a ** 2 + c ** 2)  # momento de inercia na direcao y
Ic = (m / 12) * (a ** 2 + b ** 2)  # momento de inercia na direcao z

rp = 7000  # semi eixo maior
ecc = float(0.0001)  # ecentricidade da orbita
Raan = float(0.0)  # ascencao direita do nodo ascendente
arg_per = (float(0.0))  # argumento do perigeu
true_anomaly = (float(0.0))  # anomalia verdadeira
inc = (float(15.0))  # inclinacao
mu = 398600  # constante gravitacional da terra
J2 = 1.08263e-3  # zona harmonica 2
Raio_terra = float(6371)  # raio da terra
num_orbita = 1  # numero de obitas
Is = 1367.0
Ir = 267.0
e = 1.0
ai = 1.0
T_orbita = periodo_orbital.periodo_orbital(rp)
passo = 10000
gama = 0.3
PSIP = 0.0
TETAP = 0.0
PHIP = (2 * np.pi) / T_orbita
psi0 = 0.0
teta0 = inc
phi0 = true_anomaly
prop_orb = propagador_orbital(rp, ecc, Raan, arg_per, true_anomaly, inc, 10, psi0, teta0, phi0, PSIP, TETAP, PHIP)



import plotly.express as px
fig = px.scatter_3d(prop_orb, x='X_ECI', y='Y_ECI', z='Z_ECI')
fig.show()
'''import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = plt.axes(projection="3d")
a1 = ax.scatter3D(prop_orb['X_ECI'], prop_orb['Y_ECI'], prop_orb['Z_ECI'])
ax.set_xlabel('x [km]')
ax.set_ylabel('y [km]')
ax.set_zlabel('z [km]')
plt.show()
'''