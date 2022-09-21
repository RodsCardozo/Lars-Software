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

Posi_XYZ = propagador_orbital(rp, ecc, Raan, true_anomaly, inc, arg_per, num_orbita, 1)


