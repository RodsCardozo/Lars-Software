
import numpy as np
import matplotlib.pyplot as plt
import euler_angle
import terra
import orientacao_quat
import periodo_orbital
import propagador_orbital_mk3
import sup_terra
import pandas as pd
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
T_orbita = periodo_orbital.periodo_orbital(rp)
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

orb_sat = propagador_orbital_mk3.propagador_orbital(rp, ecc, Raan, true_anomaly, inc, arg_per, num_orbita, 0)

xyz = orientacao_quat.orientacao_quat(Ia, Ib, Ic, PSIP, TETAP, PHIP, psi0, teta0, phi0, T_orbita)


ori_xyz = np.zeros((len(orb_sat), 3))
K = len(xyz) / len(orb_sat)
print(K)

Vs = np.array([1, 0, 0])

Ai = [a * c,
      b * c,
      a * c,
      b * c,
      a * b,
      a * b]

Ni = [[1, 0, 0],
      [0, 1, 0],
      [-1, 0, 0],
      [0, -1, 0],
      [0, 0, -1],
      [0, 0, 1]]

divisao = int(500)

vet_terra = terra.terra(Raio_terra, divisao)
As = sup_terra.sup_terra(Raio_terra, divisao)

Terra = pd.concat([vet_terra,As], axis=1)
print(Terra)


























'''
rhok1 = []
rhok2 = []
rhok3 = []
rhok4 = []
rhok5 = []
rhok6 = []
for i in range(0, len(Posi_XYZ), 1):

    PSI = np.arccos(np.dot(Posi_XYZ[i], Vs) / (np.linalg.norm(Posi_XYZ[i]) * np.linalg.norm(Vs)))
    QSI = np.arcsin(Raio_terra / np.linalg.norm(Raio_terra))
    A1 = rotacao_Euler_2(Ni[0], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
    A2 = rotacao_Euler_2(Ni[1], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
    A3 = rotacao_Euler_2(Ni[2], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
    A4 = rotacao_Euler_2(Ni[3], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
    A5 = rotacao_Euler_2(Ni[4], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
    A6 = rotacao_Euler_2(Ni[5], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])

    if PSI + QSI < np.pi:

        for k in range(0, len(Terra), 1):

            rhok1.append(np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A1))
            rhok2.append(np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A2))
            rhok3.append(np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A3))
            rhok4.append(np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A4))
            rhok5.append(np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A5))
            rhok6.append(np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A6))
'''