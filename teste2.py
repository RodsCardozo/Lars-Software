import numpy as np
import pandas as pd
from scipy.integrate import odeint
import periodo_orbital
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
inc = (float(51.0))  # inclinacao
mu = 398600  # constante gravitacional da terra
J2 = 1.08263e-3  # zona harmonica 2
Raio_terra = float(6371)  # raio da terra
num_orbita = 1  # numero de obitas
Is = 1367.0
Ir = 267.0
e = 1.0
ai = 1.0
T_orbita = periodo_orbital.periodo_orbital(rp)
print(T_orbita)
passo = 10000
gama = 0.3
PSIP = 0.0
TETAP = 0.0
PHIP = (2 * np.pi) / T_orbita
psi0 = 0.0
teta0 = inc
phi0 = true_anomaly
PTP = []

psi = float(np.radians(psi0))  # angulo inicial de PSI

teta = float(np.radians(teta0))  # angulo inicial de TETA

phi = float(np.radians(phi0))  # angulo inicial de PHI

Ix3 = float(Ia)  # momento de incercia no eixo x

Iy3 = float(Ib)  # momento de incercia no eixo y

Iz3 = float(Ic)  # momento de incercia no eixo z

psip = float(PSIP)  # velocidade angular do angulo PSI

tetap = float(TETAP)  # velocidade angular do angulo TETA

phip = float(PHIP)  # velocidade angular do angulo PHI

wx3_i = float(-psip * np.sin(teta) * np.cos(phi) + tetap * np.sin(phi))  # velocidade angular do corpo em x

wy3_i = float(psip * np.sin(teta) * np.sin(phi) + tetap * np.cos(phi))  # velocidade angular do corpo em y

wz3_i = float(psip * np.cos(teta) + phip)  # velocidade angular do corpo em z

q0 = float((np.cos(psi / 2) * np.cos(teta / 2) * np.cos(phi / 2) - np.sin(psi / 2) * np.cos(teta / 2) * np.sin(
    phi / 2)))  # quaternion q0

q1 = float((np.cos(psi / 2) * np.sin(teta / 2) * np.sin(phi / 2) - np.sin(psi / 2) * np.sin(teta / 2) * np.cos(
    phi / 2)))  # quaternion q1

q2 = float((np.cos(psi / 2) * np.sin(teta / 2) * np.cos(phi / 2) + np.sin(psi / 2) * np.sin(teta / 2) * np.sin(
    phi / 2)))  # quaternion q2

q3 = float((np.sin(psi / 2) * np.cos(teta / 2) * np.cos(phi / 2) + np.cos(psi / 2) * np.cos(teta / 2) * np.sin(
    phi / 2)))  # quaternion q3

qi = [q0, q1, q2, q3, wx3_i, wy3_i, wz3_i]  # condicoes iniciais da integracao


def quat(q, t, Ix3, Iy3, Iz3):  # funcao para integrar

    Ix3 = float(Ix3)

    Iy3 = float(Iy3)

    Iz3 = float(Iz3)

    q0, q1, q2, q3, wx3, wy3, wz3 = q

    dqdt = [0.5 * (q0 * 0 - q1 * wx3 - q2 * wy3 - q3 * wz3),
            0.5 * (q1 * 0 + q0 * wx3 - q3 * wy3 + q2 * wz3),
            0.5 * (q2 * 0 + q3 * wx3 + q0 * wy3 - q1 * wz3),
            0.5 * (q3 * 0 - q2 * wx3 + q1 * wy3 + q0 * wz3),
            ((Iy3 - Iz3) / Ix3) * wy3 * wz3,
            ((Iz3 - Ix3) / Iy3) * wz3 * wx3,
            ((Ix3 - Iy3) / Iz3) * wx3 * wy3]
    return dqdt

Time_step = T_orbita
passo = 50000
t = np.linspace(0, Time_step, passo)

sol2 = odeint(quat, qi, t, args=(Ix3, Iy3, Iz3))  # integracao do conjunto de EDO
print((sol2))
Psi1 = []
Teta1 = []
Phi1 = []
for i in range(0, len(t), 1):
    x = float(2 * (sol2[i][2] * sol2[i][3] - sol2[i][0] * sol2[i][1]))
    y = float(2 * (sol2[i][1] * sol2[i][3] + sol2[i][0] * sol2[i][2]))
    Psi1.append((np.arctan2(x, y)))
for i in range(0, len(t), 1):
    B = float(2 * (sol2[i][0] ** 2 + sol2[i][3] ** 2) - 1)
    Teta1.append((np.arcsin(B)))

for i in range(0, len(t), 1):
    Phix = float(2 * (sol2[i][1] * sol2[i][3] - sol2[i][0] * sol2[i][2]))
    Phiy = float(2 * (sol2[i][2] * sol2[i][3] + sol2[i][0] * sol2[i][1]))
    Phi1.append((np.arctan2(-Phiy, Phix)))

for i in range(0, len(Phi1), 1):
    PTP.append([Psi1[i], Teta1[i], Phi1[i]])

Psi_teta_phi = pd.DataFrame(PTP, columns=['Psi', 'Teta', 'Phi'])

print(Psi_teta_phi)
Psi_teta_phi.to_csv('Psi_teta_phi.csv',sep=',')