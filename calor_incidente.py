# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 10:55:14 2022

@author: T2F-7
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import euler_angle
import terra
import orientacao_quat
import periodo_orbital
import propagador_orbital_mk3
import sup_terra

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

Posi_XYZ = propagador_orbital_mk3.propagador_orbital(rp, ecc, Raan, true_anomaly, inc, arg_per, num_orbita, 0)

xyz = orientacao_quat.orientacao_quat(Ia, Ib, Ic, PSIP, TETAP, PHIP, psi0, teta0, phi0, T_orbita)

ori_xyz = np.zeros((len(Posi_XYZ), 3))
K = len(xyz) / len(Posi_XYZ)
j = 0
for i in range(0, len(xyz), int(K)):
    ori_xyz[j][0] = xyz[i][0]
    ori_xyz[j][1] = xyz[i][1]
    ori_xyz[j][2] = xyz[i][2]
    j = j + 1
T = np.linspace(0, len(Posi_XYZ), len(Posi_XYZ))
R_terra = terra.terra(Raio_terra, 10)

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

Qs1 = []
Qs2 = []
Qs3 = []
Qs3 = []
Qs4 = []
Qs5 = []
Qs6 = []

for i in range(0, len(Posi_XYZ), 1):

    PSI = np.arccos(np.dot(Posi_XYZ[i], Vs) / (np.linalg.norm(Posi_XYZ[i]) * np.linalg.norm(Vs)))
    QSI = np.arcsin(Raio_terra / float(np.sqrt(Posi_XYZ[i][0] ** 2 + Posi_XYZ[i][1] ** 2 + Posi_XYZ[i][2] ** 2)))

    if PSI + QSI < np.pi:

        A1 = rotacao_Euler_2(Ni[0], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
        A2 = rotacao_Euler_2(Ni[1], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
        A3 = rotacao_Euler_2(Ni[2], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
        A4 = rotacao_Euler_2(Ni[3], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
        A5 = rotacao_Euler_2(Ni[4], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
        A6 = rotacao_Euler_2(Ni[5], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])

        k1 = np.dot(A1, Vs)
        k2 = np.dot(A2, Vs)
        k3 = np.dot(A3, Vs)
        k4 = np.dot(A4, Vs)
        k5 = np.dot(A5, Vs)
        k6 = np.dot(A6, Vs)

        if k1 > 0:
            qs1 = ai * Is * k1
            Qs1.append(qs1)
        else:
            Qs1.append(0)
        if k2 > 0:
            qs2 = ai * Is * k2
            Qs2.append(qs2)
        else:
            Qs2.append(0)

        if k3 > 0:
            qs3 = ai * Is * k3
            Qs3.append(qs3)
        else:
            Qs3.append(0)
        if k4 > 0:
            qs4 = ai * Is * k4
            Qs4.append(qs4)
        else:
            Qs4.append(0)
        if k5 > 0:
            qs5 = ai * Is * k5
            Qs5.append(qs5)
        else:
            Qs5.append(0)
        if k6 > 0:
            qs6 = ai * Is * k6
            Qs6.append(qs6)
        else:
            Qs6.append(0)

    else:
        Qs1.append(0)
        Qs2.append(0)
        Qs3.append(0)
        Qs4.append(0)
        Qs5.append(0)
        Qs6.append(0)

Qtotal = [[Qs1],
          [Qs2],
          [Qs3],
          [Qs4],
          [Qs5],
          [Qs6]]
'''    
Q_total = []    
for i in range (0, 200, 1):
    Q_total.append(Qs1[i] + Qs2[i] + Qs3[i] + Qs4[i] + Qs5[i] + Qs6[i])
'''

'''
plt.xlabel("Ponto da orbita")
plt.ylabel("Calor incidente em cada face [W/m^2]")
plt.plot(T, Qs1, color ='green', label='N1')
plt.plot(T, Qs2, color = 'blue', label='N2')
plt.plot(T, Qs3, color = 'cyan', label='N3')
plt.plot(T, Qs4, color = 'yellow', label='N4')
plt.plot(T, Qs5, color = 'red', label='N5')
plt.plot(T, Qs6, color = 'magenta', label='N6')   
'''

divisao = int(500)
Terra = terra.terra(Raio_terra, divisao)
As = sup_terra.sup_terra(Raio_terra, divisao)
x = []
y = []
z = []
for i in range(0, len(Terra), 1):
    x.append(Terra[i][0])
    y.append(Terra[i][1])
    z.append(Terra[i][2])
'''
fig = plt.figure(figsize=(10,7))
ax = plt.axes(projection = "3d")
ax.scatter3D(x,y,z)
plt.show()   
'''

''' Albedo '''

Halb1 = 0
Halb2 = 0
Halb3 = 0
Halb4 = 0
Halb5 = 0
Halb6 = 0

Qalb1 = []
Qalb2 = []
Qalb3 = []
Qalb4 = []
Qalb5 = []
Qalb6 = []

H1 = []
H2 = []
H3 = []
H4 = []
H5 = []
H6 = []

Rhok = []
prod_vet = []

for i in range(0, len(Posi_XYZ), 1):

    PSI = np.arccos(np.dot(Posi_XYZ[i], Vs) / (np.linalg.norm(Posi_XYZ[i]) * np.linalg.norm(Vs)))
    QSI = np.arcsin(Raio_terra / float(np.sqrt(Posi_XYZ[i][0] ** 2 + Posi_XYZ[i][1] ** 2 + Posi_XYZ[i][2] ** 2)))

    A1 = rotacao_Euler_2(Ni[0], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
    A2 = rotacao_Euler_2(Ni[1], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
    A3 = rotacao_Euler_2(Ni[2], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
    A4 = rotacao_Euler_2(Ni[3], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
    A5 = rotacao_Euler_2(Ni[4], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
    A6 = rotacao_Euler_2(Ni[5], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])

    if PSI + QSI < np.pi:

        for k in range(0, len(Terra), 1):

            rhok1 = np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A1)
            rhok2 = np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A2)
            rhok3 = np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A3)
            rhok4 = np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A4)
            rhok5 = np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A5)
            rhok6 = np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A6)

            rhok = np.array(Posi_XYZ[i]) - np.array(Terra[k])

            Rhok.append(rhok)
            prod_vet.append(np.dot(rhok, Terra[k]))

            C_bek = np.dot(Vs, Terra[k]) / (np.linalg.norm(Vs) * np.linalg.norm(Terra[k]))

            if np.dot(rhok1, Terra[k]) > 0:

                C_aek1 = np.dot(rhok1, Terra[k]) / (np.linalg.norm(rhok1) * np.linalg.norm(Terra[k]))
                C_aik1 = (np.dot(-rhok1, A1)) / (np.linalg.norm(rhok1) * np.linalg.norm(A1))

                if C_aik1 > 0 and C_aek1 > 0 and C_bek > 0:

                    Balb = float(1.0)

                    Halb1 = Halb1 + (As[k] * ((C_aek1 * C_bek * C_aik1) / (np.pi * np.linalg.norm(rhok1) ** 2)) * Balb)
                else:
                    Balb = float(0.0)

                    Halb1 = Halb1 + (As[k] * ((C_aek1 * C_bek * C_aik1) / (np.pi * np.linalg.norm(rhok1) ** 2)) * Balb)

            if np.dot(rhok2, Terra[k]) > 0:

                C_aek2 = np.dot(rhok2, Terra[k]) / (np.linalg.norm(rhok2) * np.linalg.norm(Terra[k]))
                C_aik2 = (np.dot(-rhok2, A2)) / (np.linalg.norm(rhok2) * np.linalg.norm(A2))

                if C_aik2 > 0 and C_aek2 > 0 and C_bek > 0:

                    Balb = float(1.0)

                    Halb2 = Halb2 + (As[k] * ((C_aek2 * C_bek * C_aik2) / (np.pi * np.linalg.norm(rhok2) ** 2)) * Balb)
                else:
                    Balb = float(0.0)

                    Halb2 = Halb2 + (As[k] * ((C_aek2 * C_bek * C_aik2) / (np.pi * np.linalg.norm(rhok2) ** 2)) * Balb)

            if np.dot(rhok3, Terra[k]) > 0:

                C_aek3 = np.dot(rhok3, Terra[k]) / (np.linalg.norm(rhok3) * np.linalg.norm(Terra[k]))
                C_aik3 = (np.dot(-rhok3, A3)) / (np.linalg.norm(rhok3) * np.linalg.norm(A3))

                if C_aik3 > 0 and C_aek3 > 0 and C_bek > 0:

                    Balb = float(1.0)

                    Halb3 = Halb3 + (As[k] * ((C_aek3 * C_bek * C_aik3) / (np.pi * np.linalg.norm(rhok3) ** 2)) * Balb)

                else:
                    Balb = float(0.0)

                    Halb3 = Halb3 + (As[k] * ((C_aek3 * C_bek * C_aik3) / (np.pi * np.linalg.norm(rhok3) ** 2)) * Balb)

            if np.dot(rhok4, Terra[k]) > 0:

                C_aek4 = np.dot(rhok4, Terra[k]) / (np.linalg.norm(rhok4) * np.linalg.norm(Terra[k]))
                C_aik4 = (np.dot(-rhok4, A4)) / (np.linalg.norm(rhok4) * np.linalg.norm(A4))

                if C_aik4 > 0 and C_aek4 > 0 and C_bek > 0:

                    Balb = float(1.0)

                    Halb4 = Halb4 + (As[k] * ((C_aek4 * C_bek * C_aik4) / (np.pi * np.linalg.norm(rhok4) ** 2)) * Balb)
                else:
                    Balb = float(0.0)

                    Halb4 = Halb4 + (As[k] * ((C_aek4 * C_bek * C_aik4) / (np.pi * np.linalg.norm(rhok4) ** 2)) * Balb)

            if np.dot(rhok5, Terra[k]) > 0:

                C_aek5 = np.dot(rhok5, Terra[k]) / (np.linalg.norm(rhok5) * np.linalg.norm(Terra[k]))
                C_aik5 = (np.dot(-rhok5, A5)) / (np.linalg.norm(rhok5) * np.linalg.norm(A5))

                if C_aik5 > 0 and C_aek5 > 0 and C_bek > 0:

                    Balb = float(1.0)

                    Halb5 = Halb5 + (As[k] * ((C_aek5 * C_bek * C_aik5) / (np.pi * np.linalg.norm(rhok5) ** 2)) * Balb)
                else:
                    Balb = float(0.0)

                    Halb5 = Halb5 + (As[k] * ((C_aek5 * C_bek * C_aik5) / (np.pi * np.linalg.norm(rhok5) ** 2)) * Balb)

            if np.dot(rhok6, Terra[k]) > 0:

                C_aek6 = np.dot(rhok6, Terra[k]) / (np.linalg.norm(rhok6) * np.linalg.norm(Terra[k]))
                C_aik6 = (np.dot(-rhok6, A6)) / (np.linalg.norm(rhok6) * np.linalg.norm(A6))

                if C_aik6 > 0 and C_aek6 > 0 and C_bek > 0:

                    Balb = float(1.0)

                    Halb6 = Halb6 + (As[k] * ((C_aek6 * C_bek * C_aik6) / (np.pi * np.linalg.norm(rhok6) ** 2)) * Balb)
                else:
                    Balb = float(0.0)

                    Halb6 = Halb6 + (As[k] * ((C_aek6 * C_bek * C_aik6) / (np.pi * np.linalg.norm(rhok6) ** 2)) * Balb)

        H1.append(Halb1)
        Qalb1.append(ai * gama * Is * Halb1)
        Halb1 = 0

        H2.append(Halb2)
        Qalb2.append(ai * gama * Is * Halb2)
        Halb2 = 0

        H3.append(Halb3)
        Qalb3.append(ai * gama * Is * Halb3)
        Halb3 = 0

        H4.append(Halb4)
        Qalb4.append(ai * gama * Is * Halb4)
        Halb4 = 0

        H5.append(Halb5)
        Qalb5.append(ai * gama * Is * Halb5)
        Halb5 = 0

        H6.append(Halb6)
        Qalb6.append(ai * gama * Is * Halb6)
        Halb6 = 0
    else:
        Qalb1.append(0)
        Qalb2.append(0)
        Qalb3.append(0)
        Qalb4.append(0)
        Qalb5.append(0)
        Qalb6.append(0)

''' Radiacao da terra '''

Hrad1 = 0
Hrad2 = 0
Hrad3 = 0
Hrad4 = 0
Hrad5 = 0
Hrad6 = 0

Qrad1 = []
Qrad2 = []
Qrad3 = []
Qrad4 = []
Qrad5 = []
Qrad6 = []

R1 = []
R2 = []
R3 = []
R4 = []
R5 = []
R6 = []

Rrhok = []
prod_vet2 = []

for i in range(0, len(Posi_XYZ), 1):

    A1 = rotacao_Euler_2(Ni[0], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
    A2 = rotacao_Euler_2(Ni[1], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
    A3 = rotacao_Euler_2(Ni[2], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
    A4 = rotacao_Euler_2(Ni[3], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
    A5 = rotacao_Euler_2(Ni[4], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])
    A6 = rotacao_Euler_2(Ni[5], ori_xyz[i][0], ori_xyz[i][1], ori_xyz[i][2])

    for k in range(0, len(Terra), 1):

        Rhok1 = np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A1)
        Rhok2 = np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A2)
        Rhok3 = np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A3)
        Rhok4 = np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A4)
        Rhok5 = np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A5)
        Rhok6 = np.array(Posi_XYZ[i]) - np.array(Terra[k]) + np.array(A6)

        rhok = np.array(Posi_XYZ[i]) - np.array(Terra[k])

        Rrhok.append(rhok)
        prod_vet.append(np.dot(rhok, Terra[k]))

        if np.dot(Rhok1, Terra[k]) > 0:

            C_aek1 = np.dot(Rhok1, Terra[k]) / (np.linalg.norm(Rhok1) * np.linalg.norm(Terra[k]))
            C_aik1 = (np.dot(-Rhok1, A1)) / (np.linalg.norm(Rhok1) * np.linalg.norm(A1))

            if C_aik1 > 0 and C_aek1 > 0:

                Balb = float(1.0)

                Hrad1 = Hrad1 + (As[k] * ((C_aek1 * C_aik1) / (np.pi * np.linalg.norm(Rhok1) ** 2)) * Balb)
            else:
                Balb = float(0.0)

                Hrad1 = Hrad1 + (As[k] * ((C_aek1 * C_aik1) / (np.pi * np.linalg.norm(Rhok1) ** 2)) * Balb)

        if np.dot(Rhok2, Terra[k]) > 0:

            C_aek2 = np.dot(Rhok2, Terra[k]) / (np.linalg.norm(Rhok2) * np.linalg.norm(Terra[k]))
            C_aik2 = (np.dot(-Rhok2, A2)) / (np.linalg.norm(Rhok1) * np.linalg.norm(A2))

            if C_aik2 > 0 and C_aek2 > 0:

                Balb = float(1.0)

                Hrad2 = Hrad2 + (As[k] * ((C_aek2 * C_aik2) / (np.pi * np.linalg.norm(Rhok2) ** 2)) * Balb)
            else:
                Balb = float(0.0)

                Hrad2 = Hrad2 + (As[k] * ((C_aek2 * C_aik2) / (np.pi * np.linalg.norm(Rhok2) ** 2)) * Balb)

        if np.dot(Rhok3, Terra[k]) > 0:

            C_aek3 = np.dot(Rhok3, Terra[k]) / (np.linalg.norm(Rhok3) * np.linalg.norm(Terra[k]))
            C_aik3 = (np.dot(-Rhok3, A3)) / (np.linalg.norm(Rhok3) * np.linalg.norm(A1))

            if C_aik3 > 0 and C_aek3 > 0:

                Balb = float(1.0)

                Hrad3 = Hrad3 + (As[k] * ((C_aek3 * C_aik3) / (np.pi * np.linalg.norm(Rhok3) ** 2)) * Balb)
            else:
                Balb = float(0.0)

                Hrad3 = Hrad3 + (As[k] * ((C_aek3 * C_aik3) / (np.pi * np.linalg.norm(Rhok3) ** 2)) * Balb)

        if np.dot(Rhok4, Terra[k]) > 0:

            C_aek4 = np.dot(Rhok4, Terra[k]) / (np.linalg.norm(Rhok4) * np.linalg.norm(Terra[k]))
            C_aik4 = (np.dot(-Rhok4, A4)) / (np.linalg.norm(Rhok4) * np.linalg.norm(A4))

            if C_aik4 > 0 and C_aek4 > 0:

                Balb = float(1.0)

                Hrad4 = Hrad4 + (As[k] * ((C_aek4 * C_aik4) / (np.pi * np.linalg.norm(Rhok4) ** 2)) * Balb)
            else:
                Balb = float(0.0)

                Hrad4 = Hrad4 + (As[k] * ((C_aek4 * C_aik4) / (np.pi * np.linalg.norm(Rhok4) ** 2)) * Balb)

        if np.dot(Rhok5, Terra[k]) > 0:

            C_aek5 = np.dot(Rhok5, Terra[k]) / (np.linalg.norm(Rhok5) * np.linalg.norm(Terra[k]))
            C_aik5 = (np.dot(-Rhok5, A5)) / (np.linalg.norm(Rhok5) * np.linalg.norm(A5))

            if C_aik5 > 0 and C_aek5 > 0:

                Balb = float(1.0)

                Hrad5 = Hrad5 + (As[k] * ((C_aek5 * C_aik5) / (np.pi * np.linalg.norm(Rhok5) ** 2)) * Balb)
            else:
                Balb = float(0.0)

                Hrad5 = Hrad5 + (As[k] * ((C_aek5 * C_aik5) / (np.pi * np.linalg.norm(Rhok5) ** 2)) * Balb)

        if np.dot(Rhok6, Terra[k]) > 0:

            C_aek6 = np.dot(Rhok6, Terra[k]) / (np.linalg.norm(Rhok6) * np.linalg.norm(Terra[k]))
            C_aik6 = (np.dot(-Rhok6, A6)) / (np.linalg.norm(Rhok6) * np.linalg.norm(A6))

            if C_aik6 > 0 and C_aek6 > 0:

                Balb = float(1.0)

                Hrad6 = Hrad6 + (As[k] * ((C_aek6 * C_aik6) / (np.pi * np.linalg.norm(Rhok6) ** 2)) * Balb)
            else:
                Balb = float(0.0)

                Hrad6 = Hrad6 + (As[k] * ((C_aek6 * C_aik6) / (np.pi * np.linalg.norm(Rhok6) ** 2)) * Balb)

    R1.append(Hrad1)
    Qrad1.append(e * Ir * (Hrad1))
    Hrad1 = 0

    R2.append(Hrad2)
    Qrad2.append(e * Ir * (Hrad2))
    Hrad2 = 0

    R3.append(Hrad3)
    Qrad3.append(e * Ir * (Hrad3))
    Hrad3 = 0

    R4.append(Hrad4)
    Qrad4.append(e * Ir * (Hrad4))
    Hrad4 = 0

    R5.append(Hrad5)
    Qrad5.append(e * Ir * (Hrad5))
    Hrad5 = 0

    R5.append(Hrad5)
    Qrad5.append(e * Ir * (Hrad5))
    Hrad5 = 0

    R6.append(Hrad6)
    Qrad6.append(e * Ir * (Hrad6))
    Hrad6 = 0

Qt1 = []
Qt2 = []
Qt3 = []
Qt4 = []
Qt5 = []
Qt6 = []

for i in range(0, len(Qalb1), 1):
    Qt1.append(Qalb1[i] + Qs1[i] + Qrad1[i])
    Qt2.append(Qalb2[i] + Qs2[i] + Qrad2[i])
    Qt3.append(Qalb3[i] + Qs3[i] + Qrad3[i])
    Qt4.append(Qalb4[i] + Qs4[i] + Qrad4[i])
    Qt5.append(Qalb5[i] + Qs5[i] + Qrad5[i])
    Qt6.append(Qalb6[i] + Qs6[i] + Qrad6[i])

fig = plt.figure()
plt.xlabel("Ponto da orbita")
plt.ylabel("Calor incidente em cada face [W/m^2]")
plt.plot(T, Qt1, color='green', label='N1')
plt.plot(T, Qt2, color='blue', label='N2')
plt.plot(T, Qt3, color='cyan', label='N3')
plt.plot(T, Qt4, color='yellow', label='N4')
plt.plot(T, Qt5, color='red', label='N5')
plt.plot(T, Qt6, color='magenta', label='N6')
plt.legend()
plt.show()
