def posi_ori(posi, ori):
      orb_sat = posi
      xyz = ori
      K = int(len(xyz) / len(orb_sat))
      print(K)
      A = []
      for i in range(0, len(orb_sat), 1):
            x = orb_sat.iloc[i, 0]
            y = orb_sat.iloc[i, 1]
            z = orb_sat.iloc[i, 2]
            j = i*K
            psi = xyz.iloc[j, 0]
            teta = xyz.iloc[j, 1]
            phi = xyz.iloc[j, 2]

            A.append([x, y, z, psi, teta, phi])
      posi_ori = pd.DataFrame(A, columns=['X', 'Y', 'Z', 'Psi', 'Teta', 'Phi'])
      return (posi_ori)

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

Posicao_orientacao = posi_ori(orb_sat, xyz)
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
df1 = pd.DataFrame(Ni, columns=['x', 'y', 'z'])
Posicao_orientacao = pd.concat([Posicao_orientacao, df1], axis=1)
names = [['N1_X', 'N1_Y', 'N1_Z'],
         ['N2_X', 'N2_Y', 'N2_Z'],
         ['N3_X', 'N3_Y', 'N3_Z'],
         ['N4_X', 'N4_Y', 'N4_Z'],
         ['N5_X', 'N5_Y', 'N5_Z'],
         ['N6_X', 'N6_Y', 'N6_Z']]
R = []
for j in range(0, len(Ni), 1):
      for i in range(0, len(Posicao_orientacao), 1):
            A = np.array([Posicao_orientacao.iloc[j][6],
                          Posicao_orientacao.iloc[j][7],
                          Posicao_orientacao.iloc[j][8]])

            ra = Posicao_orientacao.iloc[i][3]
            inc = Posicao_orientacao.iloc[i][4]
            ome = Posicao_orientacao.iloc[i][5]
            vetor = A
            R1 = np.array([[np.cos(ra), np.sin(ra), 0],
                           [-np.sin(ra), np.cos(ra), 0],
                           [0, 0, 1]])

            R2 = np.array([[1, 0, 0],
                           [0, np.cos(inc), np.sin(inc)],
                           [0, -np.sin(inc), np.cos(inc)]])

            R3 = np.array([[np.cos(ome), np.sin(ome), 0],
                           [-np.sin(ome), np.cos(ome), 0],
                           [0, 0, 1]])
            a = R2.dot(R1)
            Tci = R3.dot(a)
            A = (np.transpose(Tci).dot(vetor))
            R1 = A[0]
            R2 = A[1]
            R3 = A[2]
            R.append([R1, R2, R3])
      df2 = pd.DataFrame(R, columns=names[j])
      R = []
      Posicao_orientacao = pd.concat([Posicao_orientacao, df2], axis=1)

divisao = int(500)

vet_terra = terra.terra(Raio_terra, divisao)
As = sup_terra.sup_terra(Raio_terra, divisao)
concat_terra = pd.concat([vet_terra, As], axis=1)
Posicao_orientacao = pd.concat([Posicao_orientacao, concat_terra], axis=1)

Posicao_orientacao['rho1_X'] = Posicao_orientacao['X'] + Posicao_orientacao['N1_X'] - Posicao_orientacao['Terra_X']
Posicao_orientacao['rho1_Y'] = Posicao_orientacao['Y'] + Posicao_orientacao['N1_Y'] - Posicao_orientacao['Terra_Y']
Posicao_orientacao['rho1_Z'] = Posicao_orientacao['Z'] + Posicao_orientacao['N1_Z'] - Posicao_orientacao['Terra_Z']

Posicao_orientacao['rho2_X'] = Posicao_orientacao['X'] + Posicao_orientacao['N2_X'] - Posicao_orientacao['Terra_X']
Posicao_orientacao['rho2_Y'] = Posicao_orientacao['Y'] + Posicao_orientacao['N2_Y'] - Posicao_orientacao['Terra_Y']
Posicao_orientacao['rho2_Z'] = Posicao_orientacao['Z'] + Posicao_orientacao['N2_Z'] - Posicao_orientacao['Terra_Z']

Posicao_orientacao['rho3_X'] = Posicao_orientacao['X'] + Posicao_orientacao['N3_X'] - Posicao_orientacao['Terra_X']
Posicao_orientacao['rho3_Y'] = Posicao_orientacao['Y'] + Posicao_orientacao['N3_Y'] - Posicao_orientacao['Terra_Y']
Posicao_orientacao['rho3_Z'] = Posicao_orientacao['Z'] + Posicao_orientacao['N3_Z'] - Posicao_orientacao['Terra_Z']

Posicao_orientacao['rho4_X'] = Posicao_orientacao['X'] + Posicao_orientacao['N4_X'] - Posicao_orientacao['Terra_X']
Posicao_orientacao['rho4_Y'] = Posicao_orientacao['Y'] + Posicao_orientacao['N4_Y'] - Posicao_orientacao['Terra_Y']
Posicao_orientacao['rho4_Z'] = Posicao_orientacao['Z'] + Posicao_orientacao['N4_Z'] - Posicao_orientacao['Terra_Z']

Posicao_orientacao['rho5_X'] = Posicao_orientacao['X'] + Posicao_orientacao['N5_X'] - Posicao_orientacao['Terra_X']
Posicao_orientacao['rho5_Y'] = Posicao_orientacao['Y'] + Posicao_orientacao['N5_Y'] - Posicao_orientacao['Terra_Y']
Posicao_orientacao['rho5_Z'] = Posicao_orientacao['Z'] + Posicao_orientacao['N5_Z'] - Posicao_orientacao['Terra_Z']

Posicao_orientacao['rho6_X'] = Posicao_orientacao['X'] + Posicao_orientacao['N6_X'] - Posicao_orientacao['Terra_X']
Posicao_orientacao['rho6_Y'] = Posicao_orientacao['Y'] + Posicao_orientacao['N6_Y'] - Posicao_orientacao['Terra_Y']
Posicao_orientacao['rho6_Z'] = Posicao_orientacao['Z'] + Posicao_orientacao['N6_Z'] - Posicao_orientacao['Terra_Z']



























#Posicao_orientacao.to_csv('posicao.csv',sep=',')
