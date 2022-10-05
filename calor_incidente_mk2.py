from mpl_toolkits.mplot3d import Axes3D


def posi_ori(posi, ori):
      import numpy as np
      import pandas as pd
      orb_sat = posi
      xyz = ori
      K = int(len(xyz) / len(orb_sat))
      print(K)
      A = []
      ''' for i in range(0, len(orb_sat), 1):
            x = orb_sat.iloc[i][0]
            y = orb_sat.iloc[i][1]
            z = orb_sat.iloc[i][2]
            j = i*K
            psi = xyz.iloc[j][0]
            teta = xyz.iloc[j][1]
            phi = xyz.iloc[j][2]'''
      ori_xyz = np.zeros((len(orb_sat), 3))
      j = 0
      i = 0
      '''for i in range(0, len(xyz), int(K)):'''
      while j <= len(orb_sat):
            x = orb_sat.iloc[j,0]
            y = orb_sat.iloc[j,1]
            z = orb_sat.iloc[j, 2]
            psi = xyz.iloc[i, 0]
            teta = xyz.iloc[i, 1]
            phi = xyz.iloc[i, 2]
            j = j + 1
            i = i + K
            A.append([x, y, z, psi, teta, phi])
      posi_ori = pd.DataFrame(A, columns=['X', 'Y', 'Z', 'Psi', 'Teta', 'Phi'])
      return (posi_ori)
def area(vertices, faces, Raio):
    import numpy as np
    a1 = vertices[faces[0][0]]*Raio
    a2 = vertices[faces[0][1]]*Raio
    a = np.array(a1) - np.array(a2)
    e = np.linalg.norm(a)
    A = ((3)**(1/2)/2)*e**2
    return A
def orientacao_quat(Ia, Ib, Ic, PSIP, TETAP, PHIP, Psi, Teta, Phi, Time_step):

    """
    & Ia = Momento de inercia na direcao x
    & Ib = Momento de inercia na direcao y
    & Ic = Momento de inercia na direcao z
    & PSIP = Velocidade Angular do angulo PSI (3)
    & TETAP = Velocidade Angular do angulo TETA (1)
    & PHIP = Velocidade Angular do angulo PHI (3)
    & PSI = Angulo PSI (3)
    & TETA = Angulo TETA (1)
    & PHI = Angulo PSI (3)
    & Time_step = Tempo total de simulacao

    """

    PTP = []

    psi = float(np.radians(Psi))  # angulo inicial de PSI

    teta = float(np.radians(Teta))  # angulo inicial de TETA

    phi = float(np.radians(Phi))  # angulo inicial de PHI

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

    passo = 200000
    t = np.linspace(0, Time_step, passo)

    sol2 = odeint(quat, qi, t, args=(Ix3, Iy3, Iz3))  # integracao do conjunto de EDO
    Psi1 = []
    Teta1 = []
    Phi1 = []
    for i in range(0, len(t), 1):
        x = float(2 * (sol2[i][2] * sol2[i][3] - sol2[i][0] * sol2[i][1]))
        y = float(2 * (sol2[i][1] * sol2[i][3] + sol2[i][0] * sol2[i][2]))
        if y == 0:
            Psi1.append((0))
        else:
            Psi1.append((np.arctan2(x,y)))
    for i in range(0, len(t), 1):
        B = float(2 * (sol2[i][0] ** 2 + sol2[i][3] ** 2) - 1)
        if B >= 1:
            Teta1.append((np.arcsin(1)))
        else:
            Teta1.append((np.arcsin(B)))

    for i in range(0, len(t), 1):
        Phix = float(2 * (sol2[i][1] * sol2[i][3] - sol2[i][0] * sol2[i][2]))
        Phiy = float(2 * (sol2[i][2] * sol2[i][3] + sol2[i][0] * sol2[i][1]))
        if Phix == 0:
            Phi1.append((0))
        else:
            Phi1.append((np.arctan2(Phiy,-Phix)))

    for i in range(0, len(Phi1), 1):
        PTP.append([Psi1[i], Teta1[i], Phi1[i]])

    Psi_teta_phi = pd.DataFrame(PTP, columns=['Psi', 'Teta', 'Phi'])

    return (Psi_teta_phi)

import numpy as np
import matplotlib.pyplot as plt
#import orientacao_quat
import periodo_orbital
import propagador_orbital_mk3
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
passo = 10000
gama = 0.3
PSIP = 0.0
TETAP = 0.0
PHIP = (2 * np.pi) / T_orbita
psi0 = 0.0
teta0 = inc
phi0 = true_anomaly
psi = []
teta = []
phi = []
print('Calculando posicoes')
orb_sat = propagador_orbital_mk3.propagador_orbital(rp, ecc, Raan, true_anomaly, inc, arg_per, num_orbita, 0)
print('Calculando a orientacao')
xyz = orientacao_quat(Ia, Ib, Ic, PSIP, TETAP, PHIP, psi0, teta0, phi0, T_orbita) #Ia, Ib, Ic, PSIP, TETAP, PHIP, Psi, Teta, Phi, Time_step
xyz['final'] = 1
xyz.to_csv('xyz.csv',sep=',')
'''K = len(xyz) / len(orb_sat)
xyz_novo = []

for i in range(0, len(xyz), int(K)):
    x_novo = (xyz[i][0])
    y_novo = (xyz[i][1])
    z_novo = (xyz[i][2])
    xyz_novo.append([(x_novo), (y_novo), (z_novo)])'''
xyz_novo = posi_ori(orb_sat, xyz)
xyz = pd.DataFrame(xyz_novo, columns=['Psi', 'Teta', 'Phi'])

Posicao_orientacao = pd.concat([orb_sat, xyz], axis=1)

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
            A = np.array([Posicao_orientacao.iloc[j, 6],
                          Posicao_orientacao.iloc[j, 7],
                          Posicao_orientacao.iloc[j, 8]])
            ra = Posicao_orientacao.iloc[i, 3]
            inc = Posicao_orientacao.iloc[i, 4]
            ome = Posicao_orientacao.iloc[i, 5]

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
            R.append([np.array(R1), np.array(R2), np.array(R3)])

      df2 = pd.DataFrame(R, columns=names[j])
      R = []
      Posicao_orientacao = pd.concat([Posicao_orientacao, df2], axis=1)
print('Calculando icosaedro')
nu = 10
vertices, faces = icosphere.icosphere(nu)
center = []
for i in range(0, len(faces), 1):
    A = vertices[faces[i][0]]
    B = vertices[faces[i][1]]
    C = vertices[faces[i][2]]
    x = A[0] + B[0] + C[0]
    y = A[1] + B[1] + C[1]
    z = A[2] + B[2] + C[2]
    center.append([x/3*Raio_terra, y/3*Raio_terra, z/3*Raio_terra])

As = area(vertices,faces, Raio_terra)
vet_terra = pd.DataFrame(center, columns=['Terra_X', 'Terra_Y', 'Terra_Z'])

Posicao_orientacao = pd.concat([Posicao_orientacao, vet_terra], axis=1)

Posicao_orientacao['final'] = 1
Posicao_orientacao.to_csv('posicao.csv',sep=',')

vetor_terra = []
for i in range(0, len(vet_terra), 1):
    vetor_terra.append([(vet_terra.iloc[i, 0]), (vet_terra.iloc[i, 1]), (vet_terra.iloc[i, 2])])
vetor_posicao = []
for i in range(0, len(orb_sat), 1):
    R = np.array([np.array(Posicao_orientacao.iloc[i][0]), np.array(Posicao_orientacao.iloc[i][1]), np.array(Posicao_orientacao.iloc[i][2])])
    R1 = R[0]
    R2 = R[1]
    R3 = R[2]
    vetor_posicao.append([R1, R2, R3])

'''Radiacao incidente do sol'''

print('Calculando radiacao solar')
Qs1 = []
Qs2 = []
Qs3 = []
Qs3 = []
Qs4 = []
Qs5 = []
Qs6 = []

for i in range(0, len(vetor_posicao), 1):

    PSI = np.arccos(np.dot(vetor_posicao[i]/np.linalg.norm(vetor_posicao[i]), Vs/np.linalg.norm(Vs)))
    QSI = np.arcsin(Raio_terra / np.linalg.norm((vetor_posicao[i])))

    if PSI + QSI < np.pi:

        A1 = np.array([np.array(Posicao_orientacao.iloc[i][9]), np.array(Posicao_orientacao.iloc[i][10]), np.array(Posicao_orientacao.iloc[i][11])])
        A2 = np.array([np.array(Posicao_orientacao.iloc[i][12]), np.array(Posicao_orientacao.iloc[i][13]), np.array(Posicao_orientacao.iloc[i][14])])
        A3 = np.array([np.array(Posicao_orientacao.iloc[i][15]), np.array(Posicao_orientacao.iloc[i][16]), np.array(Posicao_orientacao.iloc[i][17])])
        A4 = np.array([np.array(Posicao_orientacao.iloc[i][18]), np.array(Posicao_orientacao.iloc[i][19]), np.array(Posicao_orientacao.iloc[i][20])])
        A5 = np.array([np.array(Posicao_orientacao.iloc[i][21]), np.array(Posicao_orientacao.iloc[i][22]), np.array(Posicao_orientacao.iloc[i][23])])
        A6 = np.array([np.array(Posicao_orientacao.iloc[i][24]), np.array(Posicao_orientacao.iloc[i][25]), np.array(Posicao_orientacao.iloc[i][26])])

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


'''Radiacao de albedo incidente'''
print('Calculando radiacao de albedo')

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
i = 0
k = 0
A1 = np.array([(Posicao_orientacao.iloc[i, 9]), (Posicao_orientacao.iloc[i, 10]), (Posicao_orientacao.iloc[i, 11])])
A2 = np.array([(Posicao_orientacao.iloc[i, 12]), (Posicao_orientacao.iloc[i, 13]), (Posicao_orientacao.iloc[i, 14])])
A3 = np.array([(Posicao_orientacao.iloc[i, 15]), (Posicao_orientacao.iloc[i, 16]), (Posicao_orientacao.iloc[i, 17])])
A4 = np.array([(Posicao_orientacao.iloc[i, 18]), (Posicao_orientacao.iloc[i, 19]), (Posicao_orientacao.iloc[i, 20])])
A5 = np.array([(Posicao_orientacao.iloc[i, 21]), (Posicao_orientacao.iloc[i, 22]), (Posicao_orientacao.iloc[i, 23])])
A6 = np.array([(Posicao_orientacao.iloc[i, 24]), (Posicao_orientacao.iloc[i, 25]), (Posicao_orientacao.iloc[i, 26])])
rhok1 = np.array(vetor_posicao[i]) - np.array(vetor_terra[k]) + (A1)
rhok2 = np.array(vetor_posicao[i]) - np.array(vetor_terra[k]) + (A2)
rhok3 = np.array(vetor_posicao[i]) - np.array(vetor_terra[k]) + (A3)
rhok4 = np.array(vetor_posicao[i]) - np.array(vetor_terra[k]) + (A4)
rhok5 = np.array(vetor_posicao[i]) - np.array(vetor_terra[k]) + (A5)
rhok6 = np.array(vetor_posicao[i]) - np.array(vetor_terra[k]) + (A6)

print(rhok1)
print(np.dot(rhok1, vetor_terra[k]))
for i in range(0, len(vetor_posicao), 1):

    PSI = np.arccos(np.dot(vetor_posicao[i], Vs)/(np.linalg.norm(vetor_posicao[i])*np.linalg.norm(Vs)))
    QSI = np.arcsin(Raio_terra / np.linalg.norm(vetor_posicao[i]))

    if PSI < np.pi/2:


        for k in range(0, len(vetor_terra), 1):
            A1 = np.array(
                [(Posicao_orientacao.iloc[i, 9]), (Posicao_orientacao.iloc[i, 10]), (Posicao_orientacao.iloc[i, 11])])
            A2 = np.array(
                [(Posicao_orientacao.iloc[i, 12]), (Posicao_orientacao.iloc[i, 13]), (Posicao_orientacao.iloc[i, 14])])
            A3 = np.array(
                [(Posicao_orientacao.iloc[i, 15]), (Posicao_orientacao.iloc[i, 16]), (Posicao_orientacao.iloc[i, 17])])
            A4 = np.array(
                [(Posicao_orientacao.iloc[i, 18]), (Posicao_orientacao.iloc[i, 19]), (Posicao_orientacao.iloc[i, 20])])
            A5 = np.array(
                [(Posicao_orientacao.iloc[i, 21]), (Posicao_orientacao.iloc[i, 22]), (Posicao_orientacao.iloc[i, 23])])
            A6 = np.array(
                [(Posicao_orientacao.iloc[i, 24]), (Posicao_orientacao.iloc[i, 25]), (Posicao_orientacao.iloc[i, 26])])
            #rho =   np.array(vetor_posicao[i]) - np.array(vetor_terra[k]))

            rhok1 = np.array(vetor_posicao[i]) - np.array(vetor_terra[k]) + (A1)
            rhok2 = np.array(vetor_posicao[i]) - np.array(vetor_terra[k]) + (A2)
            rhok3 = np.array(vetor_posicao[i]) - np.array(vetor_terra[k]) + (A3)
            rhok4 = np.array(vetor_posicao[i]) - np.array(vetor_terra[k]) + (A4)
            rhok5 = np.array(vetor_posicao[i]) - np.array(vetor_terra[k]) + (A5)
            rhok6 = np.array(vetor_posicao[i]) - np.array(vetor_terra[k]) + (A6)


            #As = np.array([Posicao_orientacao.iloc[k][31]])
            C_bek = np.dot(Vs, vetor_terra[k])/(np.linalg.norm(Vs)*np.linalg.norm(vetor_terra[k]))

            if np.dot(rhok1, vetor_terra[k]) > 0:
                C_aek1 = np.dot(rhok1, vetor_terra[k]) / (np.linalg.norm(rhok1) * np.linalg.norm(vetor_terra[k]))
                C_aik1 = (np.dot(-rhok1, A1)) / (np.linalg.norm(rhok1) * np.linalg.norm(A1))

                if C_aik1 > 0 and C_aek1 > 0 and C_bek > 0:
                    Balb = float(1.0)
                    Halb1 = Halb1 + (As * ((C_aek1 * C_bek * C_aik1) / (np.pi * np.linalg.norm(rhok1) ** 2)) * Balb)
                else:
                    Balb = float(0.0)

                    Halb1 = Halb1 + (As * ((C_aek1 * C_bek * C_aik1) / (np.pi * np.linalg.norm(rhok1) ** 2)) * Balb)

            if np.dot(rhok2, vetor_terra[k]) > 0:

                C_aek2 = np.dot(rhok2, vetor_terra[k]) / (np.linalg.norm(rhok2) * np.linalg.norm(vetor_terra[k]))
                C_aik2 = (np.dot(-rhok2, A2)) / (np.linalg.norm(rhok2) * np.linalg.norm(A2))

                if C_aik2 > 0 and C_aek2 > 0 and C_bek > 0:

                    Balb = float(1.0)

                    Halb2 = Halb2 + (As * ((C_aek2 * C_bek * C_aik2) / (np.pi * np.linalg.norm(rhok2) ** 2)) * Balb)
                else:
                    Balb = float(0.0)

                    Halb2 = Halb2 + (As * ((C_aek2 * C_bek * C_aik2) / (np.pi * np.linalg.norm(rhok2) ** 2)) * Balb)

            if np.dot(rhok3, vetor_terra[k]) > 0:

                C_aek3 = np.dot(rhok3, vetor_terra[k]) / (np.linalg.norm(rhok3) * np.linalg.norm(vetor_terra[k]))
                C_aik3 = (np.dot(-rhok3, A3)) / (np.linalg.norm(rhok3) * np.linalg.norm(A3))

                if C_aik3 > 0 and C_aek3 > 0 and C_bek > 0:

                    Balb = float(1.0)

                    Halb3 = Halb3 + (As * ((C_aek3 * C_bek * C_aik3) / (np.pi * np.linalg.norm(rhok3) ** 2)) * Balb)

                else:
                    Balb = float(0.0)

                    Halb3 = Halb3 + (As * ((C_aek3 * C_bek * C_aik3) / (np.pi * np.linalg.norm(rhok3) ** 2)) * Balb)

            if np.dot(rhok4, vetor_terra[k]) > 0:

                C_aek4 = np.dot(rhok4, vetor_terra[k]) / (np.linalg.norm(rhok4) * np.linalg.norm(vetor_terra[k]))
                C_aik4 = (np.dot(-rhok4, A4)) / (np.linalg.norm(rhok4) * np.linalg.norm(A4))

                if C_aik4 > 0 and C_aek4 > 0 and C_bek > 0:

                    Balb = float(1.0)

                    Halb4 = Halb4 + (As * ((C_aek4 * C_bek * C_aik4) / (np.pi * np.linalg.norm(rhok4) ** 2)) * Balb)
                else:
                    Balb = float(0.0)

                    Halb4 = Halb4 + (As * ((C_aek4 * C_bek * C_aik4) / (np.pi * np.linalg.norm(rhok4) ** 2)) * Balb)

            if np.dot(rhok5, vetor_terra[k]) > 0:

                C_aek5 = np.dot(rhok5, vetor_terra[k]) / (np.linalg.norm(rhok5) * np.linalg.norm(vetor_terra[k]))
                C_aik5 = (np.dot(-rhok5, A5)) / (np.linalg.norm(rhok5) * np.linalg.norm(A5))

                if C_aik5 > 0 and C_aek5 > 0 and C_bek > 0:

                    Balb = float(1.0)

                    Halb5 = Halb5 + (As * ((C_aek5 * C_bek * C_aik5) / (np.pi * np.linalg.norm(rhok5) ** 2)) * Balb)
                else:
                    Balb = float(0.0)

                    Halb5 = Halb5 + (As * ((C_aek5 * C_bek * C_aik5) / (np.pi * np.linalg.norm(rhok5) ** 2)) * Balb)

            if np.dot(rhok6, vetor_terra[k]) > 0:

                C_aek6 = np.dot(rhok6, vetor_terra[k]) / (np.linalg.norm(rhok6) * np.linalg.norm(vetor_terra[k]))
                C_aik6 = (np.dot(-rhok6, A6)) / (np.linalg.norm(rhok6) * np.linalg.norm(A6))

                if C_aik6 > 0 and C_aek6 > 0 and C_bek > 0:

                    Balb = float(1.0)

                    Halb6 = Halb6 + (As * ((C_aek6 * C_bek * C_aik6) / (np.pi * np.linalg.norm(rhok6) ** 2)) * Balb)
                else:
                    Balb = float(0.0)

                    Halb6 = Halb6 + (As * ((C_aek6 * C_bek * C_aik6) / (np.pi * np.linalg.norm(rhok6) ** 2)) * Balb)

        Qalb1.append(ai * gama * Is * Halb1)
        Halb1 = 0

        Qalb2.append(ai * gama * Is * Halb2)
        Halb2 = 0

        Qalb3.append(ai * gama * Is * Halb3)
        Halb3 = 0

        Qalb4.append(ai * gama * Is * Halb4)
        Halb4 = 0

        Qalb5.append(ai * gama * Is * Halb5)
        Halb5 = 0

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
print('Calculando radiacao da terra')
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
i = 0

A1 = np.array([(Posicao_orientacao.iloc[i, 9]), (Posicao_orientacao.iloc[i, 10]), (Posicao_orientacao.iloc[i, 11])])
A2 = np.array([(Posicao_orientacao.iloc[i, 12]), (Posicao_orientacao.iloc[i, 13]), (Posicao_orientacao.iloc[i, 14])])
A3 = np.array([(Posicao_orientacao.iloc[i, 15]), (Posicao_orientacao.iloc[i, 16]), (Posicao_orientacao.iloc[i, 17])])
A4 = np.array([(Posicao_orientacao.iloc[i, 18]), (Posicao_orientacao.iloc[i, 19]), (Posicao_orientacao.iloc[i, 20])])
A5 = np.array([(Posicao_orientacao.iloc[i, 21]), (Posicao_orientacao.iloc[i, 22]), (Posicao_orientacao.iloc[i, 23])])
A6 = np.array([(Posicao_orientacao.iloc[i, 24]), (Posicao_orientacao.iloc[i, 25]), (Posicao_orientacao.iloc[i, 26])])
print(f'A1: {(A1)}')
print(f'A2: {(A2)}')
print(f'A3: {(A3)}')
print(f'A4: {(A4)}')
print(f'A5: {(A5)}')
print(f'A6: {(A6)}')


for i in range(0, len(vetor_posicao), 1):

    A1 = np.array([(Posicao_orientacao.iloc[i, 9]), (Posicao_orientacao.iloc[i, 10]), (Posicao_orientacao.iloc[i, 11])])
    A2 = np.array([(Posicao_orientacao.iloc[i, 12]), (Posicao_orientacao.iloc[i, 13]), (Posicao_orientacao.iloc[i, 14])])
    A3 = np.array([(Posicao_orientacao.iloc[i, 15]), (Posicao_orientacao.iloc[i, 16]), (Posicao_orientacao.iloc[i, 17])])
    A4 = np.array([(Posicao_orientacao.iloc[i, 18]), (Posicao_orientacao.iloc[i, 19]), (Posicao_orientacao.iloc[i, 20])])
    A5 = np.array([(Posicao_orientacao.iloc[i, 21]), (Posicao_orientacao.iloc[i, 22]), (Posicao_orientacao.iloc[i, 23])])
    A6 = np.array([(Posicao_orientacao.iloc[i, 24]), (Posicao_orientacao.iloc[i, 25]), (Posicao_orientacao.iloc[i, 26])])
    for k in range(0, len(vetor_terra), 1):

        rho = (np.array(vetor_posicao[i]) - np.array(vetor_terra[k]))
        Rhok1 = rho + np.array(A1)
        Rhok2 = rho + np.array(A2)
        Rhok3 = rho + np.array(A3)
        Rhok4 = rho + np.array(A4)
        Rhok5 = rho + np.array(A5)
        Rhok6 = rho + np.array(A6)

        if np.dot(Rhok1, vetor_terra[k]) > 0:

            C_aek1 = np.dot(Rhok1, vetor_terra[k]) / (np.linalg.norm(Rhok1) * np.linalg.norm(vetor_terra[k]))
            C_aik1 = (np.dot(-Rhok1, A1)) / (np.linalg.norm(Rhok1) * np.linalg.norm(A1))
            if C_aik1 > 0 and C_aek1 > 0:

                Balb = float(1.0)

                Hrad1 = Hrad1 + (As * ((C_aek1 * C_aik1) / (np.pi * np.linalg.norm(Rhok1) ** 2)) * Balb)
            else:
                Balb = float(0.0)

                Hrad1 = Hrad1 + (As * ((C_aek1 * C_aik1) / (np.pi * np.linalg.norm(Rhok1) ** 2)) * Balb)

        if np.dot(Rhok2, vetor_terra[k]) > 0:

            C_aek2 = np.dot(Rhok2, vetor_terra[k]) / (np.linalg.norm(Rhok2) * np.linalg.norm(vetor_terra[k]))
            C_aik2 = (np.dot(-Rhok2, A2)) / (np.linalg.norm(Rhok2) * np.linalg.norm(A2))

            if C_aik2 > 0 and C_aek2 > 0:

                Balb = float(1.0)

                Hrad2 = Hrad2 + (As * ((C_aek2 * C_aik2) / (np.pi * np.linalg.norm(Rhok2) ** 2)) * Balb)
            else:
                Balb = float(0.0)

                Hrad2 = Hrad2 + (As * ((C_aek2 * C_aik2) / (np.pi * np.linalg.norm(Rhok2) ** 2)) * Balb)

        if np.dot(Rhok3, vetor_terra[k]) > 0:

            C_aek3 = np.dot(Rhok3, vetor_terra[k]) / (np.linalg.norm(Rhok3) * np.linalg.norm(vetor_terra[k]))
            C_aik3 = (np.dot(-Rhok3, A3)) / (np.linalg.norm(Rhok3) * np.linalg.norm(A3))

            if C_aik3 > 0 and C_aek3 > 0:

                Balb = float(1.0)

                Hrad3 = Hrad3 + (As * ((C_aek3 * C_aik3) / (np.pi * np.linalg.norm(Rhok3) ** 2)) * Balb)
            else:
                Balb = float(0.0)

                Hrad3 = Hrad3 + (As * ((C_aek3 * C_aik3) / (np.pi * np.linalg.norm(Rhok3) ** 2)) * Balb)

        if np.dot(Rhok4, vetor_terra[k]) > 0:

            C_aek4 = np.dot(Rhok4, vetor_terra[k]) / (np.linalg.norm(Rhok4) * np.linalg.norm(vetor_terra[k]))
            C_aik4 = (np.dot(-Rhok4, A4)) / (np.linalg.norm(Rhok4) * np.linalg.norm(A4))

            if C_aik4 > 0 and C_aek4 > 0:

                Balb = float(1.0)

                Hrad4 = Hrad4 + (As * ((C_aek4 * C_aik4) / (np.pi * np.linalg.norm(Rhok4) ** 2)) * Balb)
            else:
                Balb = float(0.0)

                Hrad4 = Hrad4 + (As * ((C_aek4 * C_aik4) / (np.pi * np.linalg.norm(Rhok4) ** 2)) * Balb)

        if np.dot(Rhok5, vetor_terra[k]) > 0:

            C_aek5 = np.dot(Rhok5, vetor_terra[k]) / (np.linalg.norm(Rhok5) * np.linalg.norm(vetor_terra[k]))
            C_aik5 = (np.dot(-Rhok5, A5)) / (np.linalg.norm(Rhok5) * np.linalg.norm(A5))

            if C_aik5 > 0 and C_aek5 > 0:

                Balb = float(1.0)

                Hrad5 = Hrad5 + (As * ((C_aek5 * C_aik5) / (np.pi * np.linalg.norm(Rhok5) ** 2)) * Balb)
            else:
                Balb = float(0.0)

                Hrad5 = Hrad5 + (As * ((C_aek5 * C_aik5) / (np.pi * np.linalg.norm(Rhok5) ** 2)) * Balb)

        if np.dot(Rhok6, vetor_terra[k]) > 0:

            C_aek6 = np.dot(Rhok6, vetor_terra[k]) / (np.linalg.norm(Rhok6) * np.linalg.norm(vetor_terra[k]))
            C_aik6 = (np.dot(-Rhok6, A6)) / (np.linalg.norm(Rhok6) * np.linalg.norm(A6))

            if C_aik6 > 0 and C_aek6 > 0:

                Balb = float(1.0)

                Hrad6 = Hrad6 + (As * ((C_aek6 * C_aik6) / (np.pi * np.linalg.norm(Rhok6) ** 2)) * Balb)
            else:
                Balb = float(0.0)

                Hrad6 = Hrad6 + (As * ((C_aek6 * C_aik6) / (np.pi * np.linalg.norm(Rhok6) ** 2)) * Balb)

    Qrad1.append(e * Ir * (Hrad1))
    Hrad1 = 0

    Qrad2.append(e * Ir * (Hrad2))
    Hrad2 = 0

    Qrad3.append(e * Ir * (Hrad3))
    Hrad3 = 0

    Qrad4.append(e * Ir * (Hrad4))
    Hrad4 = 0

    Qrad5.append(e * Ir * (Hrad5))
    Hrad5 = 0

    Qrad6.append(e * Ir * (Hrad6))
    Hrad6 = 0

rad_sol = []
for i in range(0, len(Qs1), 1):
    rad_sol.append([Qs1[i], Qs2[i], Qs3[i], Qs4[i], Qs5[i], Qs6[i]])

Q_sol = pd.DataFrame(rad_sol, columns=['Qs1', 'Qs2', 'Qs3', 'Qs4', 'Qs5', 'Qs6'])

rad_alb = []
for i in range(0, len(Qalb1), 1):
    rad_alb.append([Qalb1[i], Qalb2[i], Qalb3[i], Qalb4[i], Qalb5[i], Qalb6[i]])
Q_alb = pd.DataFrame(rad_alb, columns=['Qalb1', 'Qalb2', 'Qalb3', 'Qalb4', 'Qalb5', 'Qalb6'])

rad_terra = []
for i in range(0, len(Qrad1), 1):
    rad_terra.append([Qrad1[i], Qrad2[i], Qrad3[i], Qrad4[i], Qrad5[i], Qrad6[i]])
Q_terra = pd.DataFrame(rad_terra, columns=['Qrad1', 'Qrad2', 'Qrad3', 'Qrad4', 'Qrad5', 'Qrad6'])


print('Terminando calculo')

QT = pd.concat([Q_sol, Q_alb], axis=1)
QT = pd.concat([QT, Q_terra], axis=1)
QT['N1'] = QT['Qs1'] + QT['Qalb1'] #+ QT['Qrad1']
QT['N2'] = QT['Qs2'] + QT['Qalb2'] #+ QT['Qrad2']
QT['N3'] = QT['Qs3'] + QT['Qalb3'] #+ QT['Qrad3']
QT['N4'] = QT['Qs4'] + QT['Qalb4'] #+ QT['Qrad4']
QT['N5'] = QT['Qs5'] + QT['Qalb5'] #+ QT['Qrad5']
QT['N6'] = QT['Qs6'] + QT['Qalb6'] #+ QT['Qrad6']
QT.to_csv('calor.csv',sep=',')
T = np.linspace(0, len(vetor_posicao), len(vetor_posicao))
plt.xlabel("Ponto da orbita")
plt.ylabel("Calor incidente em cada face [W/m^2]")
plt.plot(T, QT['N1'], color ='green', label='N1')
plt.plot(T, QT['N2'], color = 'blue', label='N2')
plt.plot(T, QT['N3'], color = 'cyan', label='N3')
plt.plot(T, QT['N4'], color = 'yellow', label='N4')
plt.plot(T, QT['N5'], color = 'red', label='N5')
plt.plot(T, QT['N6'], color = 'magenta', label='N6')
plt.legend()
plt.show()

'''fig = plt.figure()
ax = plt.axes(projection="3d")
a1 = ax.scatter3D(Posicao_orientacao['Terra_X'], Posicao_orientacao['Terra_Y'], Posicao_orientacao['Terra_Z'])
ax.set_xlabel('x [km]')
ax.set_ylabel('y [km]')
ax.set_zlabel('z [km]')
plt.show()'''