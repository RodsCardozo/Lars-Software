# -*- coding: utf-8 -*-
# orientacao_quat.py>
"""
    Universidade Federal de Santa Catarina
    Thermal Fluid Flow Group (T2F) - Aerospace Engineering Team
    Orbital Mechanics Division

    Título do Algoritmo: Integrador de velocidade angular por quaternions
    Autor: Rodrigo S. Cardozo
    Versão: 0.1
    Data: 08/07/2022

"""


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
    import numpy as np
    import pandas as pd
    from scipy.integrate import odeint

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

    passo = 50000
    t = np.linspace(0, Time_step, passo)

    sol2 = odeint(quat, qi, t, args=(Ix3, Iy3, Iz3))  # integracao do conjunto de EDO
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

    return (PTP)