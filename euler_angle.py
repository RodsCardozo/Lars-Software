# -*- coding: utf-8 -*-
# euler_angle.py>
"""
    Universidade Federal de Santa Catarina
    Thermal Fluid Flow Group (T2F) - Aerospace Engineering Team
    Orbital Mechanics Division

    Título do Algoritmo: Transformador de coordenada por angulos de Euler 313
    Autor: Rodrigo S. Cardozo
    Versão: 1.0
    Data: 24/08/2022

"""


def rotacao_euler_2(VETOR, RA, INC, OME):
    import numpy as np

    """
    VETOR = Vetor a ser transforado
    RA = Angulo PSI
    INC = Angulo TETA
    OME = Angulo PHI

    """
    ra = (float(RA))
    inc = (float(INC))
    ome = (float(OME))
    vetor = np.array(VETOR)

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
    R = np.transpose(Tci).dot(vetor)

    x = R[0]
    y = R[1]
    z = R[2]
    Rxyz = np.array([x, y, z])

    return Rxyz
