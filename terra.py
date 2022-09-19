# -*- coding: utf-8 -*-

# terra.py>
"""
    Universidade Federal de Santa Catarina
    Thermal Fluid Flow Group (T2F) - Aerospace Engineering Team
    Orbital Mechanics Division

    Título do Algoritmo: Vetores normais da terra
    Autor: Rodrigo S. Cardozo
    Versão: 1.0
    Data: 08/08/2022

    Divisao da terra em areas e vetores normais

"""
def terra(Raio_terra, divisao):
    import numpy as np
    import pandas as pd
    x = []
    y = []
    z = []
    R = []

    teta = np.linspace(0, 2 * np.pi, divisao)
    phi = np.linspace(0 + np.pi / divisao, np.pi - np.pi / divisao, divisao)
    for i in range(0, divisao, 1):
        for j in range(0, divisao, 1):
            x.append(Raio_terra * np.cos(teta[i]) * np.sin(phi[j]))
            y.append(Raio_terra * np.sin(teta[i]) * np.sin(phi[j]))
            z.append(Raio_terra * np.cos(phi[j]))
        R.append(np.sqrt(x[i]**2 + y[i]**2 + z[i]**2))
    posicao = []
    for i in range(0, len(x), 1):
        posicao.append(np.array([x, y, z]))
    Terra = pd.DataFrame(posicao, columns=['Terra_X', 'Terra_Y', 'Terra_Z'])
    Terra['Terra_R'] = (Terra['Terra_X']**2 + Terra['Terra_Y']**2 + Terra['Terra_Z']**2)**0.5
    return Terra

