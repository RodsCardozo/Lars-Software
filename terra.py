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
    posicao = []
    x = []
    y = []
    z = []
    teta = np.linspace(0, 2 * np.pi, divisao)
    phi = np.linspace(0 + np.pi / divisao, np.pi - np.pi / divisao, divisao)
    for i in range(0, divisao, 1):
        for j in range(0, divisao, 1):
            x.append(Raio_terra * np.cos(teta[i]) * np.sin(phi[j]))
            y.append(Raio_terra * np.sin(teta[i]) * np.sin(phi[j]))
            z.append(Raio_terra * np.cos(phi[j]))
    for i in range(0, len(x), 1):
        posicao.append([x[i], y[i], z[i]])
    Terra = pd.DataFrame(posicao, columns=['Terra_X', 'Terra_Y', 'Terra_Z'])
    Terra['Terra_R'] = np.sqrt(Terra['Terra_X']**2 + Terra['Terra_Y']**2 + Terra['Terra_Z']**2)
    return Terra

