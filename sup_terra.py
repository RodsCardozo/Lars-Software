# -*- coding: utf-8 -*-
# sup_terra.py>
"""
    Universidade Federal de Santa Catarina
    Thermal Fluid Flow Group (T2F) - Aerospace Engineering Team
    Orbital Mechanics Division

    Título do Algoritmo: Divisao da terra em elementos de area
    Autor: Rodrigo S. Cardozo
    Versão: 1.0
    Data: 08/08/2022

    Divisao da terra em areas e vetores normais

"""


def sup_terra(Raio_terra, divisao):
    import numpy as np
    import pandas as pd
    Rt = Raio_terra
    St = []
    divi = divisao + 1
    alfa = np.linspace(0, 2 * np.pi, divi)
    beta = np.linspace(0, np.pi, divi)

    for i in range(0, divi - 1, 1):
        for j in range(0, divi - 1, 1):
            t1 = alfa[1] - alfa[0]
            t2 = beta[1] - beta[0]
            St.append(np.array((Rt ** 2) * (alfa[i] + t1 - alfa[i]) * (
                    np.cos(beta[j]) - (np.cos(beta[j]) * np.cos(t2) - np.sin(beta[j]) * np.sin(t2)))))
    Sup_terra = pd.DataFrame(St, columns=['As'])
    return Sup_terra
