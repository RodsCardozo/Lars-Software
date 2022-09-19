# -*- coding: utf-8 -*-
# periodo_orbital.py>

"""
    Universidade Federal de Santa Catarina
    Thermal Fluid Flow Group (T2F) - Aerospace Engineering Team
    Orbital Mechanics Division

    Título do Algoritmo: Integrador de velocidade angular por quaternions
    Autor: Rodrigo S. Cardozo
    Versão: 0.1
    Data: 08/07/2022
"""


def periodo_orbital(Perigeu):
    """
    Perigeu = Altitude Inicial do Satelite no perigeu
    """
    import numpy as np
    mu = 398600
    T_orb = float(((2 * np.pi) / (np.sqrt(mu))) * (Perigeu ** (3 / 2)))
    return (T_orb)