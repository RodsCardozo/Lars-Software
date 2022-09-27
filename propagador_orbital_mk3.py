# -*- coding: utf-8 -*-
# propagador_orbital_mk3.py>

"""
    Universidade Federal de Santa Catarina
    Thermal Fluid Flow Group (T2F) - Aerospace Engineering Team
    Orbital Mechanics Division
    
    Título do Algoritmo: Propagador orbital a partir dos parametros de Kepler
    Autor: Rodrigo S. Cardozo
    Versão: 0.3
    Data: 23/08/2022
    
    Simulador Orbital para determincao da orientacao do satelite em orbita
    
    definicao da quantidade de radiacao incidente sobre cada face do cubesat
    
"""


def propagador_orbital(Semi_eixo, excentricidade, asce_direita, anomalia_verdadeira, inclinacao, argu_perigeu,
                       num_orbitas, plotagem):
    """
    & Semi_eixo = altitude no periapse da orbita
    & excentricidade = e
    & asce_direita = Angulo da posicao do nodo ascendente
    & anomalia_verdadeira = algulo do vetor posicao e a linha dos apses com origem no foco
    & inclinacao = inclinacao da orbita
    & argu_perigeu = Angulo da orientacao da linha dos apses
    & num_orbitas = numero de orbitas a serem simuladas
    & plotagem = 1 para plotar / 0 para nao plotar a orbita
    
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import euler_angle
    import pandas as pd

    rp = float(Semi_eixo)  # semi eixo maior
    ecc = float(excentricidade)  # ecentricidade da orbita
    Raan = np.radians(float(asce_direita))  # ascencao direita do nodo ascendente
    arg_per = np.radians(float(argu_perigeu))  # argumento do perigeu
    true_anomaly = np.radians(float(anomalia_verdadeira/2))  # anomalia verdadeira
    inc = np.radians(float(inclinacao))  # inclinacao
    mu = 398600
    J2 = 1.08263e-3
    R_terra = 6371
    num_orbita = num_orbitas

    h = float(np.sqrt(rp * mu * (1 + ecc)))  # conserv da momento linear
    T_orb = float(((2 * np.pi) / (np.sqrt(mu))) * (rp ** (3 / 2)))  # periodo de uma orbita

    alfa = float((h ** 2 / mu) * (1 / (1 + ecc * np.cos(true_anomaly))))  # coef para calcular a posicao da orbita
    beta = float(mu / h)  # coef para calcular a velocidade na orbita

    r_posi_0 = [alfa * np.cos(true_anomaly), alfa * np.sin(true_anomaly), 0]  # posicao inicial
    v_velo_0 = [beta * (-np.sin(true_anomaly)), beta * (ecc + np.cos(true_anomaly)), 0]  # velocidade inicial

    r0 = np.sqrt(np.dot(r_posi_0, r_posi_0))  # modulo posicao inicial
    v0 = np.dot(r_posi_0, v_velo_0) / r0  # modulo velocidade inicial

    # listas com a posicao e velocidades ao longo da orbita
    r_posi = []
    v_velo = []
    r_posix = []
    r_posiy = []
    r_posiz = []

    # variacao da anomalia verdadeira
    inicio = true_anomaly
    fim = (2 * np.pi + true_anomaly) * num_orbita
    passo = 50 #int(T_orb*num_orbita)
    teta = np.linspace(inicio, fim, passo)

    ''' Posicao do satelite no plano perifocal '''

    for i in range(0, len(teta), 1):

        r = float(
            (h ** 2 / mu) * (1 / (1 + (h ** 2 / (mu * r0) - 1) * np.cos(teta[i]) - ((h * v0) / mu) * np.sin(teta[i]))))

        f = float(1 - ((mu * r) / h ** 2) * (1 - np.cos(teta[i])))
        g = float(((r0 * r) / h) * np.sin(teta[i]))

        if teta[i] == 0:
            fp = 0
        else:
            fp = float((mu / h) * ((1 - np.cos(teta[i])) / np.sin(teta[i])) * ((mu / h ** 2) *
                                                                               (1 - np.cos(teta[i])) - 1 / r0 - 1 / r))

        gp = float(1 - ((mu * r0) / h ** 2) * (1 - np.cos(teta[i])))

        r_posi.append(np.dot(f, r_posi_0) + np.dot(g, v_velo_0))
        v_velo.append(np.dot(fp, r_posi_0) + np.dot(gp, v_velo_0))

        r_posix.append(f * r_posi_0[0] + g * v_velo_0[0])
        r_posiy.append(f * r_posi_0[1] + g * v_velo_0[1])
        r_posiz.append(f * r_posi_0[2] + g * v_velo_0[2])
        # print(f*gp - fp*g)

    VET_posi = []
    r_X = []
    r_Y = []
    r_Z = []

    ''' Posicao do satelite no plano Inercial '''
    a = (h ** 2 / mu) * (1 / (1 - ecc ** 2))
    RA_p = 0.0 - (3 / 2) * ((np.sqrt(mu) * J2 * R_terra ** 2) / ((1 - ecc ** 2) ** 2 * a ** (7 / 2))) * np.cos(inc)
    arg_per_p = 0.0 - (3 / 2) * ((np.sqrt(mu) * J2 * R_terra ** 2) / ((1 - ecc ** 2) ** 2 * a ** (7 / 2))) * (
                (5 / 2) * np.sin(inc) ** 2 - 2)
    orbita = []
    for i in range(0, len(teta), 1):
        vet_posi = [r_posix[i], r_posiy[i], r_posiz[i]]
        VET_posi.append(euler_angle.rotacao_euler_2(vet_posi, Raan, inc, arg_per))
        r_X.append(VET_posi[i][0])
        r_Y.append(VET_posi[i][1])
        r_Z.append(VET_posi[i][2])

        E = 2 * np.arctan(np.sqrt((1 - ecc) / (1 + ecc)) * np.tan(teta[i]))
        n = (mu ** 2 / h ** 3) * (1 - ecc ** 2) ** (3 / 2)
        del_t = (E - E * np.sin(E)) / (n)
        #print(del_t )
        #Raan = Raan + RA_p * del_t
        #arg_per = arg_per + arg_per_p * del_t
        orbita.append([r_X[i], r_Y[i], r_Z[i]])
    if plotagem == 1:
        from mpl_toolkits.mplot3d import Axes3D

        fig = plt.figure()
        ax = Axes3D(fig)
        ax.set_xlim3d(min(r_X), max(r_X))
        ax.set_ylim3d(min(r_Y), max(r_Y))
        ax.set_zlim3d(min(r_Z), max(r_Z))
        ax.scatter3D(r_X, r_Y, r_Z)
        plt.show()

    else:
        print("Para plotar a orbita rode novamente com a opçao plotagem = 1")

    posi_xyz = pd.DataFrame(orbita, columns=['orb_X', 'orb_Y', 'orb_Z'])
    return (posi_xyz)
