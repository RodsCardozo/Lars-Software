# -*- coding: utf-8 -*-
# propagador_orbital_mk4.py>
'''
Codigo para integracao das equacoes do movimento do satelite
incluindo a atitude dele. Resolver o sistema de equacoes diferenciais utilizando
o solver odint. *ver como essa biblioteca faz a integracao.
'''
def propagador_orbital(data, semi_eixo, excentricidade, Raan, argumento_perigeu, anomalia_verdadeira, inclinacao, num_orbitas, delt, psi, teta, phi, psip, tetap, phip):

    """
    & Semi_eixo = altitude no periapse da orbita
    & excentricidade = e
    & Raan= Angulo da posicao do nodo ascendente
    & argumento_perigeu = Angulo da orientacao da linha dos apses
    & anomalia_verdadeira = algulo do vetor posicao e a linha dos apses com origem no foco
    & inclinacao = inclinacao da orbita
    & num_orbitas = numero de orbitas a serem simuladas
    & psi = primeiro angulo de Euler
    & teta = segundo angulo de Euler
    & psi = terceiro angulo de Euler
    & psip = velocidade angular do primeiro angulo de Euler
    & tetap = velocidade angular do segundo angulo de Euler
    & phip = velocidade angular do terceiro angulo de Euler

    """
    import numpy as np
    import pandas as pd
    from scipy.integrate import odeint
    from datetime import datetime
    from datetime import timedelta
    from nrlmsise00 import msise_model
    from periodo_orbital import periodo_orbital
    def propagador(q, t, Rho, velocidade, massa, CD, posicao, Area_transversal):  # funcao para integrar
        import numpy as np
        mu = 398600
        J2 = 1.08263e-3
        R_terra = 6371
        rho = Rho
        r = posicao
        m = float(massa)  # massa do cubesat
        a = float(0.1)  # comprimento do sat
        b = float(0.1)  # largura do sat
        c = float(0.2)  # altura do sat
        Ix3 = (m / 12) * (b ** 2 + c ** 2)  # momento de inercia na direcao x
        Iy3 = (m / 12) * (a ** 2 + c ** 2)  # momento de inercia na direcao y
        Iz3 = (m / 12) * (a ** 2 + b ** 2)  # momento de inercia na direcao z
        h, ecc, anomalia_verdadeira, raan, inc, arg_per, q0, q1, q2, q3, wx3, wy3, wz3 = q

        dMdt = [r*((-1/(2*r))*h*rho*velocidade*((CD*Area_transversal)/m) - 1.5 * ((J2*mu*R_terra**2)/r**4) * np.sin(inc)**2*np.sin(2*(arg_per + anomalia_verdadeira))),

                (h/mu)*np.sin(anomalia_verdadeira)*((-1/(2*h))*mu*ecc*rho*velocidade*((CD*Area_transversal)/m)*np.sin(anomalia_verdadeira)
                - 1.5*((J2*mu*R_terra**2)/r**4)*(1 - 3*np.sin(inc)**2*np.sin(arg_per + anomalia_verdadeira)**2))
                + (((-1/(2*r))*h*rho*velocidade*((CD*Area_transversal)/m) - 1.5*((J2*mu*R_terra**2)/r**4)*np.sin(inc)**2*np.sin(2*(arg_per
                + anomalia_verdadeira)))/(mu*h))*((h**2 + mu*r)*np.cos(anomalia_verdadeira) + mu*ecc*r),

                (h/r**2 + ((h**2*np.cos(anomalia_verdadeira))/(mu*ecc*h))*((-1/(2*h))*mu*ecc*rho*velocidade*((CD*Area_transversal)/m)*np.sin(anomalia_verdadeira)
                - 1.5*((J2*mu*R_terra**2)/r**4)*(1 - 3*np.sin(inc)**2*np.sin(arg_per + anomalia_verdadeira)**2))
                - (r + h**2/mu)*(np.sin(anomalia_verdadeira)/(ecc*h))*((-1/(2*r))*h*rho*velocidade*((CD*Area_transversal)/m)
                - (1.5 * (J2*mu*R_terra**2)/r**4) * np.sin(inc)**2 * np.sin(2*(arg_per + anomalia_verdadeira)))),

                (r/(h*np.sin(inc)))*np.sin(arg_per + anomalia_verdadeira)*(- 1.5*((J2*mu*R_terra**2)/r**4)*np.sin(2*inc)*np.sin(arg_per + anomalia_verdadeira)),

                (r / (h)) * np.cos(arg_per + anomalia_verdadeira) * (- 1.5 * ((J2 * mu * R_terra ** 2) / r ** 4) * np.sin(2 * inc) * np.sin(arg_per + anomalia_verdadeira)),

                (-1/(ecc*h))*((h**2/mu)*np.cos(anomalia_verdadeira)*((-1/(2*h))*mu*ecc*rho*velocidade*((CD*Area_transversal)/m)*np.sin(anomalia_verdadeira)
                - 1.5*((J2*mu*R_terra**2)/r**4)*(1 - 3*np.sin(inc)**2*np.sin(arg_per + anomalia_verdadeira)**2))
                - (r + h**2/mu)*np.sin(anomalia_verdadeira)*((-1/(2*r))*h*rho*velocidade*((CD*Area_transversal)/m)
                - 1.5 * ((J2*mu*R_terra**2)/r**4) * np.sin(inc)**2 * np.sin(2*(arg_per + anomalia_verdadeira))))
                - ((r*np.sin(arg_per + anomalia_verdadeira))/(h*np.tan(inc)))*(- 1.5 * (J2*mu*R_terra**2)/r**4 * np.sin(2*inc) * np.sin(arg_per + anomalia_verdadeira)),

                0.5 * (q0 * 0 - q1 * wx3 - q2 * wy3 - q3 * wz3),
                0.5 * (q1 * 0 + q0 * wx3 - q3 * wy3 + q2 * wz3),
                0.5 * (q2 * 0 + q3 * wx3 + q0 * wy3 - q1 * wz3),
                0.5 * (q3 * 0 - q2 * wx3 + q1 * wy3 + q0 * wz3),
                ((Iy3 - Iz3) / Ix3) * wy3 * wz3,
                ((Iz3 - Ix3) / Iy3) * wz3 * wx3,
                ((Ix3 - Iy3) / Iz3) * wx3 * wy3]
        return dMdt

    def rho(data, altitude, latitude, longitude):
        densidade = msise_model(data, altitude, latitude, longitude, 150, 150, 4, lst=16)
        rho = densidade[0][5] * 1000
        return rho

    # condicoes iniciais

    Rp = float(semi_eixo)  # semi eixo maior
    rp0 = Rp
    ecc0 = float(excentricidade)  # ecentricidade da orbita
    Raan0 = np.radians(float(Raan))  # ascencao direita do nodo ascendente
    arg_per0 = np.radians(float(argumento_perigeu))  # argumento do perigeu
    true_anomaly0 = np.radians(float(anomalia_verdadeira))  # anomalia verdadeira
    inc0 = np.radians(float(inclinacao))  # inclinacao
    T_orb = periodo_orbital(Rp)
    mu = 398600
    J2 = 1.08263e-3
    R_terra = 6371
    h0 = np.sqrt(rp0*mu*(1 - ecc0))
    psi = float(np.radians(psi))  # angulo inicial de PSI
    teta = float(np.radians(teta))  # angulo inicial de TETA
    phi = float(np.radians(phi))  # angulo inicial de PHI
    psip = float(psip)  # velocidade angular do angulo PSI
    tetap = float(tetap)  # velocidade angular do angulo TETA
    phip = phip  # velocidade angular do angulo PHI
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

    # comeco da integracao

    mu = 398600
    J2 = 1.08263e-3
    R_terra = 6371
    Time_step = delt
    passo = 10000
    ini_date = data
    n = num_orbitas
    T = T_orb*n
    t = np.linspace(0, Time_step, passo)
    solution = []
    posicao = (h0**2/mu)*(1/(1-ecc0*np.cos(true_anomaly0)))
    while delt < T:
        qi = [h0, ecc0, true_anomaly0, Raan0, inc0, arg_per0, q0, q1, q2, q3, wx3_i, wy3_i, wz3_i]
        altitude = rp0 - R_terra
        Rho = rho(ini_date, altitude, 0, 0)
        velocidade = (mu/h0)*np.sqrt(np.sin(true_anomaly0)**2 + (ecc0 + np.cos(true_anomaly0))**2)
        massa = 3.0
        CD = 2.2
        Area_transversal = 0.1*0.1
        sol = odeint(propagador, qi, t, args=(Rho, velocidade, massa, CD, posicao, Area_transversal))
        solution.append(sol[passo - 1])
        h0 = sol[passo-1][0]
        ecc0 = sol[passo-1][1]
        true_anomaly0 = sol[passo-1][2]
        Raan0 = sol[passo-1][3]
        inc0 = sol[passo-1][4]
        arg_per0 = sol[passo-1][5]
        posicao = (h0**2/mu)*(1/(1-ecc0*np.cos(true_anomaly0)))
        q0 = sol[passo-1][6]
        q1 = sol[passo-1][7]
        q2 = sol[passo-1][8]
        q3 = sol[passo-1][9]
        wx3_i = sol[passo-1][10]
        wy3_i = sol[passo-1][11]
        wz3_i = sol[passo-1][12]
        delt = delt + Time_step
        final_date = timedelta(seconds=Time_step)
        ini_date = ini_date + final_date

    solucao = pd.DataFrame(solution, columns=['h', 'ecc', 'anomalia_verdadeira', 'raan', 'inc', 'arg_per', 'q0', 'q1', 'q2', 'q3', 'wx3', 'wy3', 'wz3'])
    solucao['X_perifocal'] = (solucao['h']**2/mu)*(1/(1 + solucao['ecc']*np.cos(solucao['anomalia_verdadeira'])))*np.cos(solucao['anomalia_verdadeira'])
    solucao['Y_perifocal'] = (solucao['h']**2/mu)*(1/(1 + solucao['ecc']*np.cos(solucao['anomalia_verdadeira'])))*np.sin(solucao['anomalia_verdadeira'])
    solucao['Z_perifocal'] = 0
    solucao['distancia'] = np.sqrt(solucao['X_perifocal']**2 + solucao['Y_perifocal']**2)

    df = pd.DataFrame()
    df['PHI'] = np.unwrap(np.arctan2(2*(solucao['q2']*solucao['q3'] - solucao['q0']*solucao['q1']), 2*(solucao['q1']*solucao['q3'] + solucao['q0']*solucao['q2'])))
    df['TETA'] = np.arccos(2*(solucao['q0']**2 + solucao['q3']**2) - 1)
    df['PSI'] = np.unwrap(np.arctan2(-2*(solucao['q2']*solucao['q3'] + solucao['q0']*solucao['q1']), 2*(solucao['q1']*solucao['q3'] - solucao['q0']*solucao['q2']) ))

    df['X_ECI'] = ((np.cos(solucao['raan'])*np.cos(solucao['arg_per']) - np.sin(solucao['raan'])*np.sin(solucao['arg_per'])*np.cos(solucao['inc']))*solucao['X_perifocal']

                        + (-np.cos(solucao['raan'])*np.sin(solucao['arg_per']) - np.sin(solucao['raan'])*np.cos(solucao['inc'])*np.cos(solucao['arg_per']))*solucao['Y_perifocal']

                        + np.sin(solucao['raan'])*np.sin(solucao['inc'])*solucao['Z_perifocal'])

    df['Y_ECI'] = ((np.sin(solucao['raan'])*np.cos(solucao['arg_per']) + np.cos(solucao['raan'])*np.cos(solucao['inc'])*np.sin(solucao['arg_per']))*solucao['X_perifocal']

                        + (-np.sin(solucao['raan'])*np.sin(solucao['arg_per']) + np.cos(solucao['raan'])*np.cos(solucao['inc'])*np.cos(solucao['arg_per']))*solucao['Y_perifocal']

                        - np.cos(solucao['raan'])*np.sin(solucao['inc'])*solucao['Z_perifocal'])

    df['Z_ECI'] = (np.sin(solucao['inc'])*np.sin(solucao['arg_per'])*solucao['X_perifocal']
                        + np.sin(solucao['inc'])*np.cos(solucao['arg_per'])*solucao['Y_perifocal']
                        + np.cos(solucao['inc'])*solucao['Z_perifocal'])

    '''import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    a1 = ax.scatter3D(solucao['X_perifocal'], solucao['Y_perifocal'], solucao['Z_perifocal'])
    ax.set_xlabel('x [km]')
    ax.set_ylabel('y [km]')
    ax.set_zlabel('z [km]')
    plt.show()
    solucao['final'] = 1'''
    solucao.to_csv('solver.csv',sep=',')
    return df