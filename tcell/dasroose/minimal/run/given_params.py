'''
Created on Aug 13, 2014

@author: karmel

Parameters from the published supplement.
'''
from scipy.integrate import odeint

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from tcell.dasroose.minimal.constants import Rates, DiffusionBox
from tcell.dasroose.minimal.equations import Reactions
from tcell.dasroose.minimal.species import Ras, RasGAP, Sos
if __name__ == '__main__':
    box = DiffusionBox()

    # Rates are listed in uM, from table S1
    # k_off and kcat rates are already in s^-1; k_on needs to be converted.
    rates = [box.convert_membrane_rate(0.12), 3.0,
             box.convert_membrane_rate(0.11), 0.4,
             box.convert_membrane_rate(0.05), 0.1, 0.038,
             box.convert_membrane_rate(0.07), 1.0, 0.003,
             box.convert_cytosolic_rate(1.74), 0.2, 0.1]
    print(rates)
    rates = Rates(*rates)
    print(rates.K3m, rates.K4m, rates.K5m)
    # We have used a intial Ras-GDP, and Ras-GAP concentrations of
    # 75 molecules/(um)2, 125 molecules/(um)3 respectively.
    ras = Ras(GDP=75)
    rasgap = RasGAP(125)

    # Do the simplest thing first. Solve the ODE ignoring SOS-RasGTP
    sos = Sos(num_molecules=100)
    reactions = Reactions(rates, sos.num_molecules,
                          ras.num_molecules, rasgap.num_molecules)

    # y = [Sos, RasGDP]. Solve dy/dt
    def dy_dt(y, t):
        Sos_i = y[0]
        Sos_GTP_i = y[1]
        RasGTP_i = y[2]

        f0 = reactions.dSos_dt(Sos_i, Sos_GTP_i, RasGTP_i)
        f1 = reactions.dSos_GTP_dt(Sos_i, Sos_GTP_i, RasGTP_i)
        f2 = reactions.dRasGTP_dt(Sos_i, Sos_GTP_i, RasGTP_i)
        # print(f0,f1,f2)
        return [f0, f1, f2]

    # Initial conditions
    y0 = [sos.unbound, sos.GTP, ras.GTP]
    t = np.linspace(0, 600., 10000)

    solution = odeint(dy_dt, y0, t)
    Sos_unbound = solution[:, 0]
    Sos_GTP = solution[:, 1]
    RasGTP = solution[:, 2]

    plt.figure()
    plt.plot(t, Sos_unbound, label='SOS')
    plt.plot(t, Sos_GTP, label='SOS-RasGTP')
    plt.plot(t, RasGTP, label='RasGTP')
    plt.legend()
    plt.xlabel('Time')
    plt.ylabel('Molecule count')
    plt.title('Das Minimal Model: Species presence (SOS initial: {})'.format(
        sos.num_molecules))
    plt.show()
