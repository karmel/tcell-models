'''
Created on Aug 13, 2014

@author: karmel

Parameters from the published supplement.
'''
import os

from scipy.integrate import odeint
from scipy.optimize.minpack import fsolve

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from tcell.dasroose.minimal.constants import Rates, DiffusionBox
from tcell.dasroose.minimal.equations import Reactions
from tcell.dasroose.minimal.species import Ras, RasGAP, Sos


def plot_molecule_count(rates, ras, rasgap):
    for sos_mol in range(50, 400, 50):

        sos = Sos(num_molecules=sos_mol)
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
            return [f0, f1, f2]

        # Initial conditions
        y0 = [sos.unbound, sos.GTP, ras.GTP]
        t = np.linspace(0, 600., 10000)

        #################################################
        # Molecule count by time
        #################################################
        solution = odeint(dy_dt, y0, t)
        Sos_unbound = solution[:, 0]
        Sos_GTP = solution[:, 1]
        RasGTP = solution[:, 2]

        plt.figure(figsize=[10, 8])
        plt.plot(t, Sos_unbound, label='SOS')
        plt.plot(t, Sos_GTP, label='SOS-RasGTP')
        plt.plot(t, RasGTP, label='RasGTP')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Molecule count')

        title = 'Das Minimal Model- Species presence (SOS initial={})'.format(
            sos.num_molecules)
        dirname = os.path.join(os.path.dirname(
            os.path.realpath(__file__)), 'output')
        plt.title(title)
        plt.savefig(os.path.join(dirname, title.replace(' ', '-') + '.png'))
        # plt.show()

        #################################################
        # Rate-balance plots
        #################################################

        # Chop off the first 1000 timepoints because it completely throws off
        # the scale.
        skip = 1000
        calculated_rates = np.array([dy_dt(y, 0) for y in solution[skip:]])
        dSos_dt = calculated_rates[:, 0]
        dSos_GTP_dt = calculated_rates[:, 1]
        dRasGTP_dt = calculated_rates[:, 2]
        trunc_t = t[skip:]
        
        # Rate vs Time
        plt.figure(figsize=[10, 8])
        plt.plot(trunc_t, dSos_dt, label='dSOS/dt')
        plt.plot(trunc_t, dSos_GTP_dt, label='dSOS-RasGTP/dt')
        plt.plot(trunc_t, dRasGTP_dt, label='dRasGTP/dt')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Rate')

        title = 'Das Minimal Model- Rates of change (SOS initial={})'.format(
            sos.num_molecules)
        dirname = os.path.join(os.path.dirname(
            os.path.realpath(__file__)), 'output')
        plt.title(title)
        plt.savefig(os.path.join(dirname, title.replace(' ', '-') + '.png'))
        # plt.show()

        # Rate vs RasGTP
        trunc_RasGTP = RasGTP[skip:]
        plt.figure(figsize=[10, 8])
        plt.plot(trunc_RasGTP, dSos_dt, label='dSOS/dt')
        plt.plot(trunc_RasGTP, dSos_GTP_dt, label='dSOS-RasGTP/dt')
        plt.plot(trunc_RasGTP, dRasGTP_dt, label='dRasGTP/dt')
        plt.legend()
        plt.xlabel('RasGTP Molecule Count')
        plt.ylabel('Rate')

        title = 'Das Minimal Model- Rate balance (SOS initial={})'.format(
            sos.num_molecules)
        dirname = os.path.join(os.path.dirname(
            os.path.realpath(__file__)), 'output')
        plt.title(title)
        plt.savefig(os.path.join(dirname, title.replace(' ', '-') + '.png'))
        # plt.show()


def find_steady_states(rates, ras, rasgap):
    ############################################
    # Solve for steady states
    ############################################
    sos = Sos(num_molecules=100)
    reactions = Reactions(rates, sos.num_molecules,
                          ras.num_molecules, rasgap.num_molecules)

    def system_of_eq(y):
        Sos_i = y[0]
        Sos_GTP_i = y[1]
        RasGTP_i = y[2]

        f0 = reactions.dSos_dt(Sos_i, Sos_GTP_i, RasGTP_i)
        f1 = reactions.dSos_GTP_dt(Sos_i, Sos_GTP_i, RasGTP_i)
        f2 = reactions.dRasGTP_dt(Sos_i, Sos_GTP_i, RasGTP_i)

        return [f0, f1, f2]
    print(system_of_eq([30,80,200]))
    Sos_i, Sos_GTP_i, RasGTP_i = fsolve(system_of_eq, [30, 80, 180])
    print(Sos_i, Sos_GTP_i, RasGTP_i)
    print(system_of_eq((Sos_i, Sos_GTP_i, RasGTP_i)))

if __name__ == '__main__':
    box = DiffusionBox()

    # Rates are listed in uM, from table S1
    # k_off and kcat rates are already in s^-1; k_on needs to be converted.
    rates = [box.convert_membrane_rate(0.12), 3.0,
             box.convert_membrane_rate(0.11), 0.4,
             box.convert_membrane_rate(0.05), 0.1, 0.038,
             box.convert_membrane_rate(0.07), 1.0, 0.003,
             box.convert_cytosolic_rate(1.74), 0.2, 0.1]
    rates = Rates(*rates)

    # We have used a intial Ras-GDP, and Ras-GAP concentrations of
    # 75 molecules/(um)2, 125 molecules/(um)3 respectively.
    ras = Ras(GDP=75)
    rasgap = RasGAP(125)

    plot_molecule_count(rates, ras, rasgap)
    #find_steady_states(rates, ras, rasgap)
