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
from tcell.dasroose.minimal.species import Ras, RasGAP, Sos, RasGRP


def dy_dt(y, t, reactions):
    '''
    y = [Sos unbound, Sos-RasGTP, RasGTP]. Solve dy/dt
    '''
    Sos_i = y[0]
    Sos_GTP_i = y[1]
    RasGTP_i = y[2]

    f0 = reactions.dSos_dt(Sos_i, Sos_GTP_i, RasGTP_i)
    f1 = reactions.dSos_GTP_dt(Sos_i, Sos_GTP_i, RasGTP_i)
    f2 = reactions.dRasGTP_dt(Sos_i, Sos_GTP_i, RasGTP_i)
    return [f0, f1, f2]


def solve_ode(sos, ras, rasgap, rasgrp):
    reactions = Reactions(rates, sos.num_molecules,
                          ras.num_molecules, rasgap.num_molecules,
                          rasgrp.num_molecules)

    # Initial conditions
    y0 = [sos.unbound, sos.GTP, ras.GTP]
    t = np.linspace(0, 600., 10000)

    solution = odeint(dy_dt, y0, t, args=(reactions,))
    return t, reactions, solution


def plot_molecule_count(rates, sos_counts, ras, rasgap, rasgrp,
                        output_dirname='output'):
    final_vals = {}
    for sos_mol in sos_counts:

        sos = Sos(num_molecules=sos_mol)

        t, reactions, solution = solve_ode(sos, ras, rasgap, rasgrp)

        #################################################
        # Molecule count by time
        #################################################
        final_vals[sos_mol] = solution[-1:][0]
        print('Final counts for Sos initial = {}: {}'.format(
            sos_mol, final_vals[sos_mol]))
        Sos_unbound = solution[:, 0]
        Sos_GTP = solution[:, 1]
        RasGTP = solution[:, 2]

        dirname = os.path.join(output_dirname, 'molecule_count')
        if not os.path.exists(dirname):
            os.mkdir(dirname)

        plt.figure(figsize=[10, 8])
        plt.plot(t, Sos_unbound, label='SOS')
        plt.plot(t, Sos_GTP, label='SOS-RasGTP')
        plt.plot(t, RasGTP, label='RasGTP')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Molecule count')

        title = 'Das Minimal Model- Species presence (SOS initial={})'.format(
            sos.num_molecules)
        plt.title(title)
        plt.savefig(os.path.join(dirname, title.replace(' ', '-') + '.png'))
        # plt.show()
        plt.close()

        #################################################
        # Rate-balance plots
        #################################################

        # Chop off the first timepoints because it completely throws off
        # the scale.
        skip = 100
        calculated_rates = [dy_dt(y, 0, reactions) for y in solution[skip:]]
        calculated_rates = np.array(calculated_rates)
        dSos_dt = calculated_rates[:, 0]
        dSos_GTP_dt = calculated_rates[:, 1]
        dRasGTP_dt = calculated_rates[:, 2]
        trunc_t = t[skip:]

        # Rate vs Time
        dirname = os.path.join(output_dirname, 'rates_of_change')
        if not os.path.exists(dirname):
            os.mkdir(dirname)

        plt.figure(figsize=[10, 8])
        plt.plot(trunc_t, dSos_dt, label='dSOS/dt')
        plt.plot(trunc_t, dSos_GTP_dt, label='dSOS-RasGTP/dt')
        plt.plot(trunc_t, dRasGTP_dt, label='dRasGTP/dt')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Rate')

        title = 'Das Minimal Model- Rates of change (SOS initial={})'.format(
            sos.num_molecules)
        plt.title(title)
        plt.savefig(os.path.join(dirname, title.replace(' ', '-') + '.png'))
        # plt.show()
        plt.close()

        # Rate vs RasGTP
        dirname = os.path.join(output_dirname, 'rate_balance')
        if not os.path.exists(dirname):
            os.mkdir(dirname)

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
        plt.title(title)
        plt.savefig(os.path.join(dirname, title.replace(' ', '-') + '.png'))
        # plt.show()
        plt.close()

    return final_vals


def find_steady_states(rates, sos_counts, ras, rasgap, rasgrp, guesses,
                       output_dirname='output'):
    '''
    Solve for steady states. We start with the guesses
    obtained from solving the ODE to plot it. We then bootstrap
    to get all the intermediate points.
    '''
    print(guesses)
    last = None
    steady_states = []

    for sos_mol in sos_counts:
        sos = Sos(num_molecules=sos_mol)
        reactions = Reactions(rates, sos.num_molecules,
                              ras.num_molecules, rasgap.num_molecules,
                              rasgrp.num_molecules)

        try:
            guess = guesses[sos_mol]
        except KeyError:
            guess = last
        steady_state_y = fsolve(dy_dt, guess, args=(0, reactions))
        print('Steady Y: {}'.format(steady_state_y))

        last = steady_state_y
        steady_states.append(steady_state_y)

    # Steady states by starting Sos molecules
    steady_states = np.array(steady_states)
    plt.figure(figsize=[10, 8])
    plt.plot(sos_counts, steady_states[:, 0], 'o', label='SOS')
    plt.plot(sos_counts, steady_states[:, 1], 'o', label='SOS-RasGTP')
    plt.plot(sos_counts, steady_states[:, 2], 'o', label='RasGTP')
    plt.legend()
    plt.xlabel('Initial SOS Count')
    plt.ylabel('Steady State Count')

    title = 'Das Minimal Model- Steady State RasGTP'
    dirname = os.path.join(output_dirname, 'steady_states')
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    plt.title(title)
    plt.savefig(os.path.join(dirname, title.replace(' ', '-') + '.png'))
    # plt.show()
    plt.close()

    return steady_states


def steady_states_by_rasgrp(rasgrp_steady_states, rasgrp_counts, sos_counts,
                            output_dirname='output'):
    '''
    Given the full set of steady state output from the SOS set,
    plot steady states when SOS is fixed while varying RasGRP1.
    '''
    # Right now we have rasgrp_count sets of sos_count steady states;
    # We want sos_count on the outside, so we move through columns
    # instead of rows
    rasgrp_steady_states = np.array(rasgrp_steady_states)
    for i, count in enumerate(sos_counts):
        steady_states = rasgrp_steady_states[:, i]

        # Steady states by starting Sos molecules
        plt.figure(figsize=[10, 8])
        plt.plot(rasgrp_counts, steady_states[:, 0], 'o', label='SOS')
        plt.plot(rasgrp_counts, steady_states[:, 1], 'o', label='SOS-RasGTP')
        plt.plot(rasgrp_counts, steady_states[:, 2], 'o', label='RasGTP')
        plt.legend()
        plt.xlabel('Initial RasGRP1 Count')
        plt.ylabel('Steady State Count')

        title = 'Das Minimal Model- Steady State RasGTP (SOS Initial={})'.format(
            count)
        dirname = os.path.join(output_dirname, 'steady_states')
        if not os.path.exists(dirname):
            os.mkdir(dirname)
        plt.title(title)
        plt.savefig(os.path.join(dirname, title.replace(' ', '-') + '.png'))
        # plt.show()
        plt.close()


if __name__ == '__main__':
    box = DiffusionBox()

    # Rates are listed in uM, from table S1
    # k_off and kcat rates are already in s^-1; k_on needs to be converted.
    rates = [box.convert_membrane_rate(0.12), 3.0,
             box.convert_membrane_rate(0.11), 0.4,
             box.convert_membrane_rate(0.05), 0.1, 0.038,
             box.convert_membrane_rate(0.07), 1.0, 0.003,
             box.convert_cytosolic_rate(1.74), 0.2, 0.1,
             0.01, 1.0, 0.01]
    rates = Rates(*rates)

    # We have used a intial Ras-GDP, and Ras-GAP concentrations of
    # 75 molecules/(um)2, 125 molecules/(um)3 respectively.
    ras = Ras(GDP=75)
    rasgap = RasGAP(125)

    rasgrp_final_vals, rasgrp_steady_states = [], []
    rasgrp_counts = range(0, 200, 10)
    sos_counts = range(0, 250, 10)
    for x in rasgrp_counts:

        # Vary RasGRP in the outer circle, SOS in inner.
        rasgrp = RasGRP(x)

        output_dirname = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            'output', 'rasgrp_{}'.format(x))
        if not os.path.exists(output_dirname):
            os.makedirs(output_dirname)

        final_vals = plot_molecule_count(rates, sos_counts, ras, rasgap,
                                         rasgrp, output_dirname)
        steady_states = find_steady_states(rates, sos_counts, ras, rasgap,
                                           rasgrp, final_vals, output_dirname)

        rasgrp_final_vals.append(final_vals)
        rasgrp_steady_states.append(steady_states)

    # Now we want to make the same steady state plot, but with respect
    # to varying amounts of RasGRP1, not SOS.
    # In rasgrp_steady_states, we have an outer loop of rasgrp values,
    # and we want the steady state when SOS is zero from each.
    output_dirname = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        'output', 'steady_state_by_rasgrp')
    if not os.path.exists(output_dirname):
        os.makedirs(output_dirname)
    steady_states_by_rasgrp(rasgrp_steady_states, rasgrp_counts, sos_counts,
                            output_dirname)
