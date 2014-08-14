'''
Created on Aug 13, 2014

@author: karmel

Parameters from the published supplement.
'''
from tcell.dasroose.minimal.constants import Rates
from tcell.dasroose.minimal.species import Ras, RasGAP, Sos
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.integrate import odeint
import numpy as np
from tcell.dasroose.minimal.equations import Reactions

if __name__ == '__main__':
    rates = [0.00018, 3.0, 
             0.00017, 0.4, 
             0.038, 1640,
             0.003, 9120, 
             0.1, 107,]
    
    rates = Rates(*rates)
    
    # We have used a intial Ras-GDP, and Ras-GAP concentrations of 
    # 75 molecules/(μm)2, 125 molecules/(μm)3 respectively.
    ras = Ras(GDP=75)
    rasgap = RasGAP(125)
    
    # Do the simplest thing first. Solve the ODE ignoring SOS-RasGTP
    sos = Sos(unbound=100)
    reactions = Reactions(rates, sos.concentration, 
                          ras.concentration, rasgap.concentration)
    print(sos.unbound, ras.GDP)
    # y = [Sos, RasGDP]. Solve dy/dt
    def dy_dt(y, t):
        Sos_GDP_i = y[0]
        RasGTP_i = y[1]
        print(y)
        f0 = reactions.dSos_GDP_dt(Sos_GDP_i, RasGTP_i)
        f1 = reactions.dRasGTP_dt_simple(Sos_GDP_i, RasGTP_i)
        print(f0,f1)
        return [f0, f1]
    
    # Initial conditions
    y0 = [sos.GDP, ras.GTP]
    t = np.linspace(0, 5., 1000)
    
    solution = odeint(dy_dt, y0, t)
    Sos_unbound = solution[:,0]
    RasGDP = solution[:,1]
    
    plt.figure()
    plt.plot(t, Sos_unbound, label='SOS-RasGDP')
    plt.plot(t, RasGDP, label='RasGTP')
    plt.legend()
    plt.show()
    
    