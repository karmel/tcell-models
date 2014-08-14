'''
Created on Aug 13, 2014

@author: karmel

Parameters from the published supplement.
'''
from tcell.dasroose.minimal.constants import Rates
from tcell.dasroose.minimal.species import Ras, RasGAP
import matplotlib.pyplot as plt
import seaborn as sns

if __name__ == '__main__':
    rates = [0.12, 3.0, 0.11, 0.4, 
             0.05, 0.1, 0.038,
             0.07, 1.0, 0.003,
             1.74, 0.2, 0.1]
    
    rates = Rates(*rates)
    
    # We have used a intial Ras-GDP, and Ras-GAP concentrations of 
    # 75 molecules/(μm)2, 125 molecules/(μm)3 respectively.
    ras = Ras(GDP=75)
    rasgap = RasGAP(125)
    
    