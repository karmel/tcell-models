'''
Created on Aug 13, 2014

@author: karmel
'''

class Reactions(object):
    '''
    Mean-field equations described in Experimental Procedures.
    Parameters are the concentrations of the specified species. 
    '''
    
    def dSos_dt(self, rates, Sos, Ras_GDP, Ras_GTP, Sos_GDP, Sos_GTP):
        val = -rates.k1f*Sos*Ras_GDP \
                + rates.k1r*Sos_GDP \
                - rates.k2f*Sos*Ras_GTP \
                + rates.k2r*Sos_GTP
        return val
    
    def dSos_GTP_dt(self, rates, Sos, Ras_GTP, Sos_GTP):
        val = rates.k2f*Sos*Ras_GTP - rates.k2r*Sos_GTP
        return val
        
    def dRas_GTP_dt(self, rates, Sos, Ras_GDP, Ras_GTP, Sos_GDP, Sos_GTP, Ras_Gap):
        val = -rates.k2f*Sos*Ras_GTP \
                + rates.k2r*Sos_GTP \
                + (rates.kcat3*Ras_GDP*Sos_GTP)/(rates.K3m + Ras_GDP) \
                + (rates.kcat4*Ras_GDP*Sos_GDP)/(rates.K4m + Ras_GDP) \
                - (rates.kcat5*Ras_Gap*Ras_GTP)/(rates.K5m + Ras_GTP)
        return val