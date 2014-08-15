'''
Created on Aug 13, 2014

@author: karmel
'''


class Reactions(object):

    '''
    Mean-field equations described in Experimental Procedures.
    Parameters are the concentrations of the specified species.
    '''
    rates = None
    Sos_tot = None
    Ras_tot = None
    RasGAP = None

    def __init__(self, rates, Sos_tot=None, Ras_tot=None, RasGAP=None):

        self.rates = rates
        self.Sos_tot = Sos_tot
        self.Ras_tot = Ras_tot
        self.RasGAP = RasGAP

        super().__init__()

    def dSos_dt(self, Sos, Sos_GTP, RasGTP,
                Sos_tot=None, Ras_tot=None):
        rates = self.rates

        if not Sos_tot:
            Sos_tot = self.Sos_tot
        if not Ras_tot:
            Ras_tot = self.Ras_tot
        print(Sos_tot, Sos, Sos_GTP)
        Sos_GDP = Sos_tot - Sos - Sos_GTP
        RasGDP = Ras_tot - RasGTP - Sos_GDP - Sos_GTP

        val = -rates.k1f * Sos * RasGDP \
            + rates.k1r * Sos_GDP \
            - rates.k2f * Sos * RasGTP \
            + rates.k2r * Sos_GTP
        return val

    def dRasGTP_dt(self, Sos, Sos_GTP, RasGTP,
                   Sos_tot=None, Ras_tot=None, RasGAP=None):
        rates = self.rates
        if not Sos_tot:
            Sos_tot = self.Sos_tot
        if not Ras_tot:
            Ras_tot = self.Ras_tot
        if not RasGAP:
            RasGAP = self.RasGAP

        Sos_GDP = Sos_tot - Sos_GTP
        RasGDP = Ras_tot - RasGTP - Sos_GDP - Sos_GTP
        val = -rates.k2f * Sos * RasGTP \
            + rates.k2r * Sos_GTP \
            + (rates.kcat3 * RasGDP * Sos_GTP) / (rates.K3m + RasGDP) \
            + (rates.kcat4 * RasGDP * Sos_GDP) / (rates.K4m + RasGDP) \
            - (rates.kcat5 * RasGAP * RasGTP) / (rates.K5m + RasGTP)
        return val

    def dSos_GDP_dt(self, Sos_GDP, RasGTP,
                    Sos_tot=None, Ras_tot=None):
        rates = self.rates
        if not Sos_tot:
            Sos_tot = self.Sos_tot
        if not Ras_tot:
            Ras_tot = self.Ras_tot
        Sos = Sos_tot - Sos_GDP
        RasGDP = Ras_tot - RasGTP
        val = rates.k1f * Sos * RasGDP - rates.k1r * Sos_GDP
        return val

    def dRasGTP_dt_simple(self, Sos_GDP, RasGTP,
                          Sos_tot=None, Ras_tot=None, RasGAP=None):
        '''
        Remove the SOS-RasGTP component for simplicity.
        '''
        rates = self.rates
        if not Sos_tot:
            Sos_tot = self.Sos_tot
        if not Ras_tot:
            Ras_tot = self.Ras_tot
        if not RasGAP:
            RasGAP = self.RasGAP

        Sos = Sos_tot - Sos_GDP
        RasGDP = Ras_tot - RasGTP
        print(Sos, RasGDP)
        val = (rates.kcat4 * RasGDP * Sos_GDP) / (rates.K4m + RasGDP) \
            - (rates.kcat5 * RasGAP * RasGTP) / (rates.K5m + RasGTP)
        return val
