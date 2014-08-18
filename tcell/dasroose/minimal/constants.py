'''
Created on Aug 13, 2014

@author: karmel
'''


class Rates(object):
    k1f = 0
    k1r = 0

    k2f = 0
    k2r = 0

    k3f = 0
    k3r = 0
    kcat3 = 0
    K3m = 0

    k4f = 0
    k4r = 0
    kcat4 = 0
    K4m = 0

    k5f = 0
    k5r = 0
    kcat5 = 0
    K5m = 0

    k8f = 0
    k8r = 0
    kcat8 = 0
    K8m = 0

    def __init__(self, *args, **kwargs):
        '''
        Allow for keyword-specified rate setting OR just in-order lists.
        '''
        super().__init__()

        if kwargs:
            for key, val in kwargs.items():
                setattr(self, key, val)
        else:
            props = ('k1f', 'k1r',
                     'k2f', 'k2r',
                     'k3f', 'k3r', 'kcat3',
                     'k4f', 'k4r', 'kcat4',
                     'k5f', 'k5r', 'kcat5',
                     'k8f', 'k8r', 'kcat8',)
            for i, prop in enumerate(props):
                try:
                    setattr(self, prop, args[i])
                except IndexError:
                    break

        # Calculate Michaelis Menten constants
        self.K3m = self.calc_Km(self.k3f, self.k3r, self.kcat3)
        self.K4m = self.calc_Km(self.k4f, self.k4r, self.kcat4)
        self.K5m = self.calc_Km(self.k5f, self.k5r, self.kcat5)
        self.K8m = self.calc_Km(self.k8f, self.k8r, self.kcat8)

    def calc_Km(self, kf, kr, kcat):
        '''
        Calculate Michaelis Menten constant from components.
        '''
        return (kcat + kr) / kf


class DiffusionBox(object):

    '''
    The molecules are assumed to be in a Diffusion box, evenly distributed.
    '''
    d = 0.0017  # microns
    A = 4.0  # microns^2
    V = 0.08  # microns^3

    mols_per_uM = 600

    def absolute_to_membrane(self, num_molecules):
        '''
        Given an absolute number of molecules, get the membrane concentration.
        '''
        return num_molecules / self.A

    def absolute_to_cytosolic(self, num_molecules):
        '''
        Given an absolute number of molecules, get the cytosolic concentration.
        '''
        return num_molecules / self.V

    def membrane_to_absolute(self, concentration):
        '''
        Given the membrane concentration, return the absolute
        number of molecules.
        '''
        return concentration * self.A

    def cytosolic_to_absolute(self, concentration):
        '''
        Given the cytosolic concentration, return the absolute
        number of molecules.
        '''
        return concentration * self.V

    def convert_cytosolic_rate(self, rate):
        '''
        We want to convert from uM/s rates to just s^-1 rates,
        assuming our diffusion box.

        According to the supp, 1uM = 600 mol/micron^3, and this is happening
        in the small diffusion box specified.

        So, for cytosolic rates, we want the uM rate/600/volume
        to get our adjusted rate.
        '''

        return rate / self.mols_per_uM / self.V

    def convert_membrane_rate(self, rate):
        '''
        We want to convert from uM/s rates to just s^-1 rates,
        assuming our diffusion box.

        According to the supp, 1uM = 600 mol/micron^3, and this is happening
        in the small membrane-proximal box specified.

        So, for membrane rates, we want the uM rate/600/A*d
        to get our adjusted rate.
        '''

        return rate / self.mols_per_uM / (self.A * self.d)
