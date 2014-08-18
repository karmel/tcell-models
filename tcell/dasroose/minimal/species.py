'''
Created on Aug 13, 2014

@author: karmel

The molecular species that are involved in the minimal Das-Roose model of
TCR signaling via Ras and SOS.
'''
from tcell.base.species import Molecule
from tcell.dasroose.minimal.constants import DiffusionBox


class ConcentrationMolecule(Molecule):

    location = 'membrane'

    def __init__(self, concentration=None, **kwargs):

        if concentration is not None:
            box = DiffusionBox()
            if self.location == 'cytosolic':
                self.num_molecules = box.cytosolic_to_absolute(concentration)
            else:
                self.num_molecules = box.membrane_to_absolute(concentration)
        else:
            super().__init__(**kwargs)


class GTPDependent(Molecule):
    unbound = 0
    GDP = 0
    GTP = 0

    def __init__(self, GDP=0, GTP=0, unbound=0):
        box = DiffusionBox()

        # We want absolute number of molecules.
        self.GDP = box.membrane_to_absolute(GDP)
        self.GTP = box.membrane_to_absolute(GTP)
        self.unbound = box.cytosolic_to_absolute(unbound)

        super().__init__(self.GDP + self.GTP + self.unbound)


class Sos(GTPDependent):

    '''
    SOS (Son of sevenless) catalytic domain.
    Can be bound by GTP (very active) or GDP (less inactive).
    '''

    def __init__(self, num_molecules=None, **kwargs):
        '''
        SOS starts with number of molecules, not concentration.
        '''
        if num_molecules:
            self.unbound = num_molecules
            self.num_molecules = num_molecules
        else:
            super().__init__(**kwargs)


class Ras(GTPDependent):

    '''
    GTPase Ras. Can be bound by GTP (active) or GDP (inactive).
    '''


class RasGAP(ConcentrationMolecule):

    '''
    Converts Ras-GTP back to Ras-GDP--> inhibitor.
    '''

    location = 'cytosolic'


class RasGRP(ConcentrationMolecule):

    '''
    Converts Ras-GDP to Ras-GDP, albeit at a slower rate than SOS.
    '''

    location = 'cytosolic'
