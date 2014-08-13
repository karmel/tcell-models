'''
Created on Aug 13, 2014

@author: karmel

The molecular species that are involved in the minimal Das-Roose model of 
TCR signaling via Ras and SOS.
'''
from tcell.base.species import Molecule

GTP = 1
GDP = 0

class GTPDependent(Molecule):
    bound_by = None
    
    @property
    def is_active(self):
        '''
        When bound by RasGTP, 75x more activity --> on. Otherwise, off.
        '''
        return bool(self.bound_by)
    
class Sos(GTPDependent):
    '''
    SOS (Son of sevenless) catalytic domain.
    '''
    
class Ras(GTPDependent):
    '''
    GTPase Ras. Can be bound by GTP (active) or GDP (inactive). 
    '''
    
class GAPS(Molecule):
    '''
    Converts Ras-GTP back to Ras-GDP--> inhibitor.
    '''