'''
Created on Aug 13, 2014

@author: karmel

The molecular species that are involved in the minimal Das-Roose model of 
TCR signaling via Ras and SOS.
'''
from tcell.base.species import Molecule

class GTPDependent(Molecule):
    unbound = 0
    GDP = 0
    GTP = 0
    
    def __init__(self, GDP=0, GTP=0, unbound=0):
        self.GDP = GDP
        self.GTP = GTP
        self.unbound = unbound
        
        super().__init__(GDP + GTP + unbound)
        
    def add_GDP(self, amt=1):
        self.add_mod('GDP', amt)
    def remove_GDP(self, amt=1):
        self.remove_mod('GDP', amt)
    def add_GTP(self, amt=1):
        self.add_mod('GTP', amt)
    def remove_GTP(self, amt=1):
        self.remove_mod('GTP', amt)
    
    def update_concentration(self):
        self.concentration = self.GDP + self.GTP + self.unbound
        
    def add_mod(self, mod, amt):
        '''
        Total amount stays the same, so make sure to keep everything in check.
        '''
        setattr(self, mod, getattr(self, mod) + amt)
        self.unbound -= amt
        self.update_concentration()
        
    def remove_mod(self, mod, amt):
        '''
        Total amount stays the same, so make sure to keep everything in check.
        '''
        setattr(self, mod, getattr(self, mod) - amt)
        self.unbound += amt
        self.update_concentration()
        
class Sos(GTPDependent):
    '''
    SOS (Son of sevenless) catalytic domain.
    Can be bound by GTP (very active) or GDP (less inactive).
    '''    
        
class Ras(GTPDependent):
    '''
    GTPase Ras. Can be bound by GTP (active) or GDP (inactive). 
    '''
    
class RasGAP(Molecule):
    '''
    Converts Ras-GTP back to Ras-GDP--> inhibitor.
    '''