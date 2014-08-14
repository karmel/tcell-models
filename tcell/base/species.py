'''
Created on Aug 13, 2014

@author: karmel
'''

class Molecule(object):
    concentration = None
    
    def __init__(self, concentration=None):
        super().__init()
        
        self.concentration = concentration