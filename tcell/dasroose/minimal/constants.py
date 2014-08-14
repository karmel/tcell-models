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
    
    k4f = 0
    k4r = 0
    kcat4 = 0
    
    k5f = 0
    k5r = 0
    kcat5 = 0
    
    def __init__(self, *args, **kwargs):
        '''
        Allow for keyword-specified rate setting OR just in-order lists.
        '''
        super().__init__()
        
        if kwargs:
            for key, val in kwargs.items():
                setattr(self, key, val)
        else:
            props = ('k1f', 'k1r', 'k2f', 'k2r',
                     'k3f', 'k3r', 'kcat3',
                     'k4f', 'k4r', 'kcat4',
                     'k5f', 'k5r', 'kcat5',)
            for i, prop in enumerate(props):
                try: 
                    setattr(self, key, args[i])
                except IndexError: break
            
        
    # Michaelis constants
    @property
    def K3m(self): return (self.kcat3 + self.k3r)/self.k3f
    
    @property
    def K4m(self): return (self.kcat4 + self.k4r)/self.k4f
    
    @property
    def K5m(self): return (self.kcat5 + self.k5r)/self.k5f