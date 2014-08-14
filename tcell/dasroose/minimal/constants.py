'''
Created on Aug 13, 2014

@author: karmel
'''

class Rates(object):
    k1f = 0
    k1r = 0 
    
    k2f = 0
    k2r = 0
    
    kcat3 = 0
    K3m = 0
    
    kcat4 = 0
    K4m = 0
    
    kcat5 = 0
    K5m = 0
    
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
                     'kcat3', 'K3m',
                     'kcat4', 'K4m',
                     'kcat5', 'K5m',)
            for i, prop in enumerate(props):
                try: 
                    setattr(self, prop, args[i])
                except IndexError: break
            
        
    