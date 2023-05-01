#Planetary Data
#Purpose of this file is to provide planet data :)


Gm = 6.67430e-11
G = Gm * 10**-9

earth = {   'name'  : 'Earth',
            'mass'  : 5.972e24,
            'mu'    : 5.972e24 * G,
            'radius': 6378.0,
            'J2'    : 1.081874e-3,
            'SOI'   : 926006.6608, # km
            'map'   : 'Blues'
        }

custom = {   'name'  : 'Custom',
            'mass'  : 1,
            'mu'    : 1 * G,
            'radius': 0,
            'J2'    : 0,
            'SOI'   : 926006.6608, # km
            'map'   : 'Blues'
        }

