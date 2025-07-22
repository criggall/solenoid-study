import numpy as np


def calculateF(data, p=200e-3):

    ''' Calculates focusing length for single coil using formula provided by Katsuya. Uses natural units.
    
    Inputs:
    data = pandas DataFrame containing particle tracking data from G4bl -- only uses B_z values
    p = reference momentum (GeV) -- default is 200 MeV '''

    m = 105.7e-3 # GeV
    e = 0.303
    gamma = np.sqrt(1+(p/m)**2)

    peak_Bz = np.max(data['Bz'])

    f = 2*gamma*m/(e*peak_Bz) # m
    f = f*1000 # m --> mm

    return f