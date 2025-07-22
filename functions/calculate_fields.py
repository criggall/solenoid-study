import numpy as np


def calculate_Br_analytic(data):

    ''' Calculate the radial magnetic field component via analytic approximation. '''

    Br_vals = []
    for i in range(len(data['z'])-1):
        x = data['x'].values[i]; y = data['y'].values[i]; z = data['z'].values[i]
        Bz = data['Bz'].values[i]
        r = np.sqrt(x**2+y**2)
        deltaBz = data['Bz'].values[i+1] - Bz
        deltaz = data['z'].values[i+1] - z
        if deltaz != 0:
            dBz_dz = deltaBz / deltaz
            Br = -r/2*dBz_dz
            Br_vals.append(Br)
        else:
            Br_vals.append(np.nan)

    Br_vals.append(np.nan)

    return Br_vals


def calculate_Br_rotated(data):

    ''' Calculate the radial magnetic field component via rotation. '''

    Br_vals = []; Bphi_vals = []
    for i in range(len(data['z'])-1):
        x = data['x'].values[i]; y = data['y'].values[i]; z = data['z'].values[i]
        Bx = data['Bx'].values[i]; By = data['By'].values[i]; Bz = data['Bz'].values[i]
        theta = np.arctan2(y,x)
        Br = Bx*np.cos(theta) + By*np.sin(theta)
        Bphi = -Bx*np.sin(theta) + By*np.cos(theta)
        Br_vals.append(Br)
        Bphi_vals.append(Bphi)
    Br_vals.append(np.nan)
    Bphi_vals.append(np.nan)

    return Br_vals, Bphi_vals