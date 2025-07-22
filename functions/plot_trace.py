import numpy as np
from matplotlib import pyplot as plt


def plot_solenoids(n, d, L, zshift):

    ''' Plots locations and extent of solenoids in the lattice.
    
    Inputs: 
    n = number of solenoids
    d = distance between solenoids 
    L = length of solenoid '''

    for i in range(n):
        center = i*d + zshift
        plt.axvspan(xmin=center-L/2, xmax=center+L/2, color='gray', alpha=0.2)


def plot_Bz(data):

    ''' Plots longitudinal B field component. '''

    plt.scatter(data['z'], data['Bz'], color='green', s=1)
    plt.xlabel('$z$ [mm]')
    plt.ylabel('$B_z$ [T]')


def plot_Btrans(data):

    ''' Plots transverse B field components. '''

    plt.scatter(data['z'], data['Bx'], color='red', label='Bx', s=1)
    plt.scatter(data['z'], data['By'], color='blue', label='By', s=1)
    plt.legend()
    plt.xlabel('$z$ [mm]')
    plt.ylabel('$B_z$ [T]')


def plot_Lz(data):

    ''' Plots angular momentum. '''

    plt.scatter(data['z'], data['Lz'], color='black', s=1)
    plt.xlabel('$z$ [mm]')
    plt.ylabel(r'$L_z$ [mm$\cdot$MeV/c]')


def plot_xy(data):

    ''' Plot trajectory in xy-plane. '''

    plt.scatter(data['x'], data['y'], color='black', s=1)
    plt.xlabel('$x$ [mm]')
    plt.ylabel('$y$ [mm]')