import numpy as np
from matplotlib import pyplot as plt
import os
import imageio.v2 as imageio

##### INPUTS #####

# Choose to run flipped or not flipped lattice:
polarity = 'flipped'
# polarity = 'not_flipped'

# Main directory:
main_dir = '/Users/criggall/Documents/solenoid-study/'

# Define working directory:
if polarity == 'flipped':
     main_dir = main_dir+'flipped/'
elif polarity == 'not_flipped':
     main_dir = main_dir+'not-flipped/'

# Parameter space scanned:
nominal_tilt = -0.0025*180/np.pi
tilt = np.linspace(0,nominal_tilt,10)

# Number of steps in scan:
iterations = len(tilt)

##### MATCHED REFERENCE PARTICLE DATA #####

# Read in matched reference particle data:
file_ref = main_dir+'AllTracks.txt'
data_ref = np.loadtxt(file_ref)

# Values along channel:
x_vals_ref = []; y_vals_ref = []; z_vals_ref = []
for i in range(data_ref.shape[0]):
    id = data_ref[i][8]
    if id == -2:
        x_vals_ref.append(data_ref[i][0]*0.1) # mm -> cm
        y_vals_ref.append(data_ref[i][1]*0.1)
        z = data_ref[i][2]*0.001 # mm -> m
        z_vals_ref.append(z)

##### PLOTTING FUNCTIONS #####

units = '$^{\circ}$'

# Define function to plot orbit in xy-plane:
def plot_orbit(x_vals, y_vals, z_vals, param_val, dir):
    plt.clf()
    plt.scatter(x_vals,y_vals,s=1)
    plt.title(f'solenoid tilt = {round(param_val,2)} {units}')
    plt.xlabel('x (cm)')
    plt.ylabel('y (cm)')
    if polarity == 'flipped':
        plt.xlim(-1.5,1.5)
        plt.ylim(-1.5,1.5)
    elif polarity == 'not_flipped':
        plt.xlim(-0.5,3.5)
        plt.ylim(-0.5,3.5)
    plt.savefig(dir+'orbit.png',dpi=300)
    plt.close()

# # Define function to plot dispersion along z:
# def plot_dispersion(D_vals, z_vals, param_val, dir):
#     plt.clf()
#     plt.figure(figsize=(10,4))
#     plt.scatter(z_vals, D_vals,s=1)
#     plt.title(f'solenoid tilt = {round(param_val,2)} {units}')
#     plt.xlabel('z (m)')
#     plt.ylabel('D (cm)')
#     plt.savefig(dir+'dispersion.png',dpi=300)
#     plt.close()

# Define function to plot Lz along z:
def plot_angular_momentum(Lz_vals, z_vals, param_val, dir):
    plt.clf()
    plt.figure(figsize=(10,4))
    plt.scatter(z_vals, Lz_vals,s=1)
    plt.title(f'solenoid tilt = {round(param_val,2)} {units}')
    if polarity == 'flipped':
        plt.ylim(-20,20)
    elif polarity == 'not_flipped':
        plt.ylim(-15,15)
    plt.xlabel('z (m)')
    plt.ylabel('$L_z$ (cm*MeV/c)')
    plt.savefig(dir+'angular_momentum.png',dpi=300)
    plt.close()

##### MAIN LOOP #####

full_channel_indices = []
x_total_residuals = []
y_total_residuals = []
for j in range(iterations):

    # Import data:
    dir = f'{main_dir}sol_tilt_scan/g4bl-output-sim{j+1}/'
    file = f'{dir}AllTracks.txt'
    data = np.loadtxt(file)

    # Values along channel:
    x_vals = []; y_vals = []; z_vals = []
    px_vals = []; py_vals = []
    Lz_vals = []
    for i in range(data.shape[0]):
        id = data[i][8]
        if id == -2:
            x = data[i][0]*0.1; y = data[i][1]*0.1 # mm -> cm
            x_vals.append(x); y_vals.append(y)
            z = data[i][2]*0.001 # mm --> m
            z_vals.append(z)
            px = data[i][3]; py = data[i][4]
            Lz = x*py - y*px
            Lz_vals.append(Lz)

    # # Compute dispersion:
    # D_vals = []
    # for i in range(len(z_vals)):
    #     D_val = (np.sqrt( (x_vals_ref[i] - x_vals[i])**2 + (y_vals_ref[i] - y_vals[i])**2 )) / (tilt[j]/200)
    #     D_vals.append(D_val)

    # Plot:
    plot_orbit(x_vals, y_vals, z_vals, tilt[j], dir)
    # plot_dispersion(D_vals, z_vals, tilt[j], dir)
    plot_angular_momentum(Lz_vals, z_vals, tilt[j], dir)

    # del x_vals, y_vals, z_vals, px_vals, py_vals, Lz_vals, D_vals
    del x_vals, y_vals, z_vals, px_vals, py_vals, Lz_vals

##### ANIMATIONS #####

frame_duration = 400

# List of sim directories:
out_dirs = [main_dir+f'sol_tilt_scan/g4bl-output-sim{i+1}' for i in range(iterations)]

# List of orbit plot files:
orbit_plot_paths = []
for i in range(iterations):
    orbit_plot_paths.append(out_dirs[i]+'/orbit.png')

# Create animation of orbit over scan:
orbit_plots = [imageio.imread(img) for img in orbit_plot_paths]
imageio.mimsave(main_dir+'orbit_scan_sol_tilt.gif', orbit_plots, duration=frame_duration, loop=0)

# # List of dispersion plot files:
# disp_plot_paths = []
# for i in range(iterations):
#     disp_plot_paths.append(out_dirs[i]+'/dispersion.png')

# # Create animation of dispersion over scan:
# disp_plots = [imageio.imread(img) for img in disp_plot_paths]
# imageio.mimsave(main_dir+'dispersion_scan_sol_tilt.gif', disp_plots, duration=frame_duration, loop=0)

# List of angular momentum plot files:
Lz_plot_paths = []
for i in range(iterations):
    Lz_plot_paths.append(out_dirs[i]+'/angular_momentum.png')

# Create animation of angular momentum over scan:
Lz_plots = [imageio.imread(img) for img in Lz_plot_paths]
imageio.mimsave(main_dir+'Lz_scan_sol_tilt.gif', Lz_plots, duration=frame_duration, loop=0)