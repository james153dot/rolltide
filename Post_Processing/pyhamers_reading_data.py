# Importing modules
from pyhamers.readers import samrai_reader
import numpy as np
from matplotlib import ticker
import matplotlib.pyplot as plt
from matplotlib.pyplot import gca
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors
import matplotlib as mpl
import shutil
import os
import math 
#import latex

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
mpl.rcParams.update(mpl.rcParamsDefault)
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

parent_save_dir = os.getcwd() + '/'
#problem_name    = 'IBM'

'''
# Case parameters
Re = 150.0 
Ma = 0.4
Ra = 3.0
TR = Ra

thermal_boundary = 'isothermal'  # adiabatic or isothermal

viscosity_model = 'power_law' # constant or power_law
power = 0.75  # value for power law viscosity model

gamma  =  1.4
R_universal = 8.31446261815324 # universal gas constant (J*mol^-1*K^-1)
P_initial   = 101325.0         # initial far-field pressure distribution
d           = 1.0
max_step    = 200000


species_m   = (28.97 / 1000.0)   # molar mass of the gas kg/mol
R_specific  = R_universal / species_m
c_v         = (R_specific) / (gamma - 1.0)
T_initial   = 300.0            # initial temperature distribution for all domain (K)
rho_initial = (P_initial /(T_initial * (R_specific)))   # initial density distribution for all domain
v_initial   = 0.0              # initial far-field velocity distribution in y direction
u_initial   = Ma * math.sqrt(gamma * P_initial / rho_initial)
c_mu        = (rho_initial * u_initial * d) / Re
'''
xc  = 0.0
yc  = 0.0
zc  = 0.0


directory = "C/Users/james/Downloads/Archive"
        
# Initiating samrai reader
data_reader = \
    samrai_reader.SamraiDataReader(directory + '/viz_2D_uniform_flow', \
                                   periodic_dimensions = (False, False), \
                                   upsampling_method = 'sixth_order_Lagrange')
    
num_ghosts = (0, 0, 0)
steps = data_reader.steps
print(steps)


data_reader.step = steps[0]
domain_shape = data_reader.getRefinedDomainSize()

#reading the smaller domain by giving how much D you want to read from each side
x_sub_left   = 2.0
x_sub_right  = 40.0
y_sub_bottom = 4.5
y_sub_top    = 4.5
figsize_sub  = (int(x_sub_right+x_sub_left)/2, int(y_sub_bottom + y_sub_top)/2)

x_coords, y_coords = data_reader.getCombinedCoordinatesInSubdomainFromAllLevels(num_ghosts)

#refined regions mesh size
dx  = 1/120.0 #round((x_coords[1]-x_coords[0]), 3)
dy  = 1/120.0 #round((y_coords[1]-y_coords[0]), 3)
d_ip = (np.sqrt(2.0) * dx  + 10.0**(-12.0)) 
x_a = x_coords[0] - dx/2
y_a = y_coords[0] - dy/2
x_b = x_coords[-1] + dx/2
y_b = y_coords[-1] + dy/2

domain_shape = np.asarray(domain_shape)
lo_subdomain = domain_shape*0
hi_subdomain = domain_shape*0
lo_subdomain[0] = int(abs((xc-x_sub_left-x_a)/dx))
hi_subdomain[0] = int(abs((xc+x_sub_right-x_a)/dx))
lo_subdomain[1] = int(abs((yc-y_sub_bottom-y_a)/dy))
hi_subdomain[1] = int(abs((yc+y_sub_top-y_a)/dy))

Nx = hi_subdomain[0] - lo_subdomain[0] +1
Ny = hi_subdomain[1] - lo_subdomain[1] +1

data_reader.setSubDomain((lo_subdomain, hi_subdomain))

## Reading the IB Mask to decide cell types.

var_names = ('IB mask',)
            
data_reader.readCombinedDataInSubdomainFromAllLevels( \
            var_names, num_ghosts)
            
x_coords, y_coords = data_reader.getCombinedCoordinatesInSubdomainFromAllLevels(num_ghosts)
            
data_Y  = data_reader.getData('IB mask',)
            
data_masked = np.ma.masked_where(np.isnan(data_Y), data_Y) #cell type data is now here.

cell_type = np.zeros((Nx, Ny), dtype=int)

cell_type[:, :] = data_masked[:, :, 0] 
