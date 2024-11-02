### Input File ###
import os
import math
import matplotlib.pyplot as plot

class Input:
    def __init__(self):
        
        self.parent_save_dir = os.getcwd() + '/'
        self.file_name       = 'NACA_0012'
        self.lib             = 'cupy'
        
        self.monitoring_options_list  =  []
        
        if self.lib == 'cupy':
            
            import cupy as npcp
            import numpy as np
            
        elif self.lib == 'numpy':
            
            import numpy as npcp
            
        else:
            print('Unknown library!') 

        self.restart_check  = False
        
        if self.restart_check == True:
            
            self.restart_im_count = 1             # Non-dimensional restart time step

        self.IBM_Check          = True
        self.diffusive_flux     = 'nonconservative'
   
        self.viscous_term       = 'constant_mu'
        
        self.num_species      = 1                   # number of species 
        self.num_dimension    = 2                   # number of dimension
        
        self.num_eqn = 1 + self.num_species + self.num_dimension # Number of equation (5 for 2D two species problem)
        
        ## PROBLEM GEOMETRY
         
        self.Nx = 4800                              # Resolution in x direction 
        self.Ny = 2880                               # Resolution in y direction
        
        self.x_a = -2.0                               # Physical location of starting point in x direction
        self.x_b =  6.0                              # Physical location of ending point in x direction
        self.y_a = -2.4                              # Physical location of starting point in y direction
        self.y_b =  2.4                          # Physical location of ending point in y direction
        
        self.dx  = (self.x_b - self.x_a) / self.Nx    # discretization width in x direction
        self.dy  = (self.y_b - self.y_a) / self.Ny    # discretization width in y direction
        
        self.N_cylinder_dx = 90
        self.num_layers = 1
        
        self.x   =  npcp.linspace(self.x_a + self.dx/2.0, self.x_b - self.dx/2.0, self.Nx)    # Computing the physical x location of cell centers array
        self.y   =  npcp.linspace(self.y_a + self.dy/2.0, self.y_b - self.dy/2.0, self.Ny)    # Computing the physical y location of cell centers array
        
        # Definition of nondimensional numbers for Cylinder
        self.aoa       = 8            # angle of attack in degrees 
        self.Re        = 1000.0         # Reynolds number for far field
        self.Ma        =  0.1          # Mach number
        self.TR        =  1.0          # Temperature ratio in our definition 
        self.Prandtl_1 =  0.72 
        self.gamma_1   =  7.0/5.0      # gas constant (for air 7/5)
        #self.Schmidt_1 =  1.0
        
        ## MATERIAL PROPERTIES
        self.species_m1   = 28.97/1000.0 
        self.R_universal  = 8.31446261815324  # J*mol^-1*K^-1
        
        # Initial values of problem
        self.L            = 1.0          # length of the airfoil
             
        self.T_initial    = 300.0        # initial temperature distribution for all domain (K)
        self.P_initial    = 101325.0     # initial far-field pressure distribution 1 atm: 101325 Pa
        self.rho_initial  = (self.P_initial /(self.T_initial * (self.R_universal/self.species_m1)))     # initial density distribution from ideal gas equation in kg/m^3    
        
        self.isothermal   = False
        
        self.u_initial    = self.Ma * math.sqrt(self.gamma_1 * (self.R_universal/self.species_m1) * self.T_initial)                           # initial far-field velocity distribution in x direction
        self.u_initial    = self.u_initial * npcp.cos(npcp.deg2rad(self.aoa))
        self.v_initial    = self.u_initial * npcp.sin(npcp.deg2rad(self.aoa))         # initial far-field velocity distribution in y direction
        self.c_v_1        = (self.R_universal / self.species_m1) / (self.gamma_1 - 1.0)  # constant volume heat capacity      
        self.c_p_1        = (self.R_universal / self.species_m1) + self.c_v_1            # constant pressure heat capacity
        
        self.c_mu_1       =  (self.rho_initial * self.u_initial * self.L) / self.Re   # dynamic viscosity is set to constant
        self.c_kappa_1    =  (self.c_p_1 * self.c_mu_1) / self.Prandtl_1 
        #self.c_D_1       = self.c_mu_1 / self.rho_initial / self.Schmidt_1 
                                                           
        self.Ly  =  self.y_b - self.y_a          # domain size in y direction 
     
        self.xc  =  (self.x[0] + self.x[-1]) / 8.0     # center points of the airfoil
        self.yc  =  (self.y[0] + self.y[-1]) / 2.0     # center points of the airfoil
        
        if self.lib == 'cupy':
                
            self.d_ip  = npcp.sqrt(2.0) * self.dx  + 10.0**(-12.0)                # image point distance from the immersed boundary
            self.V_IB  = npcp.zeros((4, self.Nx, self.Ny,), dtype=npcp.float64)   # dtype is given by numpy dtype argument

        elif self.lib == 'numpy':
                
            self.d_ip  = npcp.sqrt(2.0) * self.dx + 10.0**(-12.0)                 # image point distance from the immersed boundary
            self.V_IB  = npcp.zeros((4, self.Nx, self.Ny,), dtype=npcp.float64)   # dtype is given by numpy dtype argument
        
        # Immersed Body boundary conditions
        self.rho_ib  = self.rho_initial                                                 # immersed body density is given same as the initial fluid density
        self.u_ib    = 0.0                                                              # fixed cylinder
        self.v_ib    = 0.0                                                              # fixed cylinder 
        
        #Adiabatic case
        if self.isothermal == False:
            
            self.P_ib    =  self.P_initial
            
        #Isothermal case
        else:
            raise Exception('will not be implemented soon')
            self.P_ib    =  (self.Ra * self.P_initial * self.rho_ib) / self.rho_initial      # immersed body pressure depends on the Rayleigh number 
            self.T_ib    =  self.P_ib / ((self.gamma_1 - 1.0)* self.c_v_1 * self.rho_ib)     # immersed body initial temperature from ideal gas equation  

        self.V_IB[0,:,:] = self.rho_ib 
        self.V_IB[1,:,:] = self.u_ib
        self.V_IB[2,:,:] = self.v_ib
        self.V_IB[3,:,:] = self.P_ib
        
        # Boundary conditions  
        self.bdry_type_Left         = 'Dirichlet'
        self.bdry_type_Right        = 'Dirichlet'
        self.bdry_type_Bottom       = 'periodic'
        self.bdry_type_Top          = 'periodic'
        self.bdry_type_TopLeft      = 'periodic_y'
        self.bdry_type_TopRight     = 'periodic_y'
        self.bdry_type_BottomRight  = 'periodic_y'
        self.bdry_type_BottomLeft   = 'periodic_y'
        
        # Below variables are needed for Dircihlet and NeumannDirichlet Boundary Conditions. Values should be changed according to the problem.
        if self.check_dirichlet_left() or self.check_dirichlet_right() or self.check_dirichlet_top() or self.check_dirichlet_bottom():
                
            self.P_boundary    = self.P_initial    # in Pa
        
            # Use this one for constant temperature initial conditons and boundary 
            self.rho_boundary  = self.rho_initial
                
            if self.check_dirichlet_left():
                     
                self.u_boundary_left  = self.u_initial 
                self.v_boundary_left  = self.v_initial
                
            if self.check_dirichlet_right():
                    
                self.u_boundary_right = self.u_initial
                self.v_boundary_right = self.v_initial
                
            if self.check_dirichlet_bottom():
                    
                self.u_boundary_bottom = self.u_initial
                self.v_boundary_bottom = self.v_initial 
                    
            if self.check_dirichlet_top():
                    
                self.u_boundary_top = self.u_initial
                self.v_boundary_top = self.v_initial
        
        ## SOLVER INPUTS
        
        # Forcing switch
        self.Forcing_Check                    = True
        self.Gravity_X_Forcing_Check          = False
        self.Gravity_Y_Forcing_Check          = False
        self.Mean_Gravity_X_Forcing_Check     = False
        self.Mean_Gravity_Y_Forcing_Check     = False
        self.Channel_Forcing_Check            = False
        self.Sponge_Forcing_Check             = True
        self.Left_Sponge_Forcing_Check        = True
        self.Right_Sponge_Forcing_Check       = True
        self.Bottom_Sponge_Forcing_Check      = True
        self.Top_Sponge_Forcing_Check         = True
        
        if self.Gravity_X_Forcing_Check == True or self.Gravity_Y_Forcing_Check == True or self.Mean_Gravity_X_Forcing_Check == True or self.Mean_Gravity_Y_Forcing_Check == True:
            
            self.gravity   = [9.81, 9.81]
        
        ## Sponge Inputs
        if self.Sponge_Forcing_Check == True:

            self.rho_target  = self.rho_initial
            self.P_target    = self.P_initial

            if self.Left_Sponge_Forcing_Check == True:
                
                self.u_target_left    = self.u_initial
                self.v_target_left    = self.v_initial
                
                self.sponge_rate_left    = 30000.0                                          # sponge rate for left sponge region
                self.wl                  = 0.2                                             # left sponge region width 
                            
            if self.Right_Sponge_Forcing_Check == True:  
                
                self.u_target_right    = self.u_initial
                self.v_target_right    = self.v_initial
                
                self.sponge_rate_right   = 30000.00                                           # sponge rate for right sponge region
                self.wr                  =  0.2                                            # right sponge region width 
             
            if self.Top_Sponge_Forcing_Check == True:
                
                self.u_target_top    = self.u_initial
                self.v_target_top    = self.v_initial
                
                self.sponge_rate_top    = 30000.0
                self.wt                 =  0.2
            
            if self.Bottom_Sponge_Forcing_Check == True:
                
                self.u_target_bottom    = self.u_initial
                self.v_target_bottom    = self.v_initial
                
                self.sponge_rate_bottom   = 30000.0
                self.wb                   =  0.2      
        
        # Diffusive and Convective Order
        self.diffusive_scheme  = 'second-order CD'
        self.convective_scheme = 'second-order CD'
        
        if self.diffusive_scheme == 'tenth-order CD' or self.convective_scheme == 'tenth-order CD':
            
            self.boundary_order = 'tenth-order'
            self.num_ghosts     = 5
        
        elif self.diffusive_scheme == 'eighth-order CD' or self.convective_scheme == 'eighth-order CD':
                
            self.boundary_order = 'eighth-order'
            self.num_ghosts     = 4
        
        elif self.diffusive_scheme == 'sixth-order CD' or self.convective_scheme == 'sixth-order CD':
                
            self.boundary_order = 'sixth-order'
            self.num_ghosts     = 3
        
        elif self.diffusive_scheme == 'fourth-order CD' or self.convective_scheme == 'fourth-order CD':
                
            self.boundary_order = 'fourth-order'
            self.num_ghosts     = 2
        
        else:
            
            self.boundary_order = 'second-order'
            self.num_ghosts     = 1
        
        # Restart dump frequency
        self.t_print = 0.1 * (self.L / self.u_initial)
        
        # Simulation end time (dimensional)
        self.t_stop = 2000.0 * self.t_print 

        #Maximum number of steps
        self.max_step = 100000000

        #CFL
        self.CFL_num = 0.4

        # Runge Kutta Type
        self.RK_type = 2

        # Postprocess dump frequency
        self.t_postprocess =  self.t_stop
        
        # By default the region to dump in postprocess set to the original flow domain
        self.x_a_postprocess = self.x_a
        self.x_b_postprocess = self.x_b
        self.y_a_postprocess = self.y_a
        self.y_b_postprocess = self.y_b
        
        self.Nx_postprocess_start = int((self.x_a_postprocess - self.x_a) / self.dx)
        self.Nx_postprocess_end   = int((self.x_b_postprocess - self.x_a) / self.dx) + 2*self.num_ghosts
        self.Ny_postprocess_start = int((self.y_a_postprocess - self.y_a) / self.dy)
        self.Ny_postprocess_end   = int((self.y_b_postprocess - self.y_a) / self.dy) + 2*self.num_ghosts
        
    def check_dirichlet_left(self):
        
        if self.bdry_type_Left == 'Dirichlet' or self.bdry_type_Left == 'NeumannDirichlet':
            return True
        else:
            return False
        
    def check_dirichlet_right(self):
        
        if self.bdry_type_Right == 'Dirichlet' or self.bdry_type_Right == 'NeumannDirichlet':
            return True
        else:
            return False
        
    def check_dirichlet_bottom(self):
        
        if self.bdry_type_Bottom == 'Dirichlet' or self.bdry_type_Bottom == 'NeumannDirichlet':
            return True
        else:
            return False
        
    def check_dirichlet_top(self):
        
        if self.bdry_type_Top == 'Dirichlet' or self.bdry_type_Top == 'NeumannDirichlet':
            return True
        else:
            return False
    
    def check_dirichlet_front(self):
        
        if self.bdry_type_Front == 'Dirichlet' or self.bdry_type_Front == 'NeumannDirichlet':
            return True
        else:
            return False
    
    def check_dirichlet_back(self):
        
        if self.bdry_type_Back == 'Dirichlet' or self.bdry_type_Back == 'NeumannDirichlet':
            return True
        else:
            return False
    
    def check_adiabatic_left(self):
        
        if self.bdry_type_Left == 'Adiabatic':
            return True
        else:
            return False
        
    def check_adiabatic_right(self):
        
        if self.bdry_type_Right == 'Adiabatic':
            return True
        else:
            return False
    
    def check_adiabatic_bottom(self):
        
        if self.bdry_type_Bottom == 'Adiabatic':
            return True
        else:
            return False
        
    def check_adiabatic_top(self):
        
        if self.bdry_type_Top == 'Adiabatic':
            return True
        else:
            return False
    
    def check_adiabatic_front(self):
        
        if self.bdry_type_Front == 'Adiabatic':
            return True
        else:
            return False
    
    def check_adiabatic_back(self):
        
        if self.bdry_type_Back == 'Adiabatic':
            return True
        else:
            return False
        
    def check_bc_initiation_left(self):
        
        if self.check_dirichlet_left() or self.check_adiabatic_left():
            return True
        else:
            return False
    
    def check_bc_initiation_right(self):
        
        if self.check_dirichlet_right() or self.check_adiabatic_right():
            return True
        else:
            return False
    
    def check_bc_initiation_top(self):
        
        if self.check_dirichlet_top() or self.check_adiabatic_top():
            return True
        else:
            return False
    
    def check_bc_initiation_bottom(self):
        
        if self.check_dirichlet_bottom() or self.check_adiabatic_bottom():
            return True
        else:
            return False
    
    def check_bc_initiation_front(self):
        
        if self.check_dirichlet_front() or self.check_adiabatic_front():
            return True
        else:
            return False
    
    def check_bc_initiation_back(self):
        
        if self.check_dirichlet_back() or self.check_adiabatic_back():
            return True
        else:
            return False
