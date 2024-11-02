## IBM Calculations for Airfoil NACA 0012 ##
import Input
import matplotlib.pyplot as plot

#
# Find a way to attach a curcle to the edge, making it fully continous over the end.
# Right now the edge is taken away by a circle, but we want to make it so that the circle isn't discontinious.
# Currently there is a jump: we want to have tangent lines / dervivatives be the same along the entire airfoil.
# The new circle has to fit fully underneath the original airfoil.
#

class IBM_Calculations:
    def __init__(self,):
            
        self.input = Input.Input()
        
        self.lib       = self.input.lib
        
        if self.lib == 'cupy':
            print('IBM Calculations include many for loops. CuPy library may cause significant slow down! NumPy is used!!')
            
            import numpy as npcp
            import cupy as cp
            
        elif self.lib == 'numpy':
        
            import numpy as npcp
        
        else:
            print('Unknown library!')
        
        self.Nx        = self.input.Nx
        self.Ny        = self.input.Ny
        
        if self.lib == 'cupy':
                
            self.x   = cp.asnumpy(self.input.x)
            self.y   = cp.asnumpy(self.input.y)
            self.xc  = npcp.float64(cp.asnumpy(self.input.xc))  
            self.yc  = npcp.float64(cp.asnumpy(self.input.yc)) 
               
        elif self.lib == 'numpy': 
                
            self.x  = self.input.x
            self.y  = self.input.y
            self.xc  = self.input.xc
            self.yc  = self.input.yc
        
        ## PROBLEM GEOMETRY

        self.x_a       = self.input.x_a
        self.x_b       = self.input.x_b
        self.y_a       = self.input.y_a
        self.y_b       = self.input.y_b
        
        # Geometric Operations
        self.dx   = self.input.dx
        self.dy   = self.input.dy
        self.d_ip = self.input.d_ip
        
        self.N_cylinder_dx = self.input.N_cylinder_dx # cylinder 0012 - amount cutted from sharp edge to save distance.
        
        self.num_layers = self.input.num_layers
        
            # Finding number of ibm ghost cells  and deciding and storing Cell Types (Level Set Function)
            
        self.distance_to_center      =  npcp.zeros((self.Nx, self.Ny, ), dtype=npcp.float64)
        #2d array with distances Nx and Ny
        self.d_gc                    =  npcp.zeros((self.Nx, self.Ny, ), dtype= npcp.float64)
        #2d array, identical
        self.L_dx                    =  int(self.input.L / self.dx)
        self.half_L_dx               =  int(self.L_dx / 2)  
        self.xc_left                 =  int((self.xc - self.x_a - self.dx/2.0) // self.dx)
        #some left boundary?
        self.cell_type               =  npcp.ones((self.Nx, self.Ny, ), dtype= npcp.float64) # Assigning first all cells as fluids
        #all ceels are fluids right now
        print(self.L_dx, self.half_L_dx, self.xc_left, self.xc, self.x[self.xc_left - self.half_L_dx + 1])
        
        #Finding Fluid and Solid Cells
        self.num_object_cell = 0
        # of cells
        for i in range(self.L_dx):
            #loops over, need a k
            for j in range(self.Ny): #Loops over all of Ny
                
                if self.y[j] > self.yc:   # Upper half of the airfoil. 
                    
                    if (self.y[j]) < self.NACA_0012_y_chopped(self.x[self.xc_left - self.half_L_dx + i + 1 ]) :
                        #likely below y-value boundary for given x-coordiante.
                        
                        self.cell_type[self.xc_left - self.half_L_dx + i + 1, j] = -1.0 # These cells are inside the airfoil. Assign them as solid cells.
                        self.num_object_cell +=1
                        # if solid, then the cell type is = -1, then object cell +1
                else: #Lower half of the airfoil
                    
                    if (self.y[j]) > -self.NACA_0012_y_chopped(self.x[self.xc_left  - self.half_L_dx + i + 1] ) :
                        
                        self.cell_type[self.xc_left - self.half_L_dx + i  + 1, j] = -1.0 # These cells are inside the airfoil. Assign them as solid cells.
                        self.num_object_cell +=1
                        
        # what does this look like?
        plot.figure()
        plot.imshow(npcp.rot90(self.cell_type), extent=[self.x_a, self.x_b, self.y_a, self.y_b])
        
        for i in range(self.Nx):
           
            plot.plot((self.x[i] - self.dx/2, self.x[i] - self.dx/2), (self.y_a, self.y_b), 'k',linewidth = '0.2')
        
        for j in range(self.Ny):
            
            plot.plot((self.x[0] - self.dx/2, self.x[-1] + self.dx/2), (self.y[j] - self.dy/2, self.y[j] - self.dy/2),  'k',linewidth = '0.2')
        
        ax = plot.gca()
        ax.get_xaxis().set_visible(False) 
        ax.get_yaxis().set_visible(False)       
        plot.savefig('airfoil-cell-types-fluid-solid.jpg',bbox_inches = 'tight', dpi=1200)
        plot.close()      
        
        # Finding Ghost Cells within Solid Cells
    
        self.num_ibm_ghost_cell_first_layer= 0
        
        for i in range(self.Nx):
            
            for j in range(self.Ny):
                
                if self.cell_type[i, j] == -1.0:  # If it is inside the body 
                                     
                    if (self.cell_type[i +1, j] == 1.0 or self.cell_type[i -1, j] == 1.0 
                        or self.cell_type[i, j+1] == 1.0 or self.cell_type[i, j-1] == 1.0): 
                        # checks if neighoring cells are fluid, if so marks the current cell as a ghost cell.
                        # or self.cell_type[i +1, j+1] == 1.0 or self.cell_type[i +1, j-1] == 1.0 
                        # or self.cell_type[i -1, j+1] == 1.0 or self.cell_type[i -1, j-1] == 1.0): 
                                 
                            self.cell_type[i, j] = 0.0  # is ghost cell
                            # If one of the cells at its left, right, top, bottom, bottom-left, bottom-right, 
                            #top-left and top-right location is fluid then it is a ghost cell.        
                            
                            self.num_ibm_ghost_cell_first_layer +=1 # is ghost +1
                            self.num_object_cell                -=1 # not solid
        
        # plot.figure()
        # plot.imshow(npcp.rot90(self.cell_type), extent=[self.x_a, self.x_b, self.y_a, self.y_b])
        # plot.imshow(npcp.rot90(self.cell_type[351:551, 351:551]), extent=[-0.5, 1.5, -0.5, 1.5])
        
        for i in range(self.Nx):
           
            plot.plot((self.x[i] - self.dx/2, self.x[i] - self.dx/2), (self.y_a, self.y_b), 'k',linewidth = '0.2')
        
        for j in range(self.Ny):
            
            plot.plot((self.x[0] - self.dx/2, self.x[-1] + self.dx/2), (self.y[j] - self.dy/2, self.y[j] - self.dy/2),  'k',linewidth = '0.2')
        
        ax = plot.gca()
        ax.get_xaxis().set_visible(False) 
        ax.get_yaxis().set_visible(False)       
        plot.savefig('airfoil-cell-types-fluid-solid-ghost.jpg',bbox_inches = 'tight', dpi=1200)
        plot.close() 
        
        #raise Exception('Below will be considered')
    
        self.ibm_total_ghost_cell = self.num_ibm_ghost_cell_first_layer 
        # number of IBM ghost cells = number of first-layer ghost cells
        
        # Creating Boundary Points - points on the boundary
        num_bp    = 10000
        
        self.x_bp =  npcp.zeros((num_bp, ), dtype=npcp.float64)
        self.y_bp =  npcp.zeros((num_bp, ), dtype=npcp.float64)
        
        #x-center?
        xc = 1.0 - self.N_cylinder_dx * self.dx  
        R  =  self.NACA_0012_y(xc)
        
        x_bp_increment = (2.0 * (xc + R)) / num_bp
        #total length of boundary / number of boundary points
        for i in range(int(num_bp/2)+1):

            self.x_bp[i] = 0.0 + x_bp_increment * i 
            self.y_bp[i] = self.NACA_0012_y_chopped(self.x_bp[i])
        # generates first half of boundary points
        # finds the x point then the y point
        for i in range(int(num_bp/2)+1, num_bp):
            
            self.x_bp[i] = self.x_bp[num_bp  -i]
            self.y_bp[i] = -self.y_bp[num_bp -i]      
        # does the same thing but reflects to save computational power 

        print(self.num_ibm_ghost_cell_first_layer)
        
        # Finding the Boundary Intercept points for each ghost cells using boundary points created above
        self.x_bi = npcp.zeros((self.Nx, self.Ny, ), dtype=npcp.float64)
        self.y_bi = npcp.zeros((self.Nx, self.Ny, ), dtype=npcp.float64)
        
        for i in range(self.Nx):
            for j in range(self.Ny):
                # over every cell in the grid
                if self.cell_type[i,j] == 0.0:
                    # if ghost cell,
                    x_gc = self.x[i]
                    y_gc = self.y[j]
                     # store ghost cell
                    distance_gc_bp = npcp.zeros((num_bp, ), dtype=npcp.float64)
                    # find distance between ghost cell and boundary point
                    for k in range(num_bp):
                        
                        distance_gc_bp[k] = npcp.sqrt((x_gc - self.x_bp[k])**2.0 + (y_gc - self.y_bp[k])**2.0)
                        # finds the distance between ghost cell/boundary.
                    index_min = npcp.where(distance_gc_bp == npcp.min(distance_gc_bp))    
                    # some minimum distance
                    #print(index_min)
                    self.x_bi[i,j] = self.x_bp[index_min] 
                    self.y_bi[i,j] = self.y_bp[index_min]
                    #closest boundary point.
    
        '''
        # Finding the Boundary Intercept points for each ghost cells
        self.x_bi    = npcp.zeros((self.Nx, self.Ny, ), dtype=npcp.float64)
        self.y_bi    = npcp.zeros((self.Nx, self.Ny, ), dtype=npcp.float64)
               
        for i in range(self.Nx):
            
            for j in range(self.Ny):
                
                if self.cell_type[i,j] == 0.0:
                   
                    t    = 0.12
                    x_gc = self.x[i]
                    y_gc = self.y[j]
                    f = lambda x: npcp.sqrt( (x - x_gc) ** 2.0 + ((5.0 * t * ((0.2969 * x **0.5) - (0.1260 * x) - (0.3516 * x**2.0) + (0.2843 * x**3.0) - (0.1015 * x**4.0))) - y_gc) **2.0)
                    root = scipy.optimize.fminbound(f, 0.0, 1.0, xtol=1e-09)
                    
                    self.x_bi[i,j] = root
                    
                    if y_gc > self.yc:
                        
                        self.y_bi[i,j] = self.NACA_0012_y_chopped(root)
                        
                    else: 
                        
                        self.y_bi[i,j] = -self.NACA_0012_y_chopped(root)
        '''
        
        
        # Finding Image Points and bilinear interpolation coefficients 
        self.d_gc         = npcp.zeros((self.Nx, self.Ny,), dtype=npcp.float64)
        self.d_ip_dynamic = npcp.zeros((self.Nx, self.Ny,), dtype=npcp.float64)
        self.x_ip        = npcp.zeros((self.Nx, self.Ny,), dtype=npcp.float64)
        self.y_ip        = npcp.zeros((self.Nx, self.Ny,), dtype=npcp.float64)
        self.x_low       = npcp.zeros((self.Nx, self.Ny,), dtype=npcp.float64)
        self.y_low       = npcp.zeros((self.Nx, self.Ny,), dtype=npcp.float64)
        self.coef_x_ip   = npcp.zeros((self.Nx, self.Ny, ), dtype=npcp.float64)
        self.coef_y_ip   = npcp.zeros((self.Nx, self.Ny, ), dtype=npcp.float64)
        self.norm_0      = npcp.zeros((self.Nx, self.Ny, ), dtype=npcp.float64)
        self.norm_1      = npcp.zeros((self.Nx, self.Ny, ), dtype=npcp.float64)
        self.x_low_bi    = npcp.zeros((self.Nx, self.Ny, ), dtype=npcp.float64)
        self.y_low_bi    = npcp.zeros((self.Nx, self.Ny, ), dtype=npcp.float64)
        self.coef_x_bi   = npcp.zeros((self.Nx, self.Ny, ), dtype=npcp.float64)
        self.coef_y_bi   = npcp.zeros((self.Nx, self.Ny, ), dtype=npcp.float64)   
        
        # Finding and storing the image points
        for i in range(self.Nx):
            
            for j in range(self.Ny):
                
                if self.cell_type[i,j] == 0.0:
                    
                        
                    self.d_gc[i,j] = npcp.sqrt((self.x_bi[i,j] - self.x[i])**2.0 + (self.y_bi[i,j] - self.y[j])**2.0)  #finds the distance between two points          
                    # iternates over each cell, if ghost finds distance from ghost and boundary cell
                    
                    #self.alfa[i,j] = cp.arctan((self.y_bi[i,j] - self.y[j]) / (self.x_bi[i,j] - self.x[j]))
                    
                    self.norm_0[i,j] = (self.x_bi[i, j] - self.x[i]) / self.d_gc[i,j] #norm
                    self.norm_1[i,j] = (self.y_bi[i, j] - self.y[j]) / self.d_gc[i,j]
                    # find normal unit vectors for ghost cell along the x and y
                    condition = True
                    
                    while condition == True: # this finds the intersection points based on the distance d_ip and the direction.
                        
                        if self.x_bi[i, j] > self.x[i]: # if x_bi > x, then find new value
                        
                            self.x_ip[i,j] = self.x[i] + npcp.sqrt((self.d_gc[i,j] + self.d_ip)**2.0 - ((self.d_gc[i,j] + self.d_ip)*self.norm_1[i,j])**2.0)
                        
                        else:
                            
                            self.x_ip[i,j] = self.x[i] - npcp.sqrt((self.d_gc[i,j] + self.d_ip)**2.0 - ((self.d_gc[i,j] + self.d_ip)*self.norm_1[i,j])**2.0)
                        
                        if self.y_bi[i, j] > self.y[j]:
                            
                            self.y_ip[i,j] = self.y[j] + npcp.sqrt((self.d_gc[i,j] + self.d_ip)**2.0 - ((self.d_gc[i,j] + self.d_ip)*self.norm_0[i,j])**2.0)
                        
                        else:
                            
                            self.y_ip[i,j] = self.y[j] - npcp.sqrt((self.d_gc[i,j] + self.d_ip)**2.0 - ((self.d_gc[i,j] + self.d_ip)*self.norm_0[i,j])**2.0)
                        
                        '''
                        if self.x_bi[i,j] > self.x[i]:
                            
                            self.x_ip[i,j] = self.x_bi[i,j] + self.d_ip * cp.cos(self.alfa[i,j])
                            self.y_ip[i,j] = self.y_bi[i,j] + self.d_ip * cp.sin(self.alfa[i,j])
                        
                        else:
                            
                            self.x_ip[i,j] = self.x_bi[i,j] - self.d_ip * cp.cos(self.alfa[i,j])
                            self.y_ip[i,j] = self.y_bi[i,j] - self.d_ip * cp.sin(self.alfa[i,j])
                        '''
                        # calculate the indices and coefficients for the grid cells containing the intersection points and boundary intercept points.
                        self.x_low[i,j]      = (self.x_ip[i,j] - self.x_a - self.dx/2.0) // self.dx
                        self.x_low_bi[i,j]   = (self.x_bi[i,j] - self.x_a - self.dx/2.0) // self.dx
                            
                        self.coef_x_ip[i,j]  = (self.x_ip[i,j] - self.x[int(self.x_low[i,j])]) / (self.dx)    
                        self.coef_x_bi[i,j]  = (self.x_bi[i,j] - self.x[int(self.x_low_bi[i,j])]) / (self.dx)
                            
                        self.y_low[i,j]      = (self.y_ip[i,j] - self.y_a - self.dy/2.0) // self.dy
                        self.y_low_bi[i,j]   = (self.y_bi[i,j] - self.y_a - self.dy/2.0) // self.dy
                            
                        self.coef_y_ip[i,j]  = (self.y_ip[i,j] - self.y[int(self.y_low[i,j])]) / (self.dy)
                        self.coef_y_bi[i,j]  = (self.y_bi[i,j] - self.y[int(self.y_low_bi[i,j])]) / (self.dy)
                        
                        # check if the cells around intesction are fluid. If true, breaks the loop and changes d_ip_dynamic
                        if self.cell_type[int(self.x_low[i, j]), int(self.y_low[i, j])] == 1.0 and self.cell_type[int(self.x_low[i, j])+1, int(self.y_low[i, j])] == 1.0 and \
                            self.cell_type[int(self.x_low[i, j]), int(self.y_low[i, j])+1] == 1.0 and self.cell_type[int(self.x_low[i, j])+1, int(self.y_low[i, j])+1] == 1.0:
                             
                            condition = False
                            self.d_ip_dynamic[i, j] = self.d_ip
                            self.d_ip = npcp.sqrt(2.0) * self.dx + 1e-12
                        
                        else: # if false, changes the distance d_ip and prints
                            self.d_ip += npcp.sqrt(2.0) * self.dx
                            print('WARNING! cell type error')
                                                
        self.GC_int       = npcp.zeros((self.ibm_total_ghost_cell, 2), dtype = npcp.int64)          # CuPy may not support object data type arrays. If it's the case, arrays can be separated. 
        self.IP_int       = npcp.zeros((self.ibm_total_ghost_cell, 2), dtype = npcp.int64)
        self.BI_int       = npcp.zeros((self.ibm_total_ghost_cell, 2), dtype = npcp.int64)
        
        self.GC_float     = npcp.zeros((self.ibm_total_ghost_cell, 3), dtype = npcp.float64)          # CuPy may not support object data type arrays. If it's the case, arrays can be separated. 
        self.IP_float     = npcp.zeros((self.ibm_total_ghost_cell, 3), dtype = npcp.float64)
        self.BI_float     = npcp.zeros((self.ibm_total_ghost_cell, 4), dtype = npcp.float64)
        
        #stores the data for GC, IP, BI.
        
        #Storing the necessary parameters of ibm ghost cells and image points
        counter = 0
        # purely storing purposes
        for i in range(self.Nx):
            
            for j in range(self.Ny):
                
                if self.cell_type[i,j] == 0.0:
                        
                        # GC array consists of x, y, and the distance to the boundary of corresponding ghost cells
                        self.GC_int[counter, 0]     =  i          
                        self.GC_int[counter, 1]     =  j
                        
                        self.GC_float[counter, 0]   =  self.d_gc[i, j]
                        self.GC_float[counter, 1]   =  self.norm_0[i,j]
                        self.GC_float[counter, 2]   =  self.norm_1[i,j]
                        
                        # IP array consists of x_low, y_low, coefficient x, and coefficient y of corresponding image points
                        self.IP_int[counter, 0]     = int(self.x_low[i, j])
                        self.IP_int[counter, 1]     = int(self.y_low[i, j])
                        
                        self.IP_float[counter, 0]   = self.coef_x_ip[i, j]
                        self.IP_float[counter, 1]   = self.coef_y_ip[i, j]                        
                        self.IP_float[counter, 2]   = self.d_ip_dynamic[i, j]
                        
                        #BI array consists of x_low, y_low, coefficient x, coefficient y, x location, and y location of corresponding body intercept points
                        self.BI_int[counter, 0]     = int(self.x_low_bi[i, j])
                        self.BI_int[counter, 1]     = int(self.y_low_bi[i, j])
                        
                        self.BI_float[counter, 0]   = self.coef_x_bi[i, j]
                        self.BI_float[counter, 1]   = self.coef_y_bi[i, j]
                        self.BI_float[counter, 2]   = self.x_bi[i,j]
                        self.BI_float[counter, 3]   = self.y_bi[i,j]
                        
                        counter          +=1

        # Storing the necessary parameters of object cells
        
        self.OC = npcp.zeros((self.num_object_cell, 2), dtype = int)
        
        counter = 0
        
        for i in range(self.Nx):
                for j in range(self.Ny):
                        
                        if self.cell_type[i,j] == -1.0:
                                
                                self.OC[counter, 0] = i
                                self.OC[counter, 1] = j
                                counter        += 1
        
        # if there is only one layer of ghost cells, starts arrays for 1st layer and stops the rest.
        
        if self.num_layers == 1:
           
            self.GC_1_int       = npcp.zeros((self.num_ibm_ghost_cell_first_layer, 2), dtype = npcp.int64)          # CuPy may not support object data type arrays. If it's the case, arrays can be separated. 
            self.IP_1_int       = npcp.zeros((self.num_ibm_ghost_cell_first_layer, 2), dtype = npcp.int64)
            self.BI_1_int       = npcp.zeros((self.num_ibm_ghost_cell_first_layer, 2), dtype = npcp.int64)
            
            self.GC_1_float     = npcp.zeros((self.num_ibm_ghost_cell_first_layer, 3), dtype = npcp.float64)          # CuPy may not support object data type arrays. If it's the case, arrays can be separated. 
            self.IP_1_float     = npcp.zeros((self.num_ibm_ghost_cell_first_layer, 3), dtype = npcp.float64)
            self.BI_1_float     = npcp.zeros((self.num_ibm_ghost_cell_first_layer, 4), dtype = npcp.float64)
            
            counter = 0
            
            for i in range(self.Nx):
                for j in range(self.Ny):
                    
                    if self.cell_type[i,j] == 0.0:
                            
                        # GC array consists of x, y, and the distance to the boundary of corresponding ghost cells
                        self.GC_1_int[counter, 0]     =  i          
                        self.GC_1_int[counter, 1]     =  j
                        
                        self.GC_1_float[counter, 0]   =  self.d_gc[i, j]
                        self.GC_1_float[counter, 1]   =  self.norm_0[i,j]
                        self.GC_1_float[counter, 2]   =  self.norm_1[i,j]
                        
                        # IP array consists of x_low, y_low, coefficient x, and coefficient y of corresponding image points
                        self.IP_1_int[counter, 0]     = int(self.x_low[i, j])
                        self.IP_1_int[counter, 1]     = int(self.y_low[i, j])
                        
                        self.IP_1_float[counter, 0]   = self.coef_x_ip[i, j]
                        self.IP_1_float[counter, 1]   = self.coef_y_ip[i, j]  
                        self.IP_1_float[counter, 2]   = self.d_ip_dynamic[i, j]                      
                        
                        #BI array consists of x_low, y_low, coefficient x, coefficient y, x location, and y location of corresponding body intercept points
                        self.BI_1_int[counter, 0]     = int(self.x_low_bi[i, j])
                        self.BI_1_int[counter, 1]     = int(self.y_low_bi[i, j])
                        
                        self.BI_1_float[counter, 0]   = self.coef_x_bi[i, j]
                        self.BI_1_float[counter, 1]   = self.coef_y_bi[i, j]
                        self.BI_1_float[counter, 2]   = self.x_bi[i,j]
                        self.BI_1_float[counter, 3]   = self.y_bi[i,j]
                        
                        counter          += 1 
            
            self.num_ibm_ghost_cell_second_layer = False
            self.GC_2_int = False
            self.GC_2_float = False
            self.IP_2_int = False
            self.IP_2_float = False
            
            self.num_ibm_ghost_cell_third_layer = False
            self.GC_3_int = False
            self.GC_3_float = False
            self.IP_3_int = False
            self.IP_3_float = False
            
            self.num_ibm_ghost_cell_fourth_layer = False
            self.GC_4_int = False
            self.GC_4_float = False
            self.IP_4_int = False
            self.IP_4_float = False
            
            self.num_ibm_ghost_cell_fifth_layer = False
            self.GC_5_int = False
            self.GC_5_float = False
            self.IP_5_int = False
            self.IP_5_float = False
        
        if self.lib == 'cupy': #some GPU B.S.
                
            self.cell_type   =  cp.asarray(self.cell_type)
                            
            self.OC          =  cp.asarray(self.OC)
            
            self.GC_int      =  cp.asarray(self.GC_int)
            self.GC_1_int    =  cp.asarray(self.GC_1_int)
            
            self.GC_float    =  cp.asarray(self.GC_float)
            self.GC_1_float  =  cp.asarray(self.GC_1_float)
            
            self.IP_int      =  cp.asarray(self.IP_int)
            self.IP_1_int    =  cp.asarray(self.IP_1_int)
            
            self.IP_float    =  cp.asarray(self.IP_float)
            self.IP_1_float  =  cp.asarray(self.IP_1_float)
            
            self.BI_int      =  cp.asarray(self.BI_int)
            self.BI_1_int    =  cp.asarray(self.BI_1_int)
            
            self.BI_float    =  cp.asarray(self.BI_float)
            self.BI_1_float  =  cp.asarray(self.BI_1_float)
            
            
            #print(self.num_ibm_ghost_cell_first_layer)
            
            #print(self.cell_type[self.IP_1_int[:, 0], self.IP_1_int[:, 1]])
            #print(self.cell_type[self.IP_1_int[:, 0]+1, self.IP_1_int[:, 1]])
            #print(self.cell_type[self.IP_1_int[:, 0]+1, self.IP_1_int[:, 1]+1])
            #print(self.cell_type[self.IP_1_int[:, 0], self.IP_1_int[:, 1]+1])
            #print(npcp.sum(self.cell_type[self.IP_1_int[:, 0] +1, self.IP_1_int[:, 1]+2]-1.0))
            #print(npcp.sum(self.cell_type[self.IP_1_int[:, 0] +2, self.IP_1_int[:, 1]+2]-1.0))
        # NACA 0012 Geometry Function    
    def NACA_0012_y_chopped(self, x):
        # compute the y-coordiantes of the airfoi.
        import numpy as npcp
        
        xc = 1.0 - self.N_cylinder_dx * self.dx  
        R =  self.NACA_0012_y(xc)
        # xc is cutoff point
        # R = radius cutoff
        if x < xc:
            t = 0.12
            y = 5.0 * t * ((0.2969 * npcp.sqrt(x)) - (0.1260 * x) - (0.3516 * x**2.0) + (0.2843 * x**3.0) - (0.1015 * x**4.0))
        # standard from NASA I think.
        
        else:
        
            if x < xc + R:
            
                y = npcp.sqrt(R**2.0 - (x - xc)**2.0)
            #circle
            else:
            
                y = 0.0
            # zero thickness
        return y

    def NACA_0012_y(self, x):
        
        import numpy as npcp
        
        t = 0.12
        y = 5.0 * t * ((0.2969 * npcp.sqrt(x)) - (0.1260 * x) - (0.3516 * x**2.0) + (0.2843 * x**3.0) - (0.1015 * x**4.0))
        # finds the standard
        return y           
    
        plot.figure()
        plot.imshow(np.rot90(self.cell_type), extent=[-1.0, 3.0, -0.5, 0.5])
        plot.colorbar()
        plot.savefig('/mnt/c/Users/akula/IBM/airfoil_cell_types.png',bbox_inches = 'tight', dpi=600)
        plot.close()     

