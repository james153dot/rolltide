import numpy as np

class PostProcess_Coefficients:
    def __init__(self, dir_plot_data, num_of_bi, Nx, Ny, Nz, x, y, z, x_sub_left, y_sub_bottom, z_sub_back, xc, yc, zc, d, d_ip, dx, dy, dz, c_mu, num_ghosts, P_initial, u_initial, rho_initial):
        self.num_of_bi = num_of_bi
        self.dir_plot_data = dir_plot_data
        
        self.Nx = Nx
        self.Ny = Ny 
        self.Nz = Nz
        
        self.x = x
        self.y = y
        self.z = z
        
        self.xc = xc
        self.yc = yc
        self.zc = zc
        
        self.x_sub_left = x_sub_left
        self.y_sub_bottom = y_sub_bottom
        self.z_sub_back = z_sub_back
        
        self.d = d
        self.d_ip = d_ip
        
        self.dx = dx
        self.dy = dy
        self.dz = dz
        
        self.num_ghosts = num_ghosts
        self.c_mu = c_mu 
        
        self.P_inf = P_initial
        self.vel_inf = u_initial
        self.rho_inf = rho_initial
        
    def SaveCoefficients(self, parameter, value):
        dir_parameter = self.dir_plot_data + '/' + parameter + '.txt'
        with open(dir_parameter, 'a+') as f:
            f.write(str(value) + '\n')
    
    def SaveCoefficientsArray(self, parameter, array):
        dir_parameter = self.dir_plot_data + '/' + parameter + '.txt'
        with open(dir_parameter, 'w') as f:
            for i in range(array.shape[0]):
                f.write(str(array[i]) + '\n')
    
    def InitializeCoefficientArrays(self, max_step):
        self.teta = np.linspace(0.0, 360.0, num=self.num_of_bi)

        self.norm_0 = np.zeros((self.num_of_bi, ), dtype=np.float64)
        self.norm_1 = np.zeros((self.num_of_bi, ), dtype=np.float64)
        self.norm_2 = np.zeros((self.num_of_bi, ), dtype=np.float64)

        self.x_bi = np.zeros((self.num_of_bi, ), dtype=np.float64)
        self.y_bi = np.zeros((self.num_of_bi, ), dtype=np.float64)
        self.z_bi = np.zeros((self.num_of_bi, ), dtype=np.float64)

        self.x_ip = np.zeros((self.num_of_bi, ), dtype=np.float64)
        self.y_ip = np.zeros((self.num_of_bi, ), dtype=np.float64)
        self.z_ip = np.zeros((self.num_of_bi, ), dtype=np.float64)

        self.x_gc = np.zeros((self.num_of_bi, ), dtype=np.float64)
        self.y_gc = np.zeros((self.num_of_bi, ), dtype=np.float64)
        self.z_gc = np.zeros((self.num_of_bi, ), dtype=np.float64)

        self.x_low_bi_index = np.zeros((self.num_of_bi, ), dtype=np.int64)
        self.y_low_bi_index = np.zeros((self.num_of_bi, ), dtype=np.int64)
        self.z_low_bi_index = np.zeros((self.num_of_bi, ), dtype=np.int64)

        self.x_low_ip_index = np.zeros((self.num_of_bi, ), dtype=np.int64)
        self.y_low_ip_index = np.zeros((self.num_of_bi, ), dtype=np.int64)
        self.z_low_ip_index = np.zeros((self.num_of_bi, ), dtype=np.int64)

        self.x_low_gc_index = np.zeros((self.num_of_bi, ), dtype=np.int64)
        self.y_low_gc_index = np.zeros((self.num_of_bi, ), dtype=np.int64)
        self.z_low_gc_index = np.zeros((self.num_of_bi, ), dtype=np.int64)

        self.coefficients_x_bi = np.zeros((self.num_of_bi, ), dtype=np.float64)
        self.coefficients_y_bi = np.zeros((self.num_of_bi, ), dtype=np.float64)
        self.coefficients_z_bi = np.zeros((self.num_of_bi, ), dtype=np.float64)

        self.coefficients_x_ip = np.zeros((self.num_of_bi, ), dtype=np.float64)
        self.coefficients_y_ip = np.zeros((self.num_of_bi, ), dtype=np.float64)
        self.coefficients_z_ip = np.zeros((self.num_of_bi, ), dtype=np.float64)

        self.coefficients_x_gc = np.zeros((self.num_of_bi, ), dtype=np.float64)
        self.coefficients_y_gc = np.zeros((self.num_of_bi, ), dtype=np.float64)
        self.coefficients_z_gc = np.zeros((self.num_of_bi, ), dtype=np.float64)

        self.f1_vel_ip = np.zeros((3, self.num_of_bi), dtype=np.float64)
        self.f2_vel_ip = np.zeros((3, self.num_of_bi), dtype=np.float64)
        self.V_ip = np.zeros((3, self.num_of_bi), dtype=np.float64)

        self.T = np.zeros((self.num_of_bi), dtype=np.float64)
        self.tau = np.zeros((self.num_of_bi), dtype=np.float64)

        self.f1_pres = np.zeros((self.num_of_bi), dtype=np.float64)
        self.f2_pres = np.zeros((self.num_of_bi), dtype=np.float64)
        self.pressure_bi = np.zeros((self.num_of_bi), dtype=np.float64)
        self.pressure_ip = np.zeros((self.num_of_bi), dtype=np.float64)

        self.f1_density = np.zeros((self.num_of_bi), dtype=np.float64)
        self.f2_density = np.zeros((self.num_of_bi), dtype=np.float64)
        self.density_bi = np.zeros((self.num_of_bi), dtype=np.float64)
        self.density_ip = np.zeros((self.num_of_bi), dtype=np.float64)

        # local current coefficients
        self.C_P = np.zeros((self.num_of_bi,), dtype=np.float64)
        self.C_F = np.zeros((self.num_of_bi,), dtype=np.float64)
        self.C_D = np.zeros((self.num_of_bi,), dtype=np.float64)
        self.C_D_pressure = np.zeros((self.num_of_bi,), dtype=np.float64)
        self.C_D_viscous = np.zeros((self.num_of_bi,), dtype=np.float64)
        self.C_L = np.zeros((self.num_of_bi,), dtype=np.float64)
        self.Nu = np.zeros((self.num_of_bi,), dtype=np.float64)
        
        self.base_pres_coef = np.zeros((max_step, ), dtype=np.float64)
        self.pres_coef = np.zeros((max_step, ), dtype=np.float64)
        self.drag_coef = np.zeros((max_step, ), dtype=np.float64)
        self.lift_coef = np.zeros((max_step, ), dtype=np.float64)
        self.skin_coef = np.zeros((max_step, ), dtype=np.float64)
        self.nu_coef = np.zeros((max_step, ), dtype=np.float64)

        
        '''
        dir_base_pres  = self.dir_plot_data + '/base_pres_coef.txt'      
        if os.path.exists(dir_base_pres):
            os.remove(dir_base_pres)
        
        dir_pres  = self.dir_plot_data + '/pres_coef.txt'      
        if os.path.exists(dir_pres):
            os.remove(dir_pres)
            
        dir_skin  = self.dir_plot_data + '/skin_coef.txt'      
        if os.path.exists(dir_skin):
            os.remove(dir_skin)

        dir_drag  = self.dir_plot_data + '/drag_coef.txt'      
        if os.path.exists(dir_drag):
            os.remove(dir_drag)
        
        dir_lift  = self.dir_plot_data + '/lift_coef.txt'      
        if os.path.exists(dir_lift):
            os.remove(dir_lift)
        
        dir_nusselt  = self.dir_plot_data + '/nu_coef.txt'      
        if os.path.exists(dir_nusselt):
            os.remove(dir_nusselt)
        '''
        
    ## Filling the solid cells that should be filled for second order cases coefficient calculations
    def fill_Solids_for_Pressure(self, V, cell_type): #SKIPPED
    
        for i in range(self.Nx):
            
            for j in range(self.Ny):
                
                if cell_type[i, j] == 2:
                    
                    if (cell_type[i, j+1] == 1 and cell_type[i+1,j]== 1) or (cell_type[i, j+1]== 1 and cell_type[i-1,j]== 1) or (cell_type[i-1, j]== 1 and cell_type[i,j-1]== 1) or (cell_type[i+1, j]== 1 and cell_type[i,j-1]== 1):
                    
                        distance_to_center_solid  =  np.sqrt((self.x[i]-self.xc)**2.0 +(self.y[j]-self.yc)**2.0)
                        norm_0_solid              = (self.x[i] - self.xc) /  distance_to_center_solid
                        norm_1_solid              = (self.y[j] - self.yc) /  distance_to_center_solid
                        
                        if self.x[i] > self.xc:
                            
                            x_ip_solid  = self.xc + np.sqrt((self.d/2.0  + self.d_ip)**2.0 - ((self.d/2.0  + self.d_ip)*norm_1_solid)**2.0) 
                        
                        else:
                            
                            x_ip_solid  = self.xc - np.sqrt((self.d/2.0  + self.d_ip)**2.0 - ((self.d/2.0  + self.d_ip)*norm_1_solid)**2.0) 
                        
                        if self.y[j] > self.yc:
                            
                            y_ip_solid  = self.yc + np.sqrt((self.d/2.0  + self.d_ip )**2.0 - ((self.d/2.0  + self.d_ip)*norm_0_solid)**2.0)
                        
                        else:
                            
                            y_ip_solid  = self.yc - np.sqrt((self.d/2.0  + self.d_ip )**2.0 - ((self.d/2.0  + self.d_ip)*norm_0_solid)**2.0)
                        
                        x_low_ip_solid   = (x_ip_solid - self.dx/2.0 - (self.xc - self.x_sub_left)) // self.dx 
                        coef_x_ip_solid  = (x_ip_solid - self.x[int(x_low_ip_solid)])  / (self.dx)
                        y_low_ip_solid   = (y_ip_solid - self.dy/2.0 - (self.yc - self.y_sub_bottom)) // self.dy 
                        coef_y_ip_solid  = (y_ip_solid - self.y[int(y_low_ip_solid)])  / (self.dy)    

                        P_ip_f1    = (1.0 - coef_x_ip_solid) * V[3, int(x_low_ip_solid)+self.num_ghosts, int(y_low_ip_solid)+self.num_ghosts]     + coef_x_ip_solid * V[3, int(x_low_ip_solid) +self.num_ghosts+1, int(y_low_ip_solid)+self.num_ghosts ]
                        P_ip_f2    = (1.0 - coef_x_ip_solid) * V[3, int(x_low_ip_solid)+self.num_ghosts, int(y_low_ip_solid)+self.num_ghosts +1]  + coef_x_ip_solid * V[3, int(x_low_ip_solid) +self.num_ghosts+1, int(y_low_ip_solid) +self.num_ghosts+ 1]
                        P_ip       = (1.0 - coef_y_ip_solid) * P_ip_f1 + coef_y_ip_solid * P_ip_f2
                                            
                        V[3, i+self.num_ghosts, j+self.num_ghosts] = P_ip
        return V
    
    ### IGNORE ABOVE

class SphereInterpolation:
    def __init__(self, num_of_bi, radius, d_ip, xc, yc, zc, dx, dy, dz, x_sub_left, y_sub_bottom, z_sub_back, teta, phi, num_ghosts, c_mu):
        self.num_of_bi = num_of_bi
        self.radius = radius
        self.d_ip = d_ip
        self.xc = xc
        self.yc = yc
        self.zc = zc
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.x_sub_left = x_sub_left
        self.y_sub_bottom = y_sub_bottom
        self.z_sub_back = z_sub_back
        self.teta = teta
        self.phi = phi
        self.num_ghosts = num_ghosts
        self.c_mu = c_mu

        self.norm_0 = np.zeros(num_of_bi)
        self.norm_1 = np.zeros(num_of_bi)
        self.norm_2 = np.zeros(num_of_bi)

        self.x_bi = np.zeros(num_of_bi)
        self.y_bi = np.zeros(num_of_bi)
        self.z_bi = np.zeros(num_of_bi)

        self.x_ip = np.zeros(num_of_bi)
        self.y_ip = np.zeros(num_of_bi)
        self.z_ip = np.zeros(num_of_bi)

        self.x_gc = np.zeros(num_of_bi)
        self.y_gc = np.zeros(num_of_bi)
        self.z_gc = np.zeros(num_of_bi)

        self.x_low_bi_index = np.zeros(num_of_bi, dtype=int)
        self.y_low_bi_index = np.zeros(num_of_bi, dtype=int)
        self.z_low_bi_index = np.zeros(num_of_bi, dtype=int)

        self.x_low_ip_index = np.zeros(num_of_bi, dtype=int)
        self.y_low_ip_index = np.zeros(num_of_bi, dtype=int)
        self.z_low_ip_index = np.zeros(num_of_bi, dtype=int)

        self.x_low_gc_index = np.zeros(num_of_bi, dtype=int)
        self.y_low_gc_index = np.zeros(num_of_bi, dtype=int)
        self.z_low_gc_index = np.zeros(num_of_bi, dtype=int)

        self.coefficients_x_bi = np.zeros(num_of_bi)
        self.coefficients_y_bi = np.zeros(num_of_bi)
        self.coefficients_z_bi = np.zeros(num_of_bi)

        self.coefficients_x_ip = np.zeros(num_of_bi)
        self.coefficients_y_ip = np.zeros(num_of_bi)
        self.coefficients_z_ip = np.zeros(num_of_bi)

        self.coefficients_x_gc = np.zeros(num_of_bi)
        self.coefficients_y_gc = np.zeros(num_of_bi)
        self.coefficients_z_gc = np.zeros(num_of_bi)

        self.V_ip = np.zeros((3, num_of_bi))
        self.tau = np.zeros(num_of_bi)
        self.pressure_bi = np.zeros(num_of_bi)
        self.density_bi = np.zeros(num_of_bi)
        self.pressure_ip = np.zeros(num_of_bi)
        self.density_ip = np.zeros(num_of_bi)

    def get_normals(self, theta, phi):
        theta_rad = np.deg2rad(theta)
        phi_rad = np.deg2rad(phi)

        norm_0 = np.sin(theta_rad) * np.cos(phi_rad)
        norm_1 = np.sin(theta_rad) * np.sin(phi_rad)
        norm_2 = np.cos(theta_rad)

        return norm_0, norm_1, norm_2

    def calculate_coordinates(self, index, norm_0, norm_1, norm_2):
        radius_bi = self.radius
        radius_ip = self.radius + self.d_ip
        radius_gc = self.radius - 0.5 * self.d_ip

        self.x_bi[index] = self.xc + radius_bi * norm_0
        self.y_bi[index] = self.yc + radius_bi * norm_1
        self.z_bi[index] = self.zc + radius_bi * norm_2

        self.x_ip[index] = self.xc + radius_ip * norm_0
        self.y_ip[index] = self.yc + radius_ip * norm_1
        self.z_ip[index] = self.zc + radius_ip * norm_2

        self.x_gc[index] = self.xc + radius_gc * norm_0
        self.y_gc[index] = self.yc + radius_gc * norm_1
        self.z_gc[index] = self.zc + radius_gc * norm_2

    def calculate_indices_and_coefficients(self, index):
        self.x_low_bi_index[index] = int((self.x_bi[index] - self.dx / 2.0 - (self.xc - self.x_sub_left)) // self.dx)
        self.y_low_bi_index[index] = int((self.y_bi[index] - self.dy / 2.0 - (self.yc - self.y_sub_bottom)) // self.dy)
        self.z_low_bi_index[index] = int((self.z_bi[index] - self.dz / 2.0 - (self.zc - self.z_sub_back)) // self.dz)

        self.x_low_ip_index[index] = int((self.x_ip[index] - self.dx / 2.0 - (self.xc - self.x_sub_left)) // self.dx)
        self.y_low_ip_index[index] = int((self.y_ip[index] - self.dy / 2.0 - (self.yc - self.y_sub_bottom)) // self.dy)
        self.z_low_ip_index[index] = int((self.z_ip[index] - self.dz / 2.0 - (self.zc - self.z_sub_back)) // self.dz)

        self.x_low_gc_index[index] = int((self.x_gc[index] - self.dx / 2.0 - (self.xc - self.x_sub_left)) // self.dx)
        self.y_low_gc_index[index] = int((self.y_gc[index] - self.dy / 2.0 - (self.yc - self.y_sub_bottom)) // self.dy)
        self.z_low_gc_index[index] = int((self.z_gc[index] - self.dz / 2.0 - (self.zc - self.z_sub_back)) // self.dz)

        self.coefficients_x_bi[index] = (self.x_bi[index] - self.x[self.x_low_bi_index[index]]) / self.dx
        self.coefficients_y_bi[index] = (self.y_bi[index] - self.y[self.y_low_bi_index[index]]) / self.dy
        self.coefficients_z_bi[index] = (self.z_bi[index] - self.z[self.z_low_bi_index[index]]) / self.dz

        self.coefficients_x_ip[index] = (self.x_ip[index] - self.x[self.x_low_ip_index[index]]) / self.dx
        self.coefficients_y_ip[index] = (self.y_ip[index] - self.y[self.y_low_ip_index[index]]) / self.dy
        self.coefficients_z_ip[index] = (self.z_ip[index] - self.z[self.z_low_ip_index[index]]) / self.dz

        self.coefficients_x_gc[index] = (self.x_gc[index] - self.x[self.x_low_gc_index[index]]) / self.dx
        self.coefficients_y_gc[index] = (self.y_gc[index] - self.y[self.y_low_gc_index[index]]) / self.dy
        self.coefficients_z_gc[index] = (self.z_gc[index] - self.z[self.z_low_gc_index[index]]) / self.dz

    def trilinear_interpolation(self, x_idx, y_idx, z_idx, coeffs, V, component):
        cx, cy, cz = coeffs

        c00 = V[component, x_idx, y_idx, z_idx] * (1 - cx) + V[component, x_idx + 1, y_idx, z_idx] * cx
        c01 = V[component, x_idx, y_idx, z_idx + 1] * (1 - cx) + V[component, x_idx + 1, y_idx, z_idx + 1] * cx
        c10 = V[component, x_idx, y_idx + 1, z_idx] * (1 - cx) + V[component, x_idx + 1, y_idx + 1, z_idx] * cx
        c11 = V[component, x_idx, y_idx + 1, z_idx + 1] * (1 - cx) + V[component, x_idx + 1, y_idx + 1, z_idx + 1] * cx

        c0 = c00 * (1 - cy) + c10 * cy
        c1 = c01 * (1 - cy) + c11 * cy

        c = c0 * (1 - cz) + c1 * cz

        return c

    def get_tau_constant_mu(self, V):
        for i in range(self.num_of_bi):
            x_idx, y_idx, z_idx = self.x_low_ip_index[i], self.y_low_ip_index[i], self.z_low_ip_index[i]
            coeffs = (self.coefficients_x_ip[i], self.coefficients_y_ip[i], self.coefficients_z_ip[i])

            for comp in range(3):
                self.V_ip[comp, i] = self.trilinear_interpolation(x_idx, y_idx, z_idx, coeffs, V, comp)

            self.tau[i] = (self.c_mu / self.d_ip) * (-self.V_ip[0, i] * self.norm_1[i] + self.V_ip[1, i] * self.norm_0[i])

            if self.teta[i] > 0.0 and self.teta[i] < 180.0:
                self.tau[i] = np.negative(self.tau[i])

        return self.tau

    def get_tau_power_law_mu_adiabatic(self, V, T_initial, R_specific, power):
        density_bi = self.get_density_bi(V)
        pressure_bi = self.get_P_bi(V)

        T_bi = pressure_bi / (density_bi * R_specific)

        for i in range(self.num_of_bi):
            x_idx, y_idx, z_idx = self.x_low_ip_index[i], self.y_low_ip_index[i], self.z_low_ip_index[i]
            coeffs = (self.coefficients_x_ip[i], self.coefficients_y_ip[i], self.coefficients_z_ip[i])

            for comp in range(3):
                self.V_ip[comp, i] = self.trilinear_interpolation(x_idx, y_idx, z_idx, coeffs, V, comp)

            self.tau[i] = ((self.c_mu * (T_bi[i] / T_initial)**power) / self.d_ip) * (-self.V_ip[0, i] * self.norm_1[i] + self.V_ip[1, i] * self.norm_0[i])

            if self.teta[i] > 0.0 and self.teta[i] < 180.0:
                self.tau[i] = np.negative(self.tau[i])

        return self.tau

    def get_tau_power_law_mu_isothermal(self, V, TR, R_specific, power):
        for i in range(self.num_of_bi):
            x_idx, y_idx, z_idx = self.x_low_ip_index[i], self.y_low_ip_index[i], self.z_low_ip_index[i]
            coeffs = (self.coefficients_x_ip[i], self.coefficients_y_ip[i], self.coefficients_z_ip[i])

            for comp in range(3):
                self.V_ip[comp, i] = self.trilinear_interpolation(x_idx, y_idx, z_idx, coeffs, V, comp)

            self.tau[i] = ((self.c_mu * TR**power) / self.d_ip) * (-self.V_ip[0, i] * self.norm_1[i] + self.V_ip[1, i] * self.norm_0[i])

            if self.teta[i] > 0.0 and self.teta[i] < 180.0:
                self.tau[i] = np.negative(self.tau[i])

        return self.tau

    def get_P_bi(self, V):
        for i in range(self.num_of_bi):
            x_idx, y_idx, z_idx = self.x_low_bi_index[i], self.y_low_bi_index[i], self.z_low_bi_index[i]
            coeffs = (self.coefficients_x_bi[i], self.coefficients_y_bi[i], self.coefficients_z_bi[i])

            self.pressure_bi[i] = self.trilinear_interpolation(x_idx, y_idx, z_idx, coeffs, V, 3)

        return self.pressure_bi

    def get_density_bi(self, V):
        for i in range(self.num_of_bi):
            x_idx, y_idx, z_idx = self.x_low_bi_index[i], self.y_low_bi_index[i], self.z_low_bi_index[i]
            coeffs = (self.coefficients_x_bi[i], self.coefficients_y_bi[i], self.coefficients_z_bi[i])

            self.density_bi[i] = self.trilinear_interpolation(x_idx, y_idx, z_idx, coeffs, V, 0)

        return self.density_bi

    def get_P_ip(self, V):
        for i in range(self.num_of_bi):
            x_idx, y_idx, z_idx = self.x_low_ip_index[i], self.y_low_ip_index[i], self.z_low_ip_index[i]
            coeffs = (self.coefficients_x_ip[i], self.coefficients_y_ip[i], self.coefficients_z_ip[i])

            self.pressure_ip[i] = self.trilinear_interpolation(x_idx, y_idx, z_idx, coeffs, V, 3)

        return self.pressure_ip

    def get_density_ip(self, V):
        for i in range(self.num_of_bi):
            x_idx, y_idx, z_idx = self.x_low_ip_index[i], self.y_low_ip_index[i], self.z_low_ip_index[i]
            coeffs = (self.coefficients_x_ip[i], self.coefficients_y_ip[i], self.coefficients_z_ip[i])

            self.density_ip[i] = self.trilinear_interpolation(x_idx, y_idx, z_idx, coeffs, V, 0)

        return self.density_ip

    def get_Coefficients(self, pressure, tau):
        for i in range(self.num_of_bi):
            self.C_P[i] = (pressure[i] - self.P_inf) / (0.5 * self.rho_inf * (self.vel_inf**2.0))
            self.C_F[i] = tau[i] / (0.5 * self.rho_inf * (self.vel_inf**2.0))

            if self.teta[i] >= 0.0 and self.teta[i] <= 180.0:
                self.C_D[i] = np.pi * ((-pressure[i] * self.norm_0[i] + tau[i] * self.norm_1[i]) / (0.5 * (self.num_of_bi-1) * self.rho_inf * (self.vel_inf ** 2.0)))
                self.C_D_pressure[i] = np.pi * ((-pressure[i] * self.norm_0[i]) / (0.5 * (self.num_of_bi-1) * self.rho_inf * (self.vel_inf ** 2.0)))
                self.C_D_viscous[i] = np.pi * ((tau[i] * self.norm_1[i]) / (0.5 * (self.num_of_bi-1) * self.rho_inf * (self.vel_inf ** 2.0)))
                self.C_L[i] = np.pi * ((-pressure[i] * self.norm_1[i] - tau[i] * self.norm_0[i]) / (0.5 * (self.num_of_bi-1) * self.rho_inf * (self.vel_inf ** 2.0)))
            else:
                self.C_D[i] = np.pi * ((-pressure[i] * self.norm_0[i] - tau[i] * self.norm_1[i]) / (0.5 * (self.num_of_bi-1) * self.rho_inf * (self.vel_inf ** 2.0)))
                self.C_D_pressure[i] = np.pi * ((-pressure[i] * self.norm_0[i]) / (0.5 * (self.num_of_bi-1) * self.rho_inf * (self.vel_inf ** 2.0)))
                self.C_D_viscous[i] = -np.pi * ((tau[i] * self.norm_1[i]) / (0.5 * (self.num_of_bi-1) * self.rho_inf * (self.vel_inf ** 2.0)))
                self.C_L[i] = np.pi * ((-pressure[i] * self.norm_1[i] + tau[i] * self.norm_0[i]) / (0.5 * (self.num_of_bi-1) * self.rho_inf * (self.vel_inf ** 2.0)))

        return self.C_P, self.C_F, self.C_D, self.C_D_pressure, self.C_D_viscous, self.C_L

    def getTemperatureGradientIB(self, temp_bi, temp_ip):
        return (temp_ip - temp_bi) / (self.d_ip)
