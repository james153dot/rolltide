#include "flow/flow_models/single-species/FlowModelImmersedBoundaryMethodSingleSpecies.hpp"

FlowModelImmersedBoundaryMethodSingleSpecies::FlowModelImmersedBoundaryMethodSingleSpecies(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const int& num_eqn,
    const HAMERS_SHARED_PTR<ImmersedBoundaries>& immersed_boundaries,
    const HAMERS_SHARED_PTR<tbox::Database>& immersed_boundary_method_db,
    const HAMERS_SHARED_PTR<EquationOfStateMixingRules>& equation_of_state_mixing_rules):
        FlowModelImmersedBoundaryMethod(
            object_name,
            dim,
            grid_geometry,
            num_species,
            num_eqn,
            immersed_boundaries,
            immersed_boundary_method_db,
            equation_of_state_mixing_rules)
{
    /*
     * Read the body density.
     */
    
    if (immersed_boundary_method_db->keyExists("body_density"))
    {
        d_rho_body = immersed_boundary_method_db->getReal("body_density");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethodSingleSpecies::FlowModelImmersedBoundaryMethodSingleSpecies()\n"
            << "Required 'body_density' entry from input database missing."
            << std::endl);
    }
    
    /*
     * Read the body velocity.
     */
    
    if (immersed_boundary_method_db->keyExists("body_velocity"))
    {
        d_vel_body = immersed_boundary_method_db->getRealVector("body_velocity");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethodSingleSpecies::FlowModelImmersedBoundaryMethodSingleSpecies()\n"
            << "Required 'body_velocity' entry from input database missing."
            << std::endl);
    }
    
    if (static_cast<int>(d_vel_body.size()) != d_dim.getValue())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethodSingleSpecies::FlowModelImmersedBoundaryMethodSingleSpecies()\n"
            << "Size of 'body_velocity' entry from input database does not match problem dimension."
            << std::endl);
    }
    
    /*
     * Read the body pressure.
     */
    
    if (immersed_boundary_method_db->keyExists("body_pressure"))
    {
        d_p_body = immersed_boundary_method_db->getReal("body_pressure");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethodSingleSpecies::FlowModelImmersedBoundaryMethodSingleSpecies()\n"
            << "Required 'body_pressure' entry from input database missing."
            << std::endl);
    }

    /*
     * Read the body temperature.
     */
    
    if (immersed_boundary_method_db->keyExists("body_temperature"))
    {
        d_temp_body = immersed_boundary_method_db->getReal("body_temperature");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethodSingleSpecies::FlowModelImmersedBoundaryMethodSingleSpecies()\n"
            << "Required 'body_temperature' entry from input database missing."
            << std::endl);
    }
    
    std::vector<Real*> thermo_properties_ptr;
    std::vector<const Real*> thermo_properties_const_ptr;
    
    const int num_thermo_properties = d_equation_of_state_mixing_rules->
        getNumberOfSpeciesThermodynamicProperties();
    
    thermo_properties_ptr.reserve(num_thermo_properties);
    thermo_properties_const_ptr.reserve(num_thermo_properties);
    d_thermo_properties.resize(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
        thermo_properties_const_ptr.push_back(&d_thermo_properties[ti]);
    }
    
    d_equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
        thermo_properties_ptr,
        0);
    
    /*
     * Compute the values of conservative variables inside the body.
     */
    
    const Real T_body = d_equation_of_state_mixing_rules->getEquationOfState()->
        getTemperature(
            &d_rho_body,
            &d_p_body,
            thermo_properties_const_ptr);
    
    const Real epsilon_body = d_equation_of_state_mixing_rules->getEquationOfState()->
        getInternalEnergyFromTemperature(
            &d_rho_body,
            &T_body,
            thermo_properties_const_ptr);
    
    d_mom_body.resize(d_dim.getValue());
    
    if (d_dim == tbox::Dimension(1))
    {
        const Real& u_body = d_vel_body[0];
        
        d_mom_body[0] = d_rho_body*u_body;;
        
        d_E_body = d_rho_body*(epsilon_body + Real(1)/Real(2)*
            (u_body*u_body));
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const Real& u_body = d_vel_body[0];
        const Real& v_body = d_vel_body[1];
        
        d_mom_body[0] = d_rho_body*u_body;
        d_mom_body[1] = d_rho_body*v_body;
        
        d_E_body = d_rho_body*(epsilon_body + Real(1)/Real(2)*
            (u_body*u_body + v_body*v_body));
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const Real& u_body = d_vel_body[0];
        const Real& v_body = d_vel_body[1];
        const Real& w_body = d_vel_body[2];
        
        d_mom_body[0] = d_rho_body*u_body;
        d_mom_body[1] = d_rho_body*v_body;
        d_mom_body[2] = d_rho_body*w_body;
        
        d_E_body = d_rho_body*(epsilon_body + Real(1)/Real(2)*
            (u_body*u_body + v_body*v_body + w_body*w_body));
    }
}


/*
 * Set the immersed boundary method ghost cells for the cell data of conservative variables.
 */
void FlowModelImmersedBoundaryMethodSingleSpecies::setConservativeVariablesCellDataImmersedBoundaryGhosts(
    const hier::Patch& patch,
    const double data_time,
    const bool initial_time,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_var_data,
    const HAMERS_SHARED_PTR<pdat::CellData<int> >& data_mask,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_wall_distance,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_surface_normal,
    const hier::IntVector& offset_cons_var,
    const hier::IntVector& offset_IB,
    const hier::IntVector& ghostcell_dims_cons_var,
    const hier::IntVector& ghostcell_dims_IB,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims)
{
    NULL_USE(data_time);
    NULL_USE(initial_time);
    
    /*
     * Get the grid spacings and the lower coordinates.
     */
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    const double* const dx = patch_geom->getDx();
    const double* const patch_xlo = patch_geom->getXLower();
    
    /*
     * Get the data of the conservative variables.
     */
    
    const HAMERS_SHARED_PTR<pdat::CellData<Real> > data_density      = conservative_var_data[0];
    const HAMERS_SHARED_PTR<pdat::CellData<Real> > data_momentum     = conservative_var_data[1];
    const HAMERS_SHARED_PTR<pdat::CellData<Real> > data_total_energy = conservative_var_data[2];
    
    /*
     * Get the pointers to the data.
     */
    
    Real* rho = data_density->getPointer(0);
    Real* E   = data_total_energy->getPointer(0);
    
    int* mask = data_mask->getPointer(0);
    Real* dist = data_wall_distance->getPointer(0);
    
    const Real& rho_body  = d_rho_body;
    const Real& E_body    = d_E_body;
    const Real& T_body    = d_temp_body;

    
    // Get the thermodynamic properties of the species.
    std::vector<const Real*> thermo_properties_ptr;
    thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
    for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
    }
    
    const Real one  = Real(1);
    const Real half = Real(1)/Real(2);
    
    // Distance from cylinder boundary to the image point sqrt(2 + epsilon), isotropic grid cells are assuemd.
    const Real d_ip = sqrt(Real(2))*Real(dx[0]) + HAMERS_REAL_EPSILON;
    
    if (d_dim == tbox::Dimension(1))
    {
        const Real& rho_u_body = d_mom_body[0];
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_cons_var = offset_cons_var[0];
        const int offset_0_IB = offset_IB[0];
        
        // Get the pointers to the data.
        Real* rho_u = data_momentum->getPointer(0);
        
        // Real* norm_0 = data_surface_normal->getPointer(0);
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_cons_var = i + offset_0_cons_var;
            const int idx_IB = i + offset_0_IB;
            
            if (mask[idx_IB] == int(IB_MASK::IB_GHOST))
            {
                // NEED TO BE CHANGED!!!
                rho[idx_cons_var]   = rho_body;
                rho_u[idx_cons_var] = rho_u_body;
                E[idx_cons_var]     = E_body;
            }
            else if (mask[idx_IB] == int(IB_MASK::BODY))
            {
                rho[idx_cons_var]   = rho_body;
                rho_u[idx_cons_var] = rho_u_body;
                E[idx_cons_var]     = E_body;
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const Real& rho_u_body = d_mom_body[0];
        const Real& rho_v_body = d_mom_body[1];
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int offset_0_cons_var = offset_cons_var[0];
        const int offset_1_cons_var = offset_cons_var[1];
        const int ghostcell_dim_0_cons_var = ghostcell_dims_cons_var[0];
        
        const int offset_0_IB = offset_IB[0];
        const int offset_1_IB = offset_IB[1];
        const int ghostcell_dim_0_IB = ghostcell_dims_IB[0];
        
        // Get the pointers to the data.
        Real* rho_u = data_momentum->getPointer(0);
        Real* rho_v = data_momentum->getPointer(1);
        
        Real* norm_0 = data_surface_normal->getPointer(0);
        Real* norm_1 = data_surface_normal->getPointer(1);
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            // HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_cons_var = (i + offset_0_cons_var) +
                    (j + offset_1_cons_var)*ghostcell_dim_0_cons_var;
                
                const int idx_IB = (i + offset_0_IB) +
                    (j + offset_1_IB)*ghostcell_dim_0_IB;
                
                // Compute the coordinates.
                Real x[2];
                x[0] = Real(patch_xlo[0]) + (Real(i) + half)*Real(dx[0]); // local x coordinates
                x[1] = Real(patch_xlo[1]) + (Real(j) + half)*Real(dx[1]); // local y coordinates
                
                if (mask[idx_IB] == int(IB_MASK::IB_GHOST))
                {
                    // General computation of image point x and y locations.
                    const Real x_ip = x[0] + (dist[idx_IB] + d_ip)*norm_0[idx_IB];
                    const Real y_ip = x[1] + (dist[idx_IB] + d_ip)*norm_1[idx_IB];
                    
                    // x and indices of bottom left in the bilinear interpolation stencil.
                    const int ip_location_index_0 = int(floor((x_ip - patch_xlo[0] - half * Real(dx[0]))/Real(dx[0])));
                    const int ip_location_index_1 = int(floor((y_ip - patch_xlo[1] - half * Real(dx[1]))/Real(dx[1])));
                    
                    const Real x_ip_BL = Real(patch_xlo[0]) + (Real(ip_location_index_0) + half)*Real(dx[0]);
                    const Real y_ip_BL = Real(patch_xlo[1]) + (Real(ip_location_index_1) + half)*Real(dx[1]);
                    
                    const Real ip_ratio_0 = (x_ip - x_ip_BL)/Real(dx[0]);
                    const Real ip_ratio_1 = (y_ip - y_ip_BL)/Real(dx[1]);
                    
                    // AFK declaration of the indexes of the cells in interpolation
                    const int idx_cons_var_BL  = (ip_location_index_0 + offset_0_cons_var) +
                                                 (ip_location_index_1 + offset_1_cons_var) * ghostcell_dim_0_cons_var;
                    
                    const int idx_cons_var_BR  = idx_cons_var_BL + 1;
                    const int idx_cons_var_TL  = idx_cons_var_BL + ghostcell_dim_0_cons_var;
                    const int idx_cons_var_TR  = idx_cons_var_BL + ghostcell_dim_0_cons_var + 1;
                    
                    // Bilinear interpolation to find image point x velocity value.
                    const Real u_IP_BL = rho_u[idx_cons_var_BL]/rho[idx_cons_var_BL];
                    const Real u_IP_BR = rho_u[idx_cons_var_BR]/rho[idx_cons_var_BR];
                    const Real u_IP_TL = rho_u[idx_cons_var_TL]/rho[idx_cons_var_TL];
                    const Real u_IP_TR = rho_u[idx_cons_var_TR]/rho[idx_cons_var_TR];
                    
                    // const Real u_body = rho_u_body / rho[idx_cons_var];
                    
                    const Real u_f1 = (one - ip_ratio_0) * u_IP_BL + ip_ratio_0 * u_IP_BR;
                    const Real u_f2 = (one - ip_ratio_0) * u_IP_TL + ip_ratio_0 * u_IP_TR;
                    const Real u_ip = (one - ip_ratio_1) * u_f1    + ip_ratio_1 * u_f2;
                    
                    // Bilinear interpolation to find image point x velocity value.
                    const Real v_IP_BL = rho_v[idx_cons_var_BL]/rho[idx_cons_var_BL];
                    const Real v_IP_BR = rho_v[idx_cons_var_BR]/rho[idx_cons_var_BR];
                    const Real v_IP_TL = rho_v[idx_cons_var_TL]/rho[idx_cons_var_TL];
                    const Real v_IP_TR = rho_v[idx_cons_var_TR]/rho[idx_cons_var_TR];
                    
                    // const Real v_body = rho_v_body / rho[idx_cons_var];
                    
                    const Real v_f1 = (one - ip_ratio_0)*v_IP_BL + ip_ratio_0*v_IP_BR; 
                    const Real v_f2 = (one - ip_ratio_0)*v_IP_TL + ip_ratio_0*v_IP_TR;
                    const Real v_ip = (one - ip_ratio_1)*v_f1    + ip_ratio_1*v_f2;
                    
                    // Compute velocities at the ghost cell for Euler equation (no-slip for normal direction and slip for tangential direction)
                    Real u_gc = Real(0);  // u velocity of the ghost cell
                    Real v_gc = Real(0);  // v velocity of the ghost cell
                    Real vel_ip_n = Real(0);
                    Real vel_ip_t = Real(0);
                    Real vel_gc_n = Real(0);
                    Real vel_gc_t = Real(0);
                    
                    if(d_bc_type_velocity == VELOCITY_IBC::SLIP)
                    {
                        vel_ip_n =  u_ip*norm_0[idx_IB] + v_ip*norm_1[idx_IB];
                        vel_ip_t = -u_ip*norm_1[idx_IB] + v_ip*norm_0[idx_IB];
                        
                        vel_gc_n = vel_ip_n - ((d_ip + dist[idx_IB]) / d_ip) * (vel_ip_n);
                        vel_gc_t = vel_ip_t;
                        
                        u_gc = vel_gc_n*norm_0[idx_IB] - vel_gc_t*norm_1[idx_IB];
                        v_gc = vel_gc_n*norm_1[idx_IB] + vel_gc_t*norm_0[idx_IB]; 
                    }
                    else if(d_bc_type_velocity == VELOCITY_IBC::NO_SLIP)
                    {
                        u_gc = u_ip - ((d_ip + dist[idx_IB])/d_ip) * (u_ip);
                        v_gc = v_ip - ((d_ip + dist[idx_IB])/d_ip) * (v_ip);
                    }
                        
                    // Bilinear interpolation to find image point energy.
                    const Real epsilon_BL = E[idx_cons_var_BL]/rho[idx_cons_var_BL] - half*(u_IP_BL*u_IP_BL + v_IP_BL*v_IP_BL);
                    const Real epsilon_TL = E[idx_cons_var_TL]/rho[idx_cons_var_TL] - half*(u_IP_TL*u_IP_TL + v_IP_TL*v_IP_TL);
                    const Real epsilon_BR = E[idx_cons_var_BR]/rho[idx_cons_var_BR] - half*(u_IP_BR*u_IP_BR + v_IP_BR*v_IP_BR);
                    const Real epsilon_TR = E[idx_cons_var_TR]/rho[idx_cons_var_TR] - half*(u_IP_TR*u_IP_TR + v_IP_TR*v_IP_TR);
                    
                    // Compute the pressure values in the stencils.
                    const Real p_BL = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getPressure(
                            &rho[idx_cons_var_BL],
                            &epsilon_BL,
                            thermo_properties_ptr);
                    
                    const Real p_TL = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getPressure(
                            &rho[idx_cons_var_TL],
                            &epsilon_TL,
                            thermo_properties_ptr);
                    
                    const Real p_BR = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getPressure(
                            &rho[idx_cons_var_BR],
                            &epsilon_BR,
                            thermo_properties_ptr);
                    
                    const Real p_TR = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getPressure(
                            &rho[idx_cons_var_TR],
                            &epsilon_TR,
                            thermo_properties_ptr);
                    
                    // Compute the temperature values in the stencils.
                    const Real T_BL = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getTemperature(
                            &rho[idx_cons_var_BL],
                            &p_BL,
                            thermo_properties_ptr);
                    
                    const Real T_TL = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getTemperature(
                            &rho[idx_cons_var_TL],
                            &p_TL,
                            thermo_properties_ptr);
                    
                    const Real T_BR = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getTemperature(
                            &rho[idx_cons_var_BR],
                            &p_BR,
                            thermo_properties_ptr);
                    
                    const Real T_TR = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getTemperature(
                            &rho[idx_cons_var_TR],
                            &p_TR,
                            thermo_properties_ptr);
                    
                    /*
                    const Real T_BL = epsilon_BL / (Real(717.5062667));
                    const Real T_BR = epsilon_BR / (Real(717.5062667));
                    const Real T_TR = epsilon_TR / (Real(717.5062667));
                    const Real T_TL = epsilon_TL / (Real(717.5062667));
                    */
                    // Bilinear interpolation to find image point temperature.
                    const Real T_f1 = (one - ip_ratio_0)*T_BL + ip_ratio_0*T_BR;
                    const Real T_f2 = (one - ip_ratio_0)*T_TL + ip_ratio_0*T_TR;
                    const Real T_ip = (one - ip_ratio_1)*T_f1 + ip_ratio_1*T_f2;
                    
                    // AFK bilinear interpolation to find image point density value.
                    const Real rho_f1 = (one - ip_ratio_0)*rho[idx_cons_var_BL] + ip_ratio_0*rho[idx_cons_var_BR];
                    const Real rho_f2 = (one - ip_ratio_0)*rho[idx_cons_var_TL] + ip_ratio_0*rho[idx_cons_var_TR];
                    const Real rho_ip = (one - ip_ratio_1)*rho_f1               + ip_ratio_1*rho_f2;

                    Real rho_gc = Real(0); 
                    Real T_gc = Real(0);     

                    if (d_bc_type_temperature == TEMPERATURE_IBC::ADIABATIC)
                    {
                        rho_gc = rho_ip;
                        T_gc = T_ip; 
                    }
                    else if (d_bc_type_temperature == TEMPERATURE_IBC::ISOTHERMAL)
                    {
                        T_gc = T_ip - ((d_ip + dist[idx_IB])/d_ip) * (T_ip - T_body);
                        rho_gc = (rho_ip * T_ip) / T_gc;
                    }
                   
                    // Calculating the total energy at the ghost cell value.
                    const Real epsilon_gc = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getInternalEnergyFromTemperature(
                            &rho_gc,
                            &T_gc,
                            thermo_properties_ptr);
                    
                    const Real E_gc = rho_gc*(epsilon_gc + half*(u_gc*u_gc + v_gc*v_gc));
                    
                    rho[idx_cons_var]   = rho_gc;
                    rho_u[idx_cons_var] = rho[idx_cons_var]*u_gc;
                    rho_v[idx_cons_var] = rho[idx_cons_var]*v_gc;
                    E[idx_cons_var]     = E_gc;
                }
                else if (mask[idx_IB] == int(IB_MASK::BODY))
                {
                    rho[idx_cons_var]   = rho_body;
                    rho_u[idx_cons_var] = rho_u_body;
                    rho_v[idx_cons_var] = rho_v_body;
                    E[idx_cons_var]     = E_body;
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const Real& rho_u_body = d_mom_body[0];
        const Real& rho_v_body = d_mom_body[1];
        const Real& rho_w_body = d_mom_body[2];
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int offset_0_cons_var = offset_cons_var[0];
        const int offset_1_cons_var = offset_cons_var[1];
        const int offset_2_cons_var = offset_cons_var[2];
        const int ghostcell_dim_0_cons_var = ghostcell_dims_cons_var[0];
        const int ghostcell_dim_1_cons_var = ghostcell_dims_cons_var[1];
        const int ghostcell_dim_2_cons_var = ghostcell_dims_cons_var[2];
        
        const int offset_0_IB = offset_IB[0];
        const int offset_1_IB = offset_IB[1];
        const int offset_2_IB = offset_IB[2];
        const int ghostcell_dim_0_IB = ghostcell_dims_IB[0];
        const int ghostcell_dim_1_IB = ghostcell_dims_IB[1];
        const int ghostcell_dim_2_IB = ghostcell_dims_IB[2];
        
        using Real = double;

        // Get the pointers to the data.
        Real* rho_u = data_momentum->getPointer(0);
        Real* rho_v = data_momentum->getPointer(1);
        Real* rho_w = data_momentum->getPointer(2);
        
        Real* norm_0 = data_surface_normal->getPointer(0);
        Real* norm_1 = data_surface_normal->getPointer(1);
        Real* norm_2 = data_surface_normal->getPointer(2);
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_cons_var = (i + offset_0_cons_var) +
                        (j + offset_1_cons_var)*ghostcell_dim_0_cons_var +
                        (k + offset_2_cons_var)*ghostcell_dim_0_cons_var*
                            ghostcell_dim_1_cons_var;
                    
                    const int idx_IB = (i + offset_0_IB) +
                        (j + offset_1_IB)*ghostcell_dim_0_IB +
                        (k + offset_2_IB)*ghostcell_dim_0_IB*
                            ghostcell_dim_1_IB;
                    
                    Real x[3];
                    x[0] = Real(patch_xlo[0]) + (Real(i) + half)*Real(dx[0]); // local x coordinates
                    x[1] = Real(patch_xlo[1]) + (Real(j) + half)*Real(dx[1]); // local y coordinates
                    x[2] = Real(patch_xlo[2]) + (Real(k) + half)*Real(dx[2]); // local z coordinates

                    if (mask[idx_IB] == int(IB_MASK::IB_GHOST))  
                    {
                        // using Ahmet's bilinear interpolation twice + my own linear interpolation
                        

                        // Start first with the back values
                        const Real x_B_ip = x[0] + (dist[idx_IB] + d_ip)*norm_0[idx_IB];
                        const Real y_B_ip = x[1] + (dist[idx_IB] + d_ip)*norm_1[idx_IB];
                        const Real z_B_ip = x[2] + (dist[idx_IB] + d_ip)*norm_2[idx_IB];

                        // Calculate indices of bottom left back in the bilinear interpolation stencil
                        const int ip_location_B_index_0 = int(floor((x_B_ip - patch_xlo[0] - half * Real(dx[0]))/Real(dx[0])));
                        const int ip_location_B_index_1 = int(floor((y_B_ip - patch_xlo[1] - half * Real(dx[1]))/Real(dx[1])));
                        const int ip_location_B_index_2 = int(floor((z_B_ip - patch_xlo[2] - half * Real(dx[2]))/Real(dx[2])));

                        // Calculate bottom back left coordinates
                        const Real x_ip_BBL = Real(patch_xlo[0]) + (Real(ip_location_B_index_0) + half)*Real(dx[0]);
                        const Real y_ip_BBL = Real(patch_xlo[1]) + (Real(ip_location_B_index_1) + half)*Real(dx[1]);
                        const Real z_ip_BBL = Real(patch_xlo[2]) + (Real(ip_location_B_index_2) + half)*Real(dx[2]);

                        // Calculate interpolation ratios
                        const Real ip_ratio_0 = (x_B_ip - x_ip_BBL)/Real(dx[0]);
                        const Real ip_ratio_1 = (y_B_ip - y_ip_BBL)/Real(dx[1]);
                        const Real ip_ratio_2 = (z_B_ip - z_ip_BBL)/Real(dx[2]);

                        // Calculate the indexes of the cells in interpolation for the back
                        const int idx_cons_var_BBL  = (ip_location_B_index_0 + offset_0_cons_var) +
                                                      (ip_location_B_index_1 + offset_1_cons_var) * ghostcell_dim_0_cons_var +
                                                      (ip_location_B_index_2 + offset_2_cons_var) * ghostcell_dim_1_cons_var * ghostcell_dim_0_cons_var;
        
                        const int idx_cons_var_BBR  = idx_cons_var_BBL + 1;
                        const int idx_cons_var_BTL  = idx_cons_var_BBL + ghostcell_dim_0_cons_var;
                        const int idx_cons_var_BTR  = idx_cons_var_BBL + ghostcell_dim_0_cons_var + 1;

                        // Now we do the front values
                        const Real x_F_ip = x[0] + (dist[idx_IB] + d_ip)*norm_0[idx_IB];
                        const Real y_F_ip = x[1] + (dist[idx_IB] + d_ip)*norm_1[idx_IB];
                        const Real z_F_ip = x[2] + (dist[idx_IB] + d_ip)*norm_2[idx_IB];

                        // Calculate indices of bottom left front in the bilinear interpolation stencil
                        const int ip_location_F_index_0 = int(floor((x_F_ip - patch_xlo[0] - half * Real(dx[0]))/Real(dx[0])));
                        const int ip_location_F_index_1 = int(floor((y_F_ip - patch_xlo[1] - half * Real(dx[1]))/Real(dx[1])));
                        const int ip_location_F_index_2 = int(floor((z_F_ip - patch_xlo[2] - half * Real(dx[2]))/Real(dx[2])));

                        // Calculate bottom left front coordinates
                        const Real x_ip_FBL = Real(patch_xlo[0]) + (Real(ip_location_F_index_0) + half)*Real(dx[0]);
                        const Real y_ip_FBL = Real(patch_xlo[1]) + (Real(ip_location_F_index_1) + half)*Real(dx[1]);
                        const Real z_ip_FBL = Real(patch_xlo[2]) + (Real(ip_location_F_index_2) + half)*Real(dx[2]);

                        // Calculate the indexes of the cells in interpolation for the front
                        const int idx_cons_var_FBL  = (ip_location_F_index_0 + offset_0_cons_var) +
                                                      (ip_location_F_index_1 + offset_1_cons_var) * ghostcell_dim_0_cons_var +
                                                      (ip_location_F_index_2 + offset_2_cons_var) * ghostcell_dim_1_cons_var * ghostcell_dim_0_cons_var;
        
                        const int idx_cons_var_FBR  = idx_cons_var_FBL + 1;
                        const int idx_cons_var_FTL  = idx_cons_var_FBL + ghostcell_dim_0_cons_var;
                        const int idx_cons_var_FTR  = idx_cons_var_FBL + ghostcell_dim_0_cons_var + 1;

                        // Now we do the final values. We should be in the z-y plane.
                        // Therefore I used the same as before but index 0 = index 2, index 1 = index 1.
                        const Real y_ip = x[1] + (dist[idx_IB] + d_ip)*norm_1[idx_IB];
                        const Real z_ip = x[2] + (dist[idx_IB] + d_ip)*norm_2[idx_IB];

                        // Calculate indices of final values in the bilinear interpolation stencil
                        const int ip_location_index_1 = int(floor((y_ip - patch_xlo[1] - half * Real(dx[1]))/Real(dx[1])));
                        const int ip_location_index_2 = int(floor((z_ip - patch_xlo[2] - half * Real(dx[2]))/Real(dx[2])));

                        // Calculate bottom back coordinates
                        const Real y_ip_BB = Real(patch_xlo[1]) + (Real(ip_location_index_1) + half)*Real(dx[1]);
                        const Real z_ip_BB = Real(patch_xlo[2]) + (Real(ip_location_index_2) + half)*Real(dx[2]);

                        // Calculate the indexes of the cells in interpolation for the front.

                        const int idx_cons_var_BB  =  (ip_location_index_2 + offset_2_cons_var) +
                                                      (ip_location_index_1 + offset_1_cons_var) * ghostcell_dim_2_cons_var;
                        // not sure if should be ghostcell_dim_2 or 1.
                        const int idx_cons_var_BF  = idx_cons_var_FBL + ghostcell_dim_2_cons_var;
                        const int idx_cons_var_TB  = idx_cons_var_FBL + ghostcell_dim_1_cons_var;
                        const int idx_cons_var_TF  = idx_cons_var_FBL + ghostcell_dim_1_cons_var + ghostcell_dim_2_cons_var;

                        // Bilinear interpolation to find image point x velocity value.
                        const Real u_IP_BBL = rho_u[idx_cons_var_BBL]/rho[idx_cons_var_BBL];
                        const Real u_IP_BBR = rho_u[idx_cons_var_BBR]/rho[idx_cons_var_BBR];
                        const Real u_IP_BTL = rho_u[idx_cons_var_BTL]/rho[idx_cons_var_BTL];
                        const Real u_IP_BTR = rho_u[idx_cons_var_BTR]/rho[idx_cons_var_BTR];
                        const Real u_IP_FBL = rho_u[idx_cons_var_FBL]/rho[idx_cons_var_FBL];
                        const Real u_IP_FBR = rho_u[idx_cons_var_FBR]/rho[idx_cons_var_FBR];
                        const Real u_IP_FTL = rho_u[idx_cons_var_FTL]/rho[idx_cons_var_FTL];
                        const Real u_IP_FTR = rho_u[idx_cons_var_FTR]/rho[idx_cons_var_FTR];
                        
                        // Bilinear interpolation to find image point y velocity value.
                        const Real v_IP_BBL = rho_v[idx_cons_var_BBL]/rho[idx_cons_var_BBL];
                        const Real v_IP_BBR = rho_v[idx_cons_var_BBR]/rho[idx_cons_var_BBR];
                        const Real v_IP_BTL = rho_v[idx_cons_var_BTL]/rho[idx_cons_var_BTL];
                        const Real v_IP_BTR = rho_v[idx_cons_var_BTR]/rho[idx_cons_var_BTR];
                        const Real v_IP_FBL = rho_v[idx_cons_var_FBL]/rho[idx_cons_var_FBL];
                        const Real v_IP_FBR = rho_v[idx_cons_var_FBR]/rho[idx_cons_var_FBR];
                        const Real v_IP_FTL = rho_v[idx_cons_var_FTL]/rho[idx_cons_var_FTL];
                        const Real v_IP_FTR = rho_v[idx_cons_var_FTR]/rho[idx_cons_var_FTR];

                        // Bilinear interpolation to find image point z velocity value.
                        const Real w_IP_BBL = rho_v[idx_cons_var_BBL]/rho[idx_cons_var_BBL];
                        const Real w_IP_BBR = rho_v[idx_cons_var_BBR]/rho[idx_cons_var_BBR];
                        const Real w_IP_BTL = rho_v[idx_cons_var_BTL]/rho[idx_cons_var_BTL];
                        const Real w_IP_BTR = rho_v[idx_cons_var_BTR]/rho[idx_cons_var_BTR];
                        const Real w_IP_FBL = rho_v[idx_cons_var_FBL]/rho[idx_cons_var_FBL];
                        const Real w_IP_FBR = rho_v[idx_cons_var_FBR]/rho[idx_cons_var_FBR];
                        const Real w_IP_FTL = rho_v[idx_cons_var_FTL]/rho[idx_cons_var_FTL];
                        const Real w_IP_FTR = rho_v[idx_cons_var_FTR]/rho[idx_cons_var_FTR];

                        // finding the new four interpolated velocity values - x-direction
                        //              B_F = BOTTOM_FRONT            T_B = TOP_BACK  
                        const Real u_IP_BB = (1 - ip_ratio_0) * u_IP_BBL + ip_ratio_0 * u_IP_BBR;
                        const Real u_IP_BF = (1 - ip_ratio_0) * u_IP_FBL + ip_ratio_0 * u_IP_FBR;
                        const Real u_IP_TB = (1 - ip_ratio_0) * u_IP_BTL + ip_ratio_0 * u_IP_BTR;
                        const Real u_IP_TF = (1 - ip_ratio_0) * u_IP_FTL + ip_ratio_0 * u_IP_FTR;
                        // finding the front/back interpolated velocity values
                        const Real u_IP_F   = (1 - ip_ratio_1) * u_IP_TF + ip_ratio_1 * u_IP_BF;
                        const Real u_IP_B   = (1 - ip_ratio_1) * u_IP_TB + ip_ratio_1 * u_IP_BB;
                        // finding the final value
                        const Real u_IP     = (1 - ip_ratio_2) * u_IP_F  + ip_ratio_2 * u_IP_B;

                        // finding the new four interpolated velocity values   - y-direction
                        //              B_F = BOTTOM_FRONT            T_B = TOP_BACK  
                        const Real v_IP_BB = (1 - ip_ratio_0) * v_IP_BBL + ip_ratio_0 * v_IP_BBR;
                        const Real v_IP_BF = (1 - ip_ratio_0) * v_IP_FBL + ip_ratio_0 * v_IP_FBR;
                        const Real v_IP_TB = (1 - ip_ratio_0) * v_IP_BTL + ip_ratio_0 * v_IP_BTR;
                        const Real v_IP_TF = (1 - ip_ratio_0) * v_IP_FTL + ip_ratio_0 * v_IP_FTR;
                        // finding the front/back interpolated velocity values
                        const Real v_IP_F   = (1 - ip_ratio_1) * v_IP_TF + ip_ratio_1 * v_IP_BF;
                        const Real v_IP_B   = (1 - ip_ratio_1) * v_IP_TB + ip_ratio_1 * v_IP_BB;
                        // finding the final velocity value
                        const Real v_IP     = (1 - ip_ratio_2) * v_IP_F  + ip_ratio_2 * v_IP_B;


                        // finding the new four interpolated velocity values - z-direction
                        //              B_F = BOTTOM_FRONT            T_B = TOP_BACK  
                        const Real w_IP_BB = (1 - ip_ratio_0) * w_IP_BBL + ip_ratio_0 * w_IP_BBR;
                        const Real w_IP_BF = (1 - ip_ratio_0) * w_IP_FBL + ip_ratio_0 * w_IP_FBR;
                        const Real w_IP_TB = (1 - ip_ratio_0) * w_IP_BTL + ip_ratio_0 * w_IP_BTR;
                        const Real w_IP_TF = (1 - ip_ratio_0) * w_IP_FTL + ip_ratio_0 * w_IP_FTR;
                        // finding the front/back interpolated velocity values
                        const Real w_IP_F   = (1 - ip_ratio_1) * w_IP_TF + ip_ratio_1 * w_IP_BF;
                        const Real w_IP_B   = (1 - ip_ratio_1) * w_IP_TB + ip_ratio_1 * w_IP_BB;
                        // finding the final value
                        const Real w_IP     = (1 - ip_ratio_2) * w_IP_F  + ip_ratio_2 * w_IP_B;

                        // Compute velocities at the ghost cell for Euler equation (no-slip for normal direction and slip for tangential direction)
                        Real u_gc = Real(0);  // u velocity of the ghost cell
                        Real v_gc = Real(0);  // v velocity of the ghost cell
                        Real w_gc = Real(0);  // w velocity of the ghost cell
                        Real vel_ip_n = Real(0);
                        Real vel_ip_t = Real(0);
                        Real vel_gc_n = Real(0);
                        Real vel_gc_t = Real(0);

                        //NEEDS TO BE UPDATED
                        if(d_bc_type_velocity == VELOCITY_IBC::SLIP)
                        {
                            vel_ip_n =  u_IP*norm_0[idx_IB] + v_ip*norm_1[idx_IB];
                            vel_ip_t = -u_IP*norm_1[idx_IB] + v_ip*norm_0[idx_IB];
                            
                            
                            vel_gc_n = vel_ip_n - ((d_ip + dist[idx_IB]) / d_ip) * (vel_ip_n);
                            vel_gc_t = vel_ip_t;
                            
                            u_gc = vel_gc_n*norm_0[idx_IB] - vel_gc_t*norm_1[idx_IB];
                            v_gc = vel_gc_n*norm_1[idx_IB] + vel_gc_t*norm_0[idx_IB]; 
                        }
                        else if(d_bc_type_velocity == VELOCITY_IBC::NO_SLIP) 
                        {
                            u_gc = u_IP - ((d_ip + dist[idx_IB])/d_ip) * (u_IP);
                            v_gc = v_IP - ((d_ip + dist[idx_IB])/d_ip) * (v_IP);
                            w_gc = w_IP - ((d_ip + dist[idx_IB])/d_ip) * (w_IP);
                        }

                        const Real epsilon_BB = E[idx_cons_var_BB]/rho[idx_cons_var_BB] - half*(u_IP_BB*u_IP_BB + v_IP_BL*v_IP_BB + w_IP_BL*w_IP_BB);
                        const Real epsilon_BF = E[idx_cons_var_TF]/rho[idx_cons_var_TF] - half*(u_IP_BF*u_IP_BF + v_IP_BF*v_IP_BF + w_IP_BF*w_IP_BF);
                        const Real epsilon_TB = E[idx_cons_var_TB]/rho[idx_cons_var_TB] - half*(u_IP_TB*u_IP_TB + v_IP_TB*v_IP_TB + w_IP_TB*w_IP_TB);
                        const Real epsilon_TF = E[idx_cons_var_TF]/rho[idx_cons_var_TF] - half*(u_IP_TR*u_IP_TF + v_IP_TR*v_IP_TF + w_IP_TF*w_IP_TF);

                        // Compute the pressure values in the stencils.
                        const Real p_BL = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho[idx_cons_var_BB],
                                &epsilon_BL,
                                thermo_properties_ptr);
                        
                        const Real p_TL = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho[idx_cons_var_BF],
                                &epsilon_TL,
                                thermo_properties_ptr);
                        
                        const Real p_BR = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho[idx_cons_var_TB],
                                &epsilon_BR,
                                thermo_properties_ptr);
                        
                        const Real p_TR = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho[idx_cons_var_TF],
                                &epsilon_TR,
                                thermo_properties_ptr);
                        
                        // Compute the temperature values in the stencils.
                        const Real T_BL = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getTemperature(
                                &rho[idx_cons_var_BB],
                                &p_BL,
                                thermo_properties_ptr);
                        
                        const Real T_TL = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getTemperature(
                                &rho[idx_cons_var_BF],
                                &p_TL,
                                thermo_properties_ptr);
                        
                        const Real T_BR = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getTemperature(
                                &rho[idx_cons_var_TB],
                                &p_BR,
                                thermo_properties_ptr);
                        
                        const Real T_TR = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getTemperature(
                                &rho[idx_cons_var_TF],
                                &p_TR,
                                thermo_properties_ptr);

                        /*
                        const Real T_BL = epsilon_BL / (Real(717.5062667));
                        const Real T_BR = epsilon_BR / (Real(717.5062667));
                        const Real T_TR = epsilon_TR / (Real(717.5062667));
                        const Real T_TL = epsilon_TL / (Real(717.5062667));
                        */

                        // Bilinear interpolation to find image point temperature.
                        const Real T_B     = (1 - ip_ratio_1)*T_TB + ip_ratio_1*T_BB;
                        const Real T_F     = (1 - ip_ratio_1)*T_TF + ip_ratio_1*T_BF;
                        const Real T_ip = (1 - ip_ratio_2)*T_B  + ip_ratio_2*T_F;
                        
                        // AFK bilinear interpolation to find image point density value.
                        const Real rho_B = (1 - ip_ratio_1)*rho[idx_cons_var_TB] + ip_ratio_1*rho[idx_cons_var_BB];
                        const Real rho_F = (1 - ip_ratio_1)*rho[idx_cons_var_TF] + ip_ratio_1*rho[idx_cons_var_BF];
                        const Real rho_ip = (1 - ip_ratio_2)*rho_B               + ip_ratio_2*rho_F;

                        Real rho_gc = Real(0); 
                        Real T_gc = Real(0);     

                        if (d_bc_type_temperature == TEMPERATURE_IBC::ADIABATIC)
                        {
                            rho_gc = rho_ip;
                            T_gc = T_ip; 
                        }
                        else if (d_bc_type_temperature == TEMPERATURE_IBC::ISOTHERMAL)
                        {
                            T_gc = T_ip - ((d_ip + dist[idx_IB])/d_ip) * (T_ip - T_body);
                            rho_gc = (rho_ip * T_ip) / T_gc;
                        }
                    
                        // Calculating the total energy at the ghost cell value.
                        const Real epsilon_gc = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getInternalEnergyFromTemperature(
                                &rho_gc,
                                &T_gc,
                                thermo_properties_ptr);
                        
                        const Real E_gc = rho_gc*(epsilon_gc + half*(u_gc*u_gc + v_gc*v_gc + w_gc*w_gc));
                    
                        rho[idx_cons_var]   = rho_body;
                        rho_u[idx_cons_var] = rho_u_body;
                        rho_v[idx_cons_var] = rho_v_body;
                        rho_w[idx_cons_var] = rho_w_body;
                        E[idx_cons_var]     = E_body;
                    }
                    else if (mask[idx_IB] == int(IB_MASK::BODY))
                    {
                        rho[idx_cons_var]   = rho_body;
                        rho_u[idx_cons_var] = rho_u_body;
                        rho_v[idx_cons_var] = rho_v_body;
                        rho_w[idx_cons_var] = rho_w_body;
                        E[idx_cons_var]     = E_body;
                    }
                }
            }
        }
    }
}
