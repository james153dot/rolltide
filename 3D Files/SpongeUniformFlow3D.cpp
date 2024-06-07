#include "flow/flow_models/FlowModelSpecialSourceTerms.hpp"

/*
 * Add the effects of the special source terms.
 */

/* 
Edges and verticies are standardized and taken from 
https://s3-eu-west-1.amazonaws.com/ppreviews-plos-725668748/2192469/preview.jpg
https://qph.cf2.quoracdn.net/main-qimg-80e9bd98f686c2c90015bac9d6326842
*/


void
FlowModelSpecialSourceTerms::computeSpecialSourceTermsOnPatch(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& source,
    const hier::Patch& patch,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_variables,
    const std::unordered_map<std::string, Real>& monitoring_statistics_map,
    const double time,
    const double dt,
    const int RK_step_number)
{
    if ((d_project_name != "3D uniform flow") ) 
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '3D uniform flow'!\n"
            << "'project_name' = '"
            << d_project_name
            << "' is given."
            << std::endl);
    }
    
    if (d_dim != tbox::Dimension(3))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Dimension of problem should be 3!"
            << std::endl);
    }
    
    if (d_special_source_exterior == false)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The 'special_source_exterior' option should be true!"
            << std::endl);
    }
    
    TBOX_ASSERT(d_source_terms_db->keyExists("sponge_rate"));
    
    Real sponge_rate = Real(0);
    if (d_source_terms_db->keyExists("sponge_rate"))
    {
        sponge_rate = d_source_terms_db->getReal("sponge_rate");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'sponge_rate' found in data for source terms."
            << std::endl);
    }
    
    TBOX_ASSERT(d_source_terms_db->keyExists("rho_inf"));
    
    Real rho_inf = Real(0);
    if (d_source_terms_db->keyExists("rho_inf"))
    {
        rho_inf = d_source_terms_db->getReal("rho_inf");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'rho_inf' found in data for source terms."
            << std::endl);
    }
    
    TBOX_ASSERT(d_source_terms_db->keyExists("u_inf"));
    
    Real u_inf = Real(0);
    if (d_source_terms_db->keyExists("u_inf"))
    {
        u_inf = d_source_terms_db->getReal("u_inf");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'u_inf' found in data for source terms."
            << std::endl);
    }
    
    TBOX_ASSERT(d_source_terms_db->keyExists("v_inf"));
    
    Real v_inf = Real(0);
    if (d_source_terms_db->keyExists("v_inf"))
    {
        v_inf = d_source_terms_db->getReal("v_inf");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'v_inf' found in data for source terms."
            << std::endl);
    }

    TBOX_ASSERT(d_source_terms_db->keyExists("w_inf")); //newly added
    
    Real w_inf = Real(0);
    if (d_source_terms_db->keyExists("w_inf"))
    {
        w_inf = d_source_terms_db->getReal("w_inf");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'w_inf' found in data for source terms."
            << std::endl);
    }
    
    TBOX_ASSERT(d_source_terms_db->keyExists("p_inf"));
    
    Real p_inf = Real(0);
    if (d_source_terms_db->keyExists("p_inf"))
    {
        p_inf = d_source_terms_db->getReal("p_inf");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'p_inf' found in data for source terms."
            << std::endl);
    }

    TBOX_ASSERT(d_source_terms_db->keyExists("D"));
    
    Real D = Real(1);
    
    if (d_source_terms_db->keyExists("D"))
    {
        D = d_source_terms_db->getReal("D");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'D' found in data for source terms."
            << std::endl);
    }
    
    Real gamma = Real(7)/Real(5);
    
    if (d_source_terms_db->keyExists("gamma"))
    {
        gamma = d_source_terms_db->getReal("gamma");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'gamma' found in data for source terms."
            << std::endl);
    }

    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    std::vector<Real*> S;
    S.reserve(d_num_eqn);
    for (int si = 0; si < d_num_eqn; si++)
    {
        S.push_back(source->getPointer(si));
    }
    
    /*
     * Get the numbers of ghost cells source and conservative variables.
     */
    const hier::IntVector num_ghosts_source     = source->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_source = source->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_cons_var     = conservative_variables[0]->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_cons_var = conservative_variables[0]->getGhostBox().numberCells();
    
    const double* const dx = patch_geom->getDx();
    const double* const patch_xlo = patch_geom->getXLower();
    
    // Get the dimensions of box that covers the interior of Patch.
    hier::Box patch_box = patch.getBox();
    const hier::IntVector patch_dims = patch_box.numberCells();
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > density      = conservative_variables[0];
    HAMERS_SHARED_PTR<pdat::CellData<Real> > momentum     = conservative_variables[1];
    HAMERS_SHARED_PTR<pdat::CellData<Real> > total_energy = conservative_variables[2];
    
    Real* rho   = density->getPointer(0);
    Real* rho_u = momentum->getPointer(0);
    Real* rho_v = momentum->getPointer(1);
    Real* rho_w = momentum->getPointer(2);
    Real* E     = total_energy->getPointer(0);
    
    TBOX_ASSERT(d_source_terms_db != nullptr);
    
    const double* const domain_xlo = d_grid_geometry->getXLower();
    const double* const domain_xhi = d_grid_geometry->getXUpper();
    
    if (d_project_name == "3D uniform flow")
    {
        for (int k = 0; k < patch_dims[2]; k++)
        {
            for (int j = 0; j < patch_dims[1]; j++)
            {
                for (int i = 0; i < patch_dims[0]; i++)
                {
                    // Compute the linear indices.
                    const int idx_source =         (i + num_ghosts_source[0]) +
                        (j + num_ghosts_source[1]) * ghostcell_dims_source[0] +
                        (k + num_ghosts_source[2]) * ghostcell_dims_source[1] * ghostcell_dims_source[0];
                    
                    const int idx_cons_var = (i + num_ghosts_cons_var[0]) +
                        (j + num_ghosts_cons_var[1]) * ghostcell_dims_cons_var[0] +
                        (k + num_ghosts_cons_var[2]) * ghostcell_dims_cons_var[1] * ghostcell_dims_cons_var[0];
                    
                    // Compute the coordinates.
                    Real x[3];
                    x[0] = patch_xlo[0] + (Real(i) + Real(1)/Real(2))*Real(dx[0]);
                    x[1] = patch_xlo[1] + (Real(j) + Real(1)/Real(2))*Real(dx[1]);
                    x[2] = patch_xlo[1] + (Real(k) + Real(1)/Real(2))*Real(dx[2]);
                    
                    const Real half  = Real(1)/Real(2);
                    const Real two_third = Real(2)/Real(3);

                    //Left face sponge region
                    if (x[0] <= d_special_source_box_lo[0])
                    {
                        const Real u_ref = u_inf;
                        const Real v_ref = v_inf;
                        const Real w_ref = w_inf;
                        
                        const Real rho_ref = rho_inf;
                        const Real p_ref   = p_inf;
                        
                        const Real rho_u_ref = rho_ref * u_ref;
                        const Real rho_v_ref = rho_ref * v_ref;
                        const Real rho_w_ref = rho_ref * w_ref;
                        const Real E_ref     = p_ref/(gamma - Real(1)) + Real(1)/Real(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);
                        
                        const Real rho_p   = rho[idx_cons_var]   - rho_ref;
                        const Real rho_u_p = rho_u[idx_cons_var] - rho_u_ref;
                        const Real rho_v_p = rho_v[idx_cons_var] - rho_v_ref;
                        const Real rho_w_p = rho_w[idx_cons_var] - rho_w_ref;
                        const Real E_p     = E[idx_cons_var]     - E_ref;

                        // Left sponge calculations
                        /*
                        Real xi_b = std::pow((Real(1) - (x[0] - domain_xlo[0])/(d_special_source_box_lo[0] - domain_xlo[0])), Real(3));
                        //Real xi_b = (Real(1) - (x[0] - domain_xlo[0])/(d_special_source_box_lo[0] - domain_xlo[0]));
                        xi_b      *= sponge_rate; // mask value needs to be improved
                        */
                        // 0927 AFK Sponge Update

                        const Real erf_start_lo  = half*(domain_xlo[0] - d_special_source_box_lo[0]) + d_special_source_box_lo[0]; //center of erf is half of the way into sponge
                        const Real erf_offset_lo = Real(-0.5) * erf((d_special_source_box_lo[0]-erf_start_lo)/(D*Real(0.5))) + Real(0.5); //value of erf at start of sponge 
                        Real xi_b                = Real(-0.5) * erf((x[0]-erf_start_lo)/(D*Real(0.5))) + Real(0.5) - erf_offset_lo; //subtract value of erf at start of sponge to start xi_b at zero
                        xi_b                     *= sponge_rate; 

                        // Edges 4, 8, 9, 12.
                        if x[1] <= d_special_source_box_lo[1]   || x[1] >= d_special_source_box_hi[1]  ||  x[2] <= d_special_source_box_lo[2] || x[2] >= d_special_source_box_hi[2]
                        {
                            x_ib = x_ib * half;
                        }
                        }
                        // Corners A, D, E, H.
                        if ((x[1] <= d_special_source_box_lo[1] || x[1] >= d_special_source_box_hi[1]) && (x[2] <= d_special_source_box_lo[2] || x[2] >= d_special_source_box_hi[2]))
                        {
                             x_ib = x_ib * two_third; 
                        }
                        S[0][idx_source] -= dt*xi_b*rho_p;
                        S[1][idx_source] -= dt*xi_b*rho_u_p;
                        S[2][idx_source] -= dt*xi_b*rho_v_p;
                        S[3][idx_source] -= dt*xi_b*rho_w_p;
                        S[4][idx_source] -= dt*xi_b*E_p;
                    }
                    
                    // Right face sponge region
                    if (x[0] >= d_special_source_box_hi[0])
                    {                    
                        const Real u_ref = u_inf;
                        const Real v_ref = v_inf;
                        const Real w_ref = w_inf;
                        
                        const Real rho_ref = rho_inf;
                        const Real p_ref   = p_inf;

                        const Real rho_u_ref = rho_ref * u_ref;
                        const Real rho_v_ref = rho_ref * v_ref;
                        const Real rho_w_ref = rho_ref * w_ref;
                        const Real E_ref     = p_ref/(gamma - Real(1)) + Real(1)/Real(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);

                        const Real rho_p   = rho[idx_cons_var]   - rho_ref;
                        const Real rho_u_p = rho_u[idx_cons_var] - rho_u_ref;
                        const Real rho_v_p = rho_v[idx_cons_var] - rho_v_ref;
                        const Real rho_w_p = rho_w[idx_cons_var] - rho_w_ref;
                        const Real E_p     = E[idx_cons_var]     - E_ref;
                        /*
                        // Right sponge calculations
                        Real xi_b      = std::pow((x[0]-d_special_source_box_hi[0])/(domain_xhi[0]-d_special_source_box_hi[0]), Real(3)); // mask value needs to be improved 
                        //Real xi_b      = (x[0]-d_special_source_box_hi[0])/(domain_xhi[0]-d_special_source_box_hi[0]);
                        xi_b          *= sponge_rate;
                        */
                        // 0927 AFK Sponge Update
                        const Real erf_start_hi  = half*(domain_xhi[0]-d_special_source_box_hi[0]) + d_special_source_box_hi[0];         //center of erf is half of the way into sponge
                        const Real erf_offset_hi = Real(0.5) * erf((d_special_source_box_hi[0]-erf_start_hi)/(D*Real(0.5))) + Real(0.5); //value of erf at start of sponge
                        Real xi_b                = Real(0.5) * erf((x[0]-erf_start_hi)/(D*Real(0.5))) + Real(0.5) - erf_offset_hi;       //subtract value of erf at start of sponge to start xi_b at zero
                        xi_b                     *= sponge_rate; 
                        
                        // Edges 2, 6, 10, 11.
                        if x[1]   <= d_special_source_box_lo[1] || x[1] >= d_special_source_box_hi[1]  ||  x[2] <= d_special_source_box_lo[2] || x[2] >= d_special_source_box_hi[2]
                        {
                            x_ib = x_ib * half;
                        }
                        // Corners B, C, F, G.
                        if ((x[1] <= d_special_source_box_lo[1] || x[1] >= d_special_source_box_hi[1]) && (x[2] <= d_special_source_box_lo[2] || x[2] >= d_special_source_box_hi[2]))
                        {
                             x_ib = x_ib * two_third; 
                        }
                        S[0][idx_source] -= dt*xi_b*rho_p;
                        S[1][idx_source] -= dt*xi_b*rho_u_p;
                        S[2][idx_source] -= dt*xi_b*rho_v_p;
                        S[3][idx_source] -= dt*xi_b*rho_w_p;
                        S[4][idx_source] -= dt*xi_b*E_p;
                    }
                    
                    // Bottom face sponge region
                    if (x[1] <= d_special_source_box_lo[1])
                    {
                        const Real u_ref = u_inf;
                        const Real v_ref = v_inf;
                        const Real w_ref = w_inf;
                        
                        const Real rho_ref = rho_inf;
                        const Real p_ref   = p_inf;

                        const Real rho_u_ref = rho_ref * u_ref;
                        const Real rho_v_ref = rho_ref * v_ref;
                        const Real rho_w_ref = rho_ref * w_ref;
                        const Real E_ref     = p_ref/(gamma - Real(1)) + Real(1)/Real(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);

                        const Real rho_p   = rho[idx_cons_var]   - rho_ref;
                        const Real rho_u_p = rho_u[idx_cons_var] - rho_u_ref;
                        const Real rho_v_p = rho_v[idx_cons_var] - rho_v_ref;
                        const Real rho_w_p = rho_w[idx_cons_var] - rho_w_ref;
                        const Real E_p     = E[idx_cons_var]     - E_ref;

                        /*
                        // Bottom sponge calculations
                        Real yi_b = std::pow((Real(1) - (x[1] - domain_xlo[1])/(d_special_source_box_lo[1] - domain_xlo[1])), Real(3));
                        //Real yi_b = (Real(1) - (x[1] - domain_xlo[1])/(d_special_source_box_lo[1] - domain_xlo[1]));
                        yi_b     *= sponge_rate; // mask value needs to be improved 
                        */
                        // 0927 AFK Sponge Update
                        const Real erf_start_lo  = half*(domain_xlo[1] - d_special_source_box_lo[1]) + d_special_source_box_lo[1]; //center of erf is 3/4 of the way into sponge
                        const Real erf_offset_lo = Real(-0.5) * erf((d_special_source_box_lo[1]-erf_start_lo)/(D*Real(0.5))) + Real(0.5); //value of erf at start of sponge 
                        Real yi_b                = Real(-0.5) * erf((x[1]-erf_start_lo)/(D*Real(0.5))) + Real(0.5) - erf_offset_lo; //subtract value of erf at start of sponge to start xi_b at zero
                        yi_b                     *= sponge_rate; 

                        // Edges 1, 2, 3, 4.
                        if x[0]   <= d_special_source_box_lo[0] || x[0] >= d_special_source_box_hi[0]  || x[2] <= d_special_source_box_lo[2]  || x[2] >= d_special_source_box_hi[2]

                        {
                            y_ib = y_ib * half;
                        }
                        // Corners A, B, C, D.
                        if ((x[0] <= d_special_source_box_lo[0] || x[0] >= d_special_source_box_hi[0]) && (x[2] <= d_special_source_box_lo[2] || x[2] >= d_special_source_box_hi[2]))
                        {
                            y_ib = y_ib * two_third;
                        }
                        S[0][idx_source] -= dt*yi_b*rho_p;
                        S[1][idx_source] -= dt*yi_b*rho_u_p;
                        S[2][idx_source] -= dt*yi_b*rho_v_p;
                        S[3][idx_source] -= dt*yi_b*rho_w_p;
                        S[4][idx_source] -= dt*yi_b*E_p;
                    }
                    // Top sponge region
                    if (x[1] >= d_special_source_box_hi[1])
                    {                    
                        const Real u_ref = u_inf;
                        const Real v_ref = v_inf;
                        const Real w_ref = w_inf;
                        
                        const Real rho_ref = rho_inf;
                        const Real p_ref   = p_inf;

                        const Real rho_u_ref = rho_ref * u_ref;
                        const Real rho_v_ref = rho_ref * v_ref;
                        const Real rho_w_ref = rho_ref * w_ref;
                        const Real E_ref     = p_ref/(gamma - Real(1)) + Real(1)/Real(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);

                        const Real rho_p   = rho[idx_cons_var]   - rho_ref;
                        const Real rho_u_p = rho_u[idx_cons_var] - rho_u_ref;
                        const Real rho_v_p = rho_v[idx_cons_var] - rho_v_ref;
                        const Real rho_w_p = rho_w[idx_cons_var] - rho_w_ref;
                        const Real E_p     = E[idx_cons_var]     - E_ref;
                        /*
                        //Top sponge region calculations
                        Real yi_b      = std::pow((x[1]-d_special_source_box_hi[1])/(domain_xhi[1]-d_special_source_box_hi[1]), Real(3)); // mask value needs to be improved 
                        //Real yi_b      = (x[1]-d_special_source_box_hi[1])/(domain_xhi[1]-d_special_source_box_hi[1]);
                        yi_b          *= sponge_rate;
                        */
                        // 0927 AFK Sponge Update
                        const Real erf_start_hi  = half*(domain_xhi[2]-d_special_source_box_hi[2]) + d_special_source_box_hi[2]; //center of erf is half of the way into sponge
                        const Real erf_offset_hi = Real(0.5) * erf((d_special_source_box_hi[2]-erf_start_hi)/(D*Real(0.5))) + Real(0.5); //value of erf at start of sponge
                        Real zi_b                = Real(0.5) * erf((x[2]-erf_start_hi)/(D*Real(0.5))) + Real(0.5) - erf_offset_hi; //subtract value of erf at start of sponge to start xi_b at zero
                        zi_b                     *= sponge_rate; 
                        
                        // Edges 5, 6, 7, 8.
                        if x[0] <= d_special_source_box_lo[0] || x[0] >= d_special_source_box_hi[0] || x[2] <= d_special_source_box_lo[2] || x[2] >= d_special_source_box_hi[2]

                        {
                            y_ib = y_ib * half;
                        }
                        // Corners E, F, G, H.
                        if ((x[0] <= d_special_source_box_lo[0] || x[0] >= d_special_source_box_hi[0]) && (x[2] <= d_special_source_box_lo[2] || x[2] >= d_special_source_box_hi[2]))
                        {
                            y_ib = y_ib * two_third;
                        }
                        S[0][idx_source] -= dt*yi_b*rho_p;
                        S[1][idx_source] -= dt*yi_b*rho_u_p;
                        S[2][idx_source] -= dt*yi_b*rho_v_p;
                        S[3][idx_source] -= dt*yi_b*rho_w_p;
                        S[4][idx_source] -= dt*yi_b*E_p;
                    }
                    // Front sponge region
                    if (x[2] >= d_special_source_box_hi[2])
                    {                    
                        const Real u_ref = u_inf;
                        const Real v_ref = v_inf;
                        const Real w_ref = w_inf;
                        
                        const Real rho_ref = rho_inf;
                        const Real p_ref   = p_inf;

                        const Real rho_u_ref = rho_ref * u_ref;
                        const Real rho_v_ref = rho_ref * v_ref;
                        const Real rho_w_ref = rho_ref * w_ref;
                        const Real E_ref     = p_ref/(gamma - Real(1)) + Real(1)/Real(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);

                        const Real rho_p   = rho[idx_cons_var]   - rho_ref;
                        const Real rho_u_p = rho_u[idx_cons_var] - rho_u_ref;
                        const Real rho_v_p = rho_v[idx_cons_var] - rho_v_ref;
                        const Real rho_w_p = rho_w[idx_cons_var] - rho_w_ref;
                        const Real E_p     = E[idx_cons_var]     - E_ref;

                        /*
                        //Top sponge region calculations
                        Real yi_b      = std::pow((x[1]-d_special_source_box_hi[1])/(domain_xhi[1]-d_special_source_box_hi[1]), Real(3)); // mask value needs to be improved 
                        //Real yi_b      = (x[1]-d_special_source_box_hi[1])/(domain_xhi[1]-d_special_source_box_hi[1]);
                        yi_b          *= sponge_rate;
                        */
                        // 0927 AFK Sponge Update
                        const Real erf_start_hi  = half*(domain_xhi[2]-d_special_source_box_hi[2]) + d_special_source_box_hi[2]; //center of erf is half of the way into sponge
                        const Real erf_offset_hi = Real(0.5) * erf((d_special_source_box_hi[2]-erf_start_hi)/(D*Real(0.5))) + Real(0.5); //value of erf at start of sponge
                        Real zi_b                = Real(0.5) * erf((x[2]-erf_start_hi)/(D*Real(0.5))) + Real(0.5) - erf_offset_hi; //subtract value of erf at start of sponge to start xi_b at zero
                        zi_b                     *= sponge_rate; 
                        
                        // Edges 3, 7, 11, 12.
                        if x[0] <= d_special_source_box_lo[0]   || x[0] >= d_special_source_box_hi[0]  || x[1] <= d_special_source_box_lo[1]  || x[1] >= d_special_source_box_hi[1]
                        {
                            z_ib = z_ib * half;
                        }
                        // Corners A, B, E, F.
                        if ((x[0] <= d_special_source_box_lo[0] || x[0] >= d_special_source_box_hi[0]) && (x[1] <= d_special_source_box_lo[1] || x[1] >= d_special_source_box_hi[1]))
                        {
                            z_ib = z_ib * two_third;
                        }
                        S[0][idx_source] -= dt*zi_b*rho_p;
                        S[1][idx_source] -= dt*zi_b*rho_u_p;
                        S[2][idx_source] -= dt*zi_b*rho_v_p;
                        S[3][idx_source] -= dt*zi_b*rho_w_p;
                        S[4][idx_source] -= dt*zi_b*E_p;
                    }
                    // Back Sponge Region
                    if (x[2] <= d_special_source_box_lo[2])
                    {                    
                        const Real u_ref = u_inf;
                        const Real v_ref = v_inf;
                        const Real w_ref = w_inf;
                        
                        const Real rho_ref = rho_inf;
                        const Real p_ref   = p_inf;

                        const Real rho_u_ref = rho_ref * u_ref;
                        const Real rho_v_ref = rho_ref * v_ref;
                        const Real rho_w_ref = rho_ref * w_ref;
                        const Real E_ref     = p_ref/(gamma - Real(1)) + Real(1)/Real(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);

                        const Real rho_p   = rho[idx_cons_var]   - rho_ref;
                        const Real rho_u_p = rho_u[idx_cons_var] - rho_u_ref;
                        const Real rho_v_p = rho_v[idx_cons_var] - rho_v_ref;
                        const Real rho_w_p = rho_w[idx_cons_var] - rho_w_ref;
                        const Real E_p     = E[idx_cons_var]     - E_ref;

                        /*
                        //Top sponge region calculations
                        Real yi_b      = std::pow((x[1]-d_special_source_box_hi[1])/(domain_xhi[1]-d_special_source_box_hi[1]), Real(3)); // mask value needs to be improved 
                        //Real yi_b      = (x[1]-d_special_source_box_hi[1])/(domain_xhi[1]-d_special_source_box_hi[1]);
                        yi_b          *= sponge_rate;
                        */
                        // 0927 AFK Sponge Update
                        const Real erf_start_lo  = half*(domain_xlo[2]-d_special_source_box_lo[2]) + d_special_source_box_lo[2]; //center of erf is half of the way into sponge
                        const Real erf_offset_lo = Real(-0.5) * erf((d_special_source_box_lo[2]-erf_start_lo)/(D*Real(0.5))) + Real(0.5); //value of erf at start of sponge
                        Real zi_b                = Real(-0.5) * erf((x[2]-erf_start_lo)/(D*Real(0.5))) + Real(0.5) - erf_offset_lo; //subtract value of erf at start of sponge to start xi_b at zero
                        zi_b                     *= sponge_rate; 
                        
                        // Edges 1, 5, 9, 10.
                        if x[0] <= d_special_source_box_lo[0]   || x[0] >= d_special_source_box_hi[0]  || x[1] <= d_special_source_box_lo[1]  || x[1] >= d_special_source_box_hi[1]
                        {
                            z_ib = z_ib * half;
                        }
                        // Corners C, D, G, H.
                        if ((x[0] <= d_special_source_box_lo[0] || x[0] >= d_special_source_box_hi[0]) && (x[1] <= d_special_source_box_lo[1] || x[1] >= d_special_source_box_hi[1]))
                        {
                            z_ib = z_ib * two_third;
                        }
                        S[0][idx_source] -= dt*zi_b*rho_p;
                        S[1][idx_source] -= dt*zi_b*rho_u_p;
                        S[2][idx_source] -= dt*zi_b*rho_v_p;
                        S[3][idx_source] -= dt*zi_b*rho_w_p;
                        S[4][idx_source] -= dt*zi_b*E_p;
                    }
                }
            }
        }
    }
}


void
FlowModelSpecialSourceTerms::putToRestart(const HAMERS_SHARED_PTR<tbox::Database>& restart_source_terms_db)
{
    putToRestartBase(restart_source_terms_db);
    
    Real sponge_rate = Real(0);
    if (d_source_terms_db->keyExists("sponge_rate"))
    {
        sponge_rate = d_source_terms_db->getReal("sponge_rate");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'sponge_rate' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putReal("sponge_rate", sponge_rate);

    Real rho_inf = Real(0);
    if (d_source_terms_db->keyExists("rho_inf"))
    {
        rho_inf = d_source_terms_db->getReal("rho_inf");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'rho_inf' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putReal("rho_inf", rho_inf);
    
    Real u_inf = Real(0);
    if (d_source_terms_db->keyExists("u_inf"))
    {
        u_inf = d_source_terms_db->getReal("u_inf");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'u_inf' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putReal("u_inf", u_inf);

    Real v_inf = Real(0);
    if (d_source_terms_db->keyExists("v_inf"))
    {
        v_inf = d_source_terms_db->getReal("v_inf");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'v_inf' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putReal("v_inf", v_inf);

    Real w_inf = Real(0);
    if (d_source_terms_db->keyExists("w_inf"))
    {
        w_inf = d_source_terms_db->getReal("w_inf");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'w_inf' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putReal("w_inf", w_inf);

    Real p_inf = Real(0);
    if (d_source_terms_db->keyExists("p_inf"))
    {
        p_inf = d_source_terms_db->getReal("p_inf");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'p_inf' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putReal("p_inf", p_inf);
    
    Real D = Real(1);
    if (d_source_terms_db->keyExists("D"))
    {
        D = d_source_terms_db->getReal("D");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'D' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putReal("D", D);

    Real gamma = Real(7)/Real(5);
    if (d_source_terms_db->keyExists("gamma"))
    {
        gamma = d_source_terms_db->getReal("gamma");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'gamma' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putReal("gamma", gamma);

}
