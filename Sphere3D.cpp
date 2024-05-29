#include "util/immersed_boundaries/ImmersedBoundaries.hpp"
/* writes down all of util/immersed_boundaries/whatever */
void
ImmersedBoundaries::setImmersedBoundaryVariablesOnPatch(
    const hier::Patch& patch,
/* every single patch for which the variables are used */
    const double data_time,
    const bool initial_time,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const HAMERS_SHARED_PTR<pdat::CellData<int> >& data_mask,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_wall_distance,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_surface_normal)
{
    NULL_USE(data_time);
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
/* sets patch geometry to said cartesian geometry*/
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);

/* patch gemetry = not null*/
#endif
    
    const double* const dx = patch_geom->getDx(); // dx = cell size 
    const double* const patch_xlo = patch_geom->getXLower();
    
    const hier::IntVector num_ghosts = data_mask->getGhostCellWidth();
    const hier::IntVector ghostcell_dims = data_mask->getGhostBox().numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(num_ghosts == data_wall_distance->getGhostCellWidth());
    TBOX_ASSERT(num_ghosts == data_surface_normal->getGhostCellWidth());
#endif
    
    /*
     * Get the pointers to the data.
     */
    int* mask      = data_mask->getPointer(0);
    double* dist   = data_wall_distance->getPointer(0);
    double* norm_0 = data_surface_normal->getPointer(0);
    double* norm_1 = data_surface_normal->getPointer(1);
    double* norm_2 = data_surface_normal->getPointer(2);
    
    /*
     * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
     */
    
    const int domain_lo_0 = domain_lo[0];
    const int domain_lo_1 = domain_lo[1];
    const int domain_lo_2 = domain_lo[2];
    const int domain_dim_0 = domain_dims[0];
    const int domain_dim_1 = domain_dims[1];
    const int domain_dim_2 = domain_dims[2];
    
    const int num_ghosts_0 = num_ghosts[0];
    const int num_ghosts_1 = num_ghosts[1]; 
    const int num_ghosts_2 = num_ghosts[2];
    const int ghostcell_dim_0 = ghostcell_dims[0];
    
    /************************************************
     * Set the immersed boundary variables from here.
     ************************************************/
    
    /*
     * Set the parameters of the sphere here.
     */
    
    const double half = double(1)/double(2);
    
    /*
     * These will be read from the input file.
     */
    
    Real radius_c = half; //double(20); AFK radius 
    Real x_c = Real(1); //half;  AFK center
    Real y_c = Real(1); //half;  AFK center 
    Real z_c = Real(1); //half;  AFK center

    if (d_initial_conditions_db != nullptr) // only if the inital condition files are unavilable. 
            {
                TBOX_ASSERT(d_initial_conditions_db->keyExists("x_c"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("y_c"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("z_c"));
                x_c     = d_initial_conditions_db->getReal("x_c");
                y_c     = d_initial_conditions_db->getReal("y_c");
                z_c     = d_initial_conditions_db->getReal("z_c");
                radius_c = d_initial_conditions_db->getReal("radius"); // radius_c ?
            }


    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++) // y
    {
        HAMERS_PRAGMA_SIMD // unncessary??
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++) // x
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++) // z
                // Compute the linear index. Get from Dr. Aslangil later
                {
                const int idx = (i + num_ghosts_0) +
                    (j + num_ghosts_1)*ghostcell_dim_0;
                
                // Compute the coordinates.
                double x[3]; // was 2, now 3???
                x[0] = patch_xlo[0] + (double(i) + double(1)/double(2))*dx[0]; // x coordinates of the point.
                x[1] = patch_xlo[1] + (double(j) + double(1)/double(2))*dx[1]; // y coordinates of the point.
                x[2] = patch_xlo[2] + (double(k) + double(1)/double(2))*dx[2]; // z coordinates of the point.

                // Distance from the sphere center.
                const double radius = sqrt(pow(x[0] - x_c, 2) + pow(x[1] - y_c, 2) + pow(x[2] - z_c, 2));
                // Angle between x axis and a line passing through center and current cell.
                const double theta  = atan2(x[1] - y_c, x[0] - x_c)*(180.0/M_PI);
                // Angle between the positive z-axis and the vector formed by the origin and the point. 0 < phi < pi
                const double phi = acos(x[2]/(sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2))))*(180.0/M_PI);
                
                if (radius < radius_c) // Condition that should be satisfied to be in sphere.
                {
                    double x_p = double(0);
                    double y_p = double(0);
                    double z_p = double(0);
                
                    if (x[2] > z_c)
                    {   
                        z_p = radius_c*cos(phi); 
                                        // positive z : (x,y,+)

                        if (x[0] > x_c) // positive x : (+,y,+)
                        {
                            x_p = radius_c*sin(phi)*cos(theta);
                        }
                        else            // negative x : (-,y,+)
                        {
                            x_p = -radius_c*sin(phi)*cos(theta);
                        }
                        
                        if (x[1] > y_c) // positive y : octant I    and octant II : (+,+,+) / (-,+,+)
                        {
                            y_p = radius_c*sin(phi)*sin(theta);
                        }
                        else            // negative y : octant III  and octant IV : (+,-,+) / (-,-,+)
                        {
                            y_p = -radius_c*sin(phi)*sin(theta);
                        }
                    }
                    else
                    {
                        z_p = -radius_c*cos(phi); 
                                        // negative z : (x,y,-)

                        if (x[0] > x_c) // positive x : (+,y,-)
                        {
                            x_p = radius_c*sin(phi)*cos(theta);
                        }
                        else            // negative x : (-,y,-)
                        {
                            x_p = -radius_c*sin(phi)*cos(theta);
                        }
                        
                        if (x[1] > y_c) // positive y : octant V    and octant VI  :  (+,+,-) / (-,+,-)
                        {
                            y_p = radius_c*sin(phi)*sin(theta);
                        }
                        else            // negative y : octant VII  and octant VIII : (-,-,-) / (+,-,-)
                        {
                            y_p = -radius_c*sin(phi)*sin(theta);
                        }
                    }

                    // is this conservative enough??/
                    if ((fabs(x_p - x[0]) < (double(d_num_immersed_boundary_ghosts[0]))*dx[0]) // ghost cell or body cell
                        (fabs(y_p - x[1]) < (double(d_num_immersed_boundary_ghosts[1]))*dx[1])
                        (fabs(z_p - x[2]) < (double(d_num_immersed_boundary_ghosts[2]))*dx[2]))

                    {
                        mask[idx]   = int(IB_MASK::IB_GHOST); /* basically takes the cell and finds the distance to see if its a ghost cell. makes it a ghost cell*/
                        dist[idx]   = radius_c - radius; // wall distance
                        norm_0[idx] = (x[0] - x_c)/radius; // normal vector from sphere center to cell center 
                        norm_1[idx] = (x[1] - y_c)/radius;
                        norm_2[idx] = (x[2] - z_c)/radius;

                    }
                    else
                    {
                        mask[idx]   = int(IB_MASK::BODY); /* if not its a body cell */
                        dist[idx]   = double(0);
                        norm_0[idx] = double(0);
                        norm_1[idx] = double(0);
                        norm_2[idx] = double(0);
                    }
                }
                else // fluid cell and can ignore.
                {
                    mask[idx]   = int(IB_MASK::FLUID);
                    dist[idx]   = double(0);
                    norm_0[idx] = double(0);
                    norm_1[idx] = double(0);
                    norm_2[idx] = double(0);
                }
                }
        }   
    }
}
