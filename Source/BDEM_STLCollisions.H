amrex::ParallelFor(np,[=]
       AMREX_GPU_DEVICE (int i) noexcept
{
    Real pos[3],pos_t[3],dist;
    Real mat[9],transf_pts[9];
    int min_tri_id=0;
    bool success,collided;
    Real centpvec[3],norm[3];

    ParticleType& p = pstruct[i];
    pos[XDIR] = p.pos(0);
    pos[YDIR] = p.pos(1);
    pos[ZDIR] = p.pos(2);
    pos_t[XDIR] = p.pos(0);
    pos_t[YDIR] = p.pos(1);
    pos_t[ZDIR] = p.pos(2);
    int min_tr = -1;

    Real mindist=BIGVAL;
    if( (p.pos(0) + p.rdata(realData::radius) > x_lo_bb ) && 
        (p.pos(0) - p.rdata(realData::radius) < x_hi_bb ) &&
        (p.pos(1) + p.rdata(realData::radius) > y_lo_bb ) &&
        (p.pos(1)- p.rdata(realData::radius) < y_hi_bb ) &&
        (p.pos(2)+ p.rdata(realData::radius) > z_lo_bb ) &&
        (p.pos(2)- p.rdata(realData::radius) < z_hi_bb )
    )
    {
        int ndistcalcs;

        const Real movedir[3] = {x_movedir, y_movedir, z_movedir};
        const Real movecenter[3] = {x_movecenter, y_movecenter, z_movecenter};
        const Real grid_bbmin[3] = {x_grid_bbmin, y_grid_bbmin, z_grid_bbmin };
        const Real  grid_bbmax[3] = {x_grid_bbmax, y_grid_bbmax, z_grid_bbmax };
        const Real grid_delta[3] = {x_grid_delta, y_grid_delta, z_grid_delta };
        const int grid_size[3] = {x_grid_size, y_grid_size, z_grid_size};

        get_closest_local_tri(
            pos,
            mindist,
            min_tri_id,
            time - dt,
            movetype,
            movedir,
            movecenter,
            movevel,
            grid_bbmin,
            grid_bbmax,
            grid_delta,
            grid_size,
            tris_per_cell,
            tris_in_grid,
            cell_start,
            tri_pts      
        );

        min_tr = min_tri_id;

        if(min_tr != -1)
        {
            centpvec[0]=pos[0]-
                0.3333*(tri_pts[min_tr*9+0]+tri_pts[min_tr*9+3]+tri_pts[min_tr*9+6]);
            centpvec[1]=pos[1]-
                0.3333*(tri_pts[min_tr*9+1]+tri_pts[min_tr*9+4]+tri_pts[min_tr*9+7]);
            centpvec[2]=pos[2]-
                0.3333*(tri_pts[min_tr*9+2]+tri_pts[min_tr*9+5]+tri_pts[min_tr*9+8]);

            norm[0]=tri_normals[min_tr*3+0];
            norm[1]=tri_normals[min_tr*3+1];
            norm[2]=tri_normals[min_tr*3+2];

            Real dotpdt=norm[0]*centpvec[0] + norm[1]*centpvec[1] + norm[2]*centpvec[2];
            
            if(mindist <p.rdata(realData::radius))
            {
                norm[0]=-tri_normals[min_tr*3+0];
                norm[1]=-tri_normals[min_tr*3+1];
                norm[2]=-tri_normals[min_tr*3+2];
                
                Real fn[3],ft[3],rollfriqtorq_by_d[3],tow[3],contact_vect[3];

                // Compute wall velocity
                Real vw[3] = {0.,0.,0.};

                if(movetype == 1)
                {
                    vw[0] = movedir[0]*movevel;
                    vw[1] = movedir[1]*movevel;
                    vw[2] = movedir[2]*movevel;
                }
                else if (movetype == 2)
                {
                    // Get contact point
                    Real contact_point[3];
                    contact_point[0] = pos[0] + norm[0]*dotpdt - movecenter[0];
                    contact_point[1] = pos[1] + norm[1]*dotpdt - movecenter[1];
                    contact_point[2] = pos[2] + norm[2]*dotpdt - movecenter[2];
                
                    // Get angular velocity vector
                    Real omega[3];
                    omega[0] = movevel*movedir[0];
                    omega[1] = movevel*movedir[1];
                    omega[2] = movevel*movedir[2];

                    // Get velocity
                    CrossProd(omega,contact_point,vw);
                }

                get_wall_forces
                (
                    mindist,
                    norm,
                    p,
                    dt,
                    fn,
                    ft,
                    rollfriqtorq_by_d,
                    pos,
                    contact_law, 
                    liquid_bridging, 
                    movetype, 
                    vw
                );
    
    #ifndef AMREX_USE_GPU

                pressure[min_tr] += sqrt(DotProd(fn,fn));
                shear_stress[3*min_tr]      += ft[0];
                shear_stress[3*min_tr + 1]  += ft[1];
                shear_stress[3*min_tr + 2]  += ft[2];
    #else
                Gpu::Atomic::Add( &(pressure[min_tr]), sqrt(fn[0]*fn[0] + fn[1]*fn[1] + fn[2]*fn[2]) );
                Gpu::Atomic::Add( &(shear_stress[3*min_tr]), ft[0] );
                Gpu::Atomic::Add( &(shear_stress[3*min_tr +1]), ft[1] );
                Gpu::Atomic::Add( &(shear_stress[3*min_tr +2]), ft[2] );
    #endif
                if(solve_fibrillation == 1)
                {
                    Real vp[3] = { p.rdata(realData::xvel), p.rdata(realData::yvel), p.rdata(realData::zvel) };
                    Real dx[3] = { norm[0]*mindist + TINYVAL, norm[1]*mindist + TINYVAL, norm[2]*mindist + TINYVAL};
                    Real fc[3]; // No cohesion with the wall... it doesn't fibrillize
                    Real p2xi = 0.;
                    solve_fibrils_equation
                    (
                        dt,
                        vp,
                        vw,
                        dx,
                        p.rdata(realData::fraction_of_fibrils),
                        p2xi,
                        0.,
                        fc
                    );
                }

            } 
        }
    }
});
