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
    if( (p.pos(0) + p.rdata(realData::radius) > stlptr->bbox_lo[0]) && 
        (p.pos(0) - p.rdata(realData::radius) < stlptr->bbox_hi[0]) &&
        (p.pos(1) + p.rdata(realData::radius) > stlptr->bbox_lo[1]) &&
        (p.pos(1)- p.rdata(realData::radius) < stlptr->bbox_hi[1]) &&
        (p.pos(2)+ p.rdata(realData::radius) > stlptr->bbox_lo[2]) &&
        (p.pos(2)- p.rdata(realData::radius) < stlptr->bbox_hi[2]) )
    {
        int ndistcalcs;

    if ( movetype == 0 ) // This seems to work ok only for static stl files
    {
        stlptr->get_closest_local_tri(
            pos,
            mindist,
            min_tri_id,
            time - dt,
            movetype,
            movedir,
            movecenter,
            movevel
        );

        min_tr = min_tri_id;
    }
    else
    {
        stlptr->searchtriangulation(0,stlptr->num_tri-1,
                    stlptr->sorted_indexarray,
                    pos,mindist,min_tri_id,ndistcalcs);
        mindist=std::sqrt(mindist);
        min_tr=stlptr->sorted_indexarray[min_tri_id];
    }
//        std::printf("\nDistance check complete %i", min_tri_id);

        if(min_tr != -1)
        {
            centpvec[0]=pos[0]-
                0.3333*(stlptr->tri_pts[min_tr*9+0]+stlptr->tri_pts[min_tr*9+3]+stlptr->tri_pts[min_tr*9+6]);
            centpvec[1]=pos[1]-
                0.3333*(stlptr->tri_pts[min_tr*9+1]+stlptr->tri_pts[min_tr*9+4]+stlptr->tri_pts[min_tr*9+7]);
            centpvec[2]=pos[2]-
                0.3333*(stlptr->tri_pts[min_tr*9+2]+stlptr->tri_pts[min_tr*9+5]+stlptr->tri_pts[min_tr*9+8]);

            norm[0]=stlptr->tri_normals[min_tr*3+0];
            norm[1]=stlptr->tri_normals[min_tr*3+1];
            norm[2]=stlptr->tri_normals[min_tr*3+2];

            Real dotpdt=norm[0]*centpvec[0] + norm[1]*centpvec[1] + norm[2]*centpvec[2];
            mindist*=amrex::Math::copysign(1.0,dotpdt);   
            if(mindist <p.rdata(realData::radius))
            {
                norm[0]=-stlptr->tri_normals[min_tr*3+0];
                norm[1]=-stlptr->tri_normals[min_tr*3+1];
                norm[2]=-stlptr->tri_normals[min_tr*3+2];
                
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
                // stlptr->pressure[min_tr] += sqrt(fn[0]*fn[0] + fn[1]*fn[1] + fn[2]*fn[2]);
                // stlptr->shear_stress[3*min_tr]      += ft[0];
                // stlptr->shear_stress[3*min_tr + 1]  += ft[1];
                // stlptr->shear_stress[3*min_tr + 2]  += ft[2];
                if(solve_fibrillation == 1)
                {
                    Real vp[3] = { p.rdata(realData::xvel), p.rdata(realData::yvel), p.rdata(realData::zvel) };
                    Real dx[3] = { norm[0]*mindist + TINYVAL, norm[1]*mindist + TINYVAL, norm[2]*mindist + TINYVAL};
                    solve_fibrils_equation
                    (
                        dt,
                        vp,
                        vw,
                        dx,
                        p.rdata(realData::fraction_of_fibrils)
                    );
                }

            } 
        }
    }
});
