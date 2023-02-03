amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
{
    ParticleType& p = pstruct[i];
    for(int pc = 0; pc<p.idata(intData::num_comp_sphere); pc++){
        Real rp = p.rdata(realData::radius);
        
        Real pos[THREEDIM];
        Real pos_body[THREEDIM];
        get_inertial_body_fixed_pos(p, pc, pos, pos_body);
    
        bool collided=false;
        Real walldist=1e20;
        Real normal[AMREX_SPACEDIM]={1.0,0.0,0.0};
        Real fn[3],ft[3],rollfriqtorq_by_d[3],contact_vect[3];
        Real tow[3];
        int dir=XDIR;

        if (x_lo_bc==SOFTWALL_BC and pos[XDIR] < (plo[XDIR]+rp))
        {
            dir=XDIR;
            walldist=pos[dir]-plo[dir];
            collided=true;
            normal[XDIR]=0.0;
            normal[YDIR]=0.0;
            normal[ZDIR]=0.0;
            normal[dir]=-1.0;
        }
        else if (x_hi_bc==SOFTWALL_BC and pos[XDIR] > (phi[XDIR]-rp))
        {
            dir=XDIR;
            walldist=phi[dir]-pos[dir];
            collided=true;
            normal[XDIR]=0.0;
            normal[YDIR]=0.0;
            normal[ZDIR]=0.0;
            normal[dir]=1.0;
        }
        else if (y_lo_bc==SOFTWALL_BC and pos[YDIR] < (plo[YDIR]+rp))
        {
            dir=YDIR;
            walldist=pos[dir]-plo[dir];
            collided=true;
            normal[XDIR]=0.0;
            normal[YDIR]=0.0;
            normal[ZDIR]=0.0;
            normal[dir]=-1.0;
        }
        else if (y_hi_bc==SOFTWALL_BC and pos[YDIR] > (phi[YDIR]-rp))
        {
            dir=YDIR;
            walldist=phi[dir]-pos[dir];
            collided=true;
            normal[XDIR]=0.0;
            normal[YDIR]=0.0;
            normal[ZDIR]=0.0;
            normal[dir]=1.0;
        }
        else if (z_lo_bc==SOFTWALL_BC and pos[ZDIR] < (plo[ZDIR]+rp))
        {
            dir=ZDIR;
            walldist=pos[dir]-plo[dir];
            collided=true;
            normal[XDIR]=0.0;
            normal[YDIR]=0.0;
            normal[ZDIR]=0.0;
            normal[dir]=-1.0;
        }
        else if (z_hi_bc==SOFTWALL_BC and pos[ZDIR] > (phi[ZDIR]-rp))
        {
            dir=ZDIR;
            walldist=phi[dir]-pos[dir];
            collided=true;
            normal[XDIR]=0.0;
            normal[YDIR]=0.0;
            normal[ZDIR]=0.0;
            normal[dir]=1.0;
        }
        else
        {
            collided=false;
        }
        
        if(collided)
        {
            get_wall_forces(walldist,normal,p,dt,fn,ft,rollfriqtorq_by_d,pos,contact_vect);

            // each particle updates its force (no need for atomics)
            p.rdata(realData::fx) += fn[XDIR] + ft[XDIR];
            p.rdata(realData::fy) += fn[YDIR] + ft[YDIR];
            p.rdata(realData::fz) += fn[ZDIR] + ft[ZDIR];

            Real ftot[THREEDIM] = {fn[XDIR]+ft[XDIR], fn[YDIR]+ft[YDIR], fn[ZDIR]+ft[ZDIR]};
            Real torq[THREEDIM];
            crosspdt(contact_vect, ftot, torq);
            p.rdata(realData::taux) += torq[XDIR] + walldist*rollfriqtorq_by_d[XDIR];
            p.rdata(realData::tauy) += torq[YDIR] + walldist*rollfriqtorq_by_d[YDIR];
            p.rdata(realData::tauz) += torq[ZDIR] + walldist*rollfriqtorq_by_d[ZDIR];

            //doing adiabatic wall for now at Cartesian walls
            //if(do_heat_transfer)
            //{
            //    Real walltemp=get_wall_temperature(p.pos(0),p.pos(1),p.pos(2),walltemp_vardir,walltemp_polynomial);
            //    p.rdata(realData::temperature)=0.5*(walltemp+p.rdata(realData::temperature));
            //}
        }
    }
});
