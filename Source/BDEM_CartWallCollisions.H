amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
{
    ParticleType& p = pstruct[i];
    Real rp = p.rdata(realData::radius);

    bool collided=false;
    Real walldist=1e20;
    Real normal[AMREX_SPACEDIM]={1.0,0.0,0.0};
    Real fn[3],ft[3],nFr[3];
    Real tow[3];
    int dir=XDIR;

    if (x_lo_bc==SOFTWALL_BC and p.pos(XDIR) < (plo[XDIR]+rp))
    {
        dir=XDIR;
        walldist=p.pos(dir)-plo[dir];
        collided=true;
        normal[XDIR]=0.0;
        normal[YDIR]=0.0;
        normal[ZDIR]=0.0;
        normal[dir]=-1.0;
    }
    else if (x_hi_bc==SOFTWALL_BC and p.pos(XDIR) > (phi[XDIR]-rp))
    {
        dir=XDIR;
        walldist=phi[dir]-p.pos(dir);
        collided=true;
        normal[XDIR]=0.0;
        normal[YDIR]=0.0;
        normal[ZDIR]=0.0;
        normal[dir]=1.0;
    }
    else if (y_lo_bc==SOFTWALL_BC and p.pos(YDIR) < (plo[YDIR]+rp))
    {
        dir=YDIR;
        walldist=p.pos(dir)-plo[dir];
        collided=true;
        normal[XDIR]=0.0;
        normal[YDIR]=0.0;
        normal[ZDIR]=0.0;
        normal[dir]=-1.0;
    }
    else if (y_hi_bc==SOFTWALL_BC and p.pos(YDIR) > (phi[YDIR]-rp))
    {
        dir=YDIR;
        walldist=phi[dir]-p.pos(dir);
        collided=true;
        normal[XDIR]=0.0;
        normal[YDIR]=0.0;
        normal[ZDIR]=0.0;
        normal[dir]=1.0;
    }
    else if (z_lo_bc==SOFTWALL_BC and p.pos(ZDIR) < (plo[ZDIR]+rp))
    {
        dir=ZDIR;
        walldist=p.pos(dir)-plo[dir];
        collided=true;
        normal[XDIR]=0.0;
        normal[YDIR]=0.0;
        normal[ZDIR]=0.0;
        normal[dir]=-1.0;
    }
    else if (z_hi_bc==SOFTWALL_BC and p.pos(ZDIR) > (phi[ZDIR]-rp))
    {
        dir=ZDIR;
        walldist=phi[dir]-p.pos(dir);
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
        get_wall_forces(walldist,normal,p,dt,fn,ft,nFr);

        // each particle updates its force (no need for atomics)
        p.rdata(realData::fx) += fn[XDIR] + ft[XDIR];
        p.rdata(realData::fy) += fn[YDIR] + ft[YDIR];
        p.rdata(realData::fz) += fn[ZDIR] + ft[ZDIR];

        crosspdt(normal, ft, tow);

        p.rdata(realData::taux) += walldist*(tow[XDIR] + nFr[XDIR]);
        p.rdata(realData::tauy) += walldist*(tow[YDIR] + nFr[YDIR]);
        p.rdata(realData::tauz) += walldist*(tow[ZDIR] + nFr[ZDIR]);

        //doing adiabatic wall for now at Cartesian walls
        //if(do_heat_transfer)
        //{
        //    Real walltemp=get_wall_temperature(p.pos(0),p.pos(1),p.pos(2),walltemp_vardir,walltemp_polynomial);
        //    p.rdata(realData::temperature)=0.5*(walltemp+p.rdata(realData::temperature));
        //}
    }

});
