amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE( int i ) noexcept {
    ParticleType &p = pstruct[i];
    Real rp         = p.rdata( realData::radius );

    Real pos[THREEDIM];
    pos[XDIR] = p.pos( 0 );
    pos[YDIR] = p.pos( 1 );
    pos[ZDIR] = p.pos( 2 );

    bool collided               = false;
    Real walldist               = 1e20;
    Real normal[AMREX_SPACEDIM] = { 1.0, 0.0, 0.0 };
    Real fn[3], ft[3], rollfriqtorq_by_d[3], contact_vect[3];
    Real tow[3];
    int dir = XDIR;

    if ( x_lo_bc == SOFTWALL_BC and pos[XDIR] < ( plo[XDIR] + rp ) )
    {
        dir          = XDIR;
        walldist     = pos[dir] - plo[dir];
        collided     = true;
        normal[XDIR] = 0.0;
        normal[YDIR] = 0.0;
        normal[ZDIR] = 0.0;
        normal[dir]  = -1.0;
        get_wall_forces(
            walldist, normal, p, dt, fn, ft, rollfriqtorq_by_d, pos, contact_law, liquid_bridging, 0, NULL, 0 );
    }
    else
    {
        p.rdata(realData::wall_toverlap    ) = 0.;
        p.rdata(realData::wall_toverlap + 1) = 0.;
        p.rdata(realData::wall_toverlap + 2) = 0.;
    }

    if ( x_hi_bc == SOFTWALL_BC and pos[XDIR] > ( phi[XDIR] - rp ) )
    {
        dir          = XDIR;
        walldist     = phi[dir] - pos[dir];
        collided     = true;
        normal[XDIR] = 0.0;
        normal[YDIR] = 0.0;
        normal[ZDIR] = 0.0;
        normal[dir]  = 1.0;
        get_wall_forces(
            walldist, normal, p, dt, fn, ft, rollfriqtorq_by_d, pos, contact_law, liquid_bridging, 0, NULL, 1 );
    }
    else
    {
        p.rdata(realData::wall_toverlap + 3) = 0.;
        p.rdata(realData::wall_toverlap + 4) = 0.;
        p.rdata(realData::wall_toverlap + 5) = 0.;
    }
    
    if ( y_lo_bc == SOFTWALL_BC and pos[YDIR] < ( plo[YDIR] + rp ) )
    {
        dir          = YDIR;
        walldist     = pos[dir] - plo[dir];
        collided     = true;
        normal[XDIR] = 0.0;
        normal[YDIR] = 0.0;
        normal[ZDIR] = 0.0;
        normal[dir]  = -1.0;
        get_wall_forces(
            walldist, normal, p, dt, fn, ft, rollfriqtorq_by_d, pos, contact_law, liquid_bridging, 0, NULL, 2 );
    }
    else
    {
        p.rdata(realData::wall_toverlap + 6) = 0.;
        p.rdata(realData::wall_toverlap + 7) = 0.;
        p.rdata(realData::wall_toverlap + 8) = 0.;
    }
    
    if ( y_hi_bc == SOFTWALL_BC and pos[YDIR] > ( phi[YDIR] - rp ) )
    {
        dir          = YDIR;
        walldist     = phi[dir] - pos[dir];
        collided     = true;
        normal[XDIR] = 0.0;
        normal[YDIR] = 0.0;
        normal[ZDIR] = 0.0;
        normal[dir]  = 1.0;
        get_wall_forces(
            walldist, normal, p, dt, fn, ft, rollfriqtorq_by_d, pos, contact_law, liquid_bridging, 0, NULL, 3 );
    }
    else
    {
        p.rdata(realData::wall_toverlap + 9) = 0.;
        p.rdata(realData::wall_toverlap + 10) = 0.;
        p.rdata(realData::wall_toverlap + 11) = 0.;
    }
    
    if ( z_lo_bc == SOFTWALL_BC and pos[ZDIR] < ( plo[ZDIR] + rp ) )
    {
        dir          = ZDIR;
        walldist     = pos[dir] - plo[dir];
        collided     = true;
        normal[XDIR] = 0.0;
        normal[YDIR] = 0.0;
        normal[ZDIR] = 0.0;
        normal[dir]  = -1.0;
        get_wall_forces(
            walldist, normal, p, dt, fn, ft, rollfriqtorq_by_d, pos, contact_law, liquid_bridging, 0, NULL, 4 );
    }
    else
    {
        p.rdata(realData::wall_toverlap + 12) = 0.;
        p.rdata(realData::wall_toverlap + 13) = 0.;
        p.rdata(realData::wall_toverlap + 14) = 0.;
    }
    
    if ( z_hi_bc == SOFTWALL_BC and pos[ZDIR] > ( phi[ZDIR] - rp ) )
    {
        dir          = ZDIR;
        walldist     = phi[dir] - pos[dir];
        collided     = true;
        normal[XDIR] = 0.0;
        normal[YDIR] = 0.0;
        normal[ZDIR] = 0.0;
        normal[dir]  = 1.0;
        get_wall_forces(
            walldist, normal, p, dt, fn, ft, rollfriqtorq_by_d, pos, contact_law, liquid_bridging, 0, NULL, 5 );
    }
    else
    {
        p.rdata(realData::wall_toverlap + 15) = 0.;
        p.rdata(realData::wall_toverlap + 16) = 0.;
        p.rdata(realData::wall_toverlap + 17) = 0.;
    }
    

    if ( collided )
    {

        // each particle updates its force (no need for atomics)

        // doing adiabatic wall for now at Cartesian walls
        // if(do_heat_transfer)
        //{
        //     Real
        //     walltemp=get_wall_temperature(p.pos(0),p.pos(1),p.pos(2),walltemp_vardir,walltemp_polynomial);
        //     p.rdata(realData::temperature)=0.5*(walltemp+p.rdata(realData::temperature));
        // }
    }
} );
