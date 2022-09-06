const auto phiarr = lsmfab->array(mfi);

amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
{
    ParticleType& p = pstruct[i];
    Real rp = p.rdata(realData::radius);

    Real ls_value = interp_level_set(p, ls_refinement, phiarr, plo, dxi);

    if (ls_value < rp)
    {
        Real normal[3],fn[3],ft[3],nFr[3];

        level_set_normal(p, ls_refinement, &normal[0], phiarr, plo, dxi);
        
        /*if(amrex::Math::abs(ls_value)>rp)
        {
            amrex::Print()<<"mindist,rad from lsets:"<<ls_value<<"\t"<<rp
            <<"\t"<<p.pos(0)<<"\t"<<p.pos(1)
            <<"\t"<<p.pos(2)<<"\t"<<normal[0]<<"\t"<<normal[1]<<"\t"<<normal[2]<<"\n";
        }*/

        normal[0] *= -1;
        normal[1] *= -1;
        normal[2] *= -1;

        get_wall_forces(ls_value,normal,p,dt,fn,ft,nFr);

        // each particle updates its force (no need for atomics)
        p.rdata(realData::fx) += fn[XDIR] + ft[XDIR];
        p.rdata(realData::fy) += fn[YDIR] + ft[YDIR];
        p.rdata(realData::fz) += fn[ZDIR] + ft[ZDIR];

        Real tow[3];

        cross_product(normal, ft, tow);

        p.rdata(realData::taux) += ls_value*(tow[XDIR] + nFr[XDIR]);
        p.rdata(realData::tauy) += ls_value*(tow[YDIR] + nFr[YDIR]);
        p.rdata(realData::tauz) += ls_value*(tow[ZDIR] + nFr[ZDIR]);

        if(do_heat_transfer)
        {
            Real walltemp=get_wall_temperature(p.pos(0),p.pos(1),p.pos(2),walltemp_vardir,walltemp_polynomial);
            p.rdata(realData::temperature)=0.5*(walltemp+p.rdata(realData::temperature));
        }

    }

});
