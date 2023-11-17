const auto phiarr = lsmfab->array(mfi);

amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
{
    ParticleType& p = pstruct[i];
    Real rp = p.rdata(aos_realData::radius);

    Real ppos_inert[THREEDIM];
    ppos_inert[XDIR] = p.pos(0);
    ppos_inert[YDIR] = p.pos(1);
    ppos_inert[ZDIR] = p.pos(2);

    Real ls_value = get_levelset_value(ppos_inert, ls_refinement, phiarr, plo, dx);

    if (ls_value < rp)
    {
        Real normal[3],fn[3],ft[3],rollfriqtorq_by_d[3],contact_vect[3];

        get_levelset_normal(ppos_inert,ls_refinement,phiarr,plo,dx,normal);

        /*if(amrex::Math::abs(ls_value)>rp)
          {
          amrex::Print()<<"mindist,rad from lsets:"<<ls_value<<"\t"<<rp
          <<"\t"<<p.pos(0)<<"\t"<<p.pos(1)
          <<"\t"<<p.pos(2)<<"\t"<<normal[0]<<"\t"<<normal[1]<<"\t"<<normal[2]<<"\n";
          }*/

        normal[0] *= -1;
        normal[1] *= -1;
        normal[2] *= -1;

        Real ftau_arr[6];
        get_wall_forces(ls_value,normal,p,dt,fn,ft,rollfriqtorq_by_d,ppos_inert,contact_law, liquid_bridging, ftau_arr);

        fx_arr[i] += ftau_arr[0];
        fy_arr[i] += ftau_arr[1];
        fz_arr[i] += ftau_arr[2];
        taux_arr[i] += ftau_arr[3];
        tauy_arr[i] += ftau_arr[4];
        tauz_arr[i] += ftau_arr[5];

        if(do_heat_transfer)
        {
            Real walltemp=get_wall_temperature(p.pos(0),p.pos(1),p.pos(2),walltemp_vardir,walltemp_polynomial);
            p.rdata(aos_realData::temperature)=0.5*(walltemp+p.rdata(aos_realData::temperature));
        }

    }
});
