amrex::ParallelFor(np,[=]
                  AMREX_GPU_DEVICE (int i) noexcept
{
    ParticleType& p1 = pstruct[i];

    for (const auto& p2 : nbor_data.getNeighbors(i))
    {                
        Real dist_x = p2.pos(0) - p1.pos(0);
        Real dist_y = p2.pos(1) - p1.pos(1);
        Real dist_z = p2.pos(2) - p1.pos(2);

        Real dist2 = dist_x*dist_x +
        dist_y*dist_y +
        dist_z*dist_z;

        Real r_lm = p1.rdata(realData::radius) + p2.rdata(realData::radius);

        if ( dist2 <= (r_lm - TINYVAL)*(r_lm - TINYVAL) and (p1.id() != p2.id()))
        {
            Real distmag     = sqrt(dist2);

            AMREX_ASSERT(distmag >= TINYVAL);

            Real distmag_inv = one/distmag;

            Real normal[THREEDIM];
            normal[XDIR] = dist_x*distmag_inv;
            normal[YDIR] = dist_y*distmag_inv;
            normal[ZDIR] = dist_z*distmag_inv;

            Real overlap_n = r_lm - distmag;
            Real vrel_trans_norm;
            Real vrel_t[THREEDIM];
            Real vRoll[THREEDIM];

            get_relative_velocities(p1, p2, vrel_trans_norm, vrel_t, normal, distmag);

            Real m1 = p1.rdata(realData::mass);
            Real m2 = p2.rdata(realData::mass);
            Real meff = m1*m2/(m1+m2);
            Real rad1 = p1.rdata(realData::radius);
            Real rad2 = p2.rdata(realData::radius);
            Real Reff = rad1*rad2/(rad1+rad2);

            Real kn_dem = DEM::k_n;
            Real eta_n_dem = two*sqrt(DEM::k_n*meff)*fabs(log(DEM::e_n));
            eta_n_dem /= sqrt(PI*PI + log(DEM::e_n)*log(DEM::e_n));

            // tangential components not implemented yet
            Real kt_dem = DEM::k_t;
            Real eta_t_dem = two*sqrt(DEM::k_t*meff)*fabs(log(DEM::e_t));
            eta_t_dem /= sqrt(PI*PI + log(DEM::e_t)*log(DEM::e_t));

            Real fn[THREEDIM];
            Real ft[THREEDIM];
            Real overlap_t[THREEDIM];
            Real mag_overlap_t;

            // calculate the normal contact force
            fn[XDIR] = -(kn_dem*overlap_n*normal[XDIR]
                         + eta_n_dem*vrel_trans_norm*normal[XDIR]);
            fn[YDIR] = -(kn_dem*overlap_n*normal[YDIR]
                         + eta_n_dem*vrel_trans_norm*normal[YDIR]);
            fn[ZDIR] = -(kn_dem*overlap_n*normal[ZDIR]
                         + eta_n_dem*vrel_trans_norm*normal[ZDIR]);


            Real tangent[THREEDIM];

            // calculate the tangential overlap
            overlap_t[XDIR] = dt*vrel_t[XDIR];
            overlap_t[YDIR] = dt*vrel_t[YDIR];
            overlap_t[ZDIR] = dt*vrel_t[ZDIR];
            mag_overlap_t = sqrt(dotpdt(overlap_t, overlap_t));



            if (mag_overlap_t > zero)
            {
                Real fnmd = DEM::mu * sqrt(dotpdt(fn, fn));
                Real ftsp[THREEDIM];
                Real ftn[THREEDIM];
                Real mag_ftsp;
                Real mag_ftn;
                tangent[XDIR] = overlap_t[XDIR]/mag_overlap_t;
                tangent[YDIR] = overlap_t[YDIR]/mag_overlap_t;
                tangent[ZDIR] = overlap_t[ZDIR]/mag_overlap_t;
                // calculate the tangential contact force
                Real overlap_t_tan = dotpdt(overlap_t,tangent);
                Real vrel_t_tan = dotpdt(vrel_t,tangent);

                ftsp[XDIR] = -(kt_dem*overlap_t_tan*tangent[XDIR]
                               + eta_t_dem*vrel_t_tan*tangent[XDIR]);
                ftsp[YDIR] = -(kt_dem*overlap_t_tan*tangent[YDIR]
                               + eta_t_dem*vrel_t_tan*tangent[YDIR]);
                ftsp[ZDIR] = -(kt_dem*overlap_t_tan*tangent[ZDIR]
                               + eta_t_dem*vrel_t_tan*tangent[ZDIR]);
                ftn[XDIR] = -fnmd * tangent[XDIR];
                ftn[YDIR] = -fnmd * tangent[YDIR];
                ftn[ZDIR] = -fnmd * tangent[ZDIR];
                mag_ftsp = sqrt(dotpdt(ftsp,ftsp));
                mag_ftn = sqrt(dotpdt(ftn,ftn));
                if (mag_ftsp < mag_ftn)
                {
                    ft[XDIR] = ftsp[XDIR];
                    ft[YDIR] = ftsp[YDIR];
                    ft[ZDIR] = ftsp[ZDIR];
                }
                else
                {
                    ft[XDIR] = ftn[XDIR];
                    ft[YDIR] = ftn[YDIR];
                    ft[ZDIR] = ftn[ZDIR];
                }
            }  
            else 
            {
                ft[XDIR] = zero;
                ft[YDIR] = zero;
                ft[ZDIR] = zero;
            }

            Real omega_p1[THREEDIM];
            omega_p1[XDIR] = p1.rdata(realData::xangvel);
            omega_p1[YDIR] = p1.rdata(realData::yangvel);
            omega_p1[ZDIR] = p1.rdata(realData::zangvel);
            Real mag_omegap1 = sqrt(dotpdt(omega_p1,omega_p1));
            Real nFr[THREEDIM];
            Real mag_fn = sqrt(dotpdt(fn,fn));
            if (mag_omegap1 > zero)
            {
                Real vRatio[THREEDIM];
                Real Froll[THREEDIM];
                vRatio[XDIR] = omega_p1[XDIR]/mag_omegap1;
                vRatio[YDIR] = omega_p1[YDIR]/mag_omegap1;
                vRatio[ZDIR] = omega_p1[ZDIR]/mag_omegap1;
                Froll[XDIR] = -DEM::muR * mag_fn * vRatio[XDIR];
                Froll[YDIR] = -DEM::muR * mag_fn * vRatio[YDIR];
                Froll[ZDIR] = -DEM::muR * mag_fn * vRatio[ZDIR];
                crosspdt(normal, Froll, nFr);
            }
            else
            {
                nFr[XDIR] = zero;
                nFr[YDIR] = zero;
                nFr[ZDIR] = zero;
            }

            // each particle updates its force (no need for atomics)
            p1.rdata(realData::fx) += fn[XDIR] + ft[XDIR];
            p1.rdata(realData::fy) += fn[YDIR] + ft[YDIR];
            p1.rdata(realData::fz) += fn[ZDIR] + ft[ZDIR];

            Real part1_r = p1.rdata(realData::radius);
            Real part2_r = p2.rdata(realData::radius);

            Real dist_cl = 0.5 * (distmag + (part1_r*part1_r - part2_r*part2_r)*distmag_inv);
            dist_cl = distmag - dist_cl;

            Real tow[THREEDIM];

            crosspdt(normal, ft, tow);

            p1.rdata(realData::taux) += dist_cl*(tow[XDIR]+nFr[XDIR]);
            p1.rdata(realData::tauy) += dist_cl*(tow[YDIR]+nFr[YDIR]);
            p1.rdata(realData::tauz) += dist_cl*(tow[ZDIR]+nFr[ZDIR]);

            Real avgtemp=0.5*(p1.rdata(realData::temperature) + p2.rdata(realData::temperature));
            p1.rdata(realData::temperature)=avgtemp;
        }
    }

});
