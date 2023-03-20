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

        Real rad_sum = p1.rdata(realData::radius) + p2.rdata(realData::radius);

        if ( dist2 <= (rad_sum - TINYVAL)*(rad_sum - TINYVAL) and (p1.id() != p2.id()))
        {
            Real distmag     = sqrt(dist2);

            AMREX_ASSERT(distmag >= TINYVAL);

            Real distmag_inv = one/distmag;

            Real normal[THREEDIM];
            normal[XDIR] = dist_x*distmag_inv;
            normal[YDIR] = dist_y*distmag_inv;
            normal[ZDIR] = dist_z*distmag_inv;

            Real overlap_n = rad_sum - distmag;
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

            // Contact law: default - linear, 1 - hertz
            Real kn_dem;
            Real eta_n_dem;
            Real kt_dem;
            Real eta_t_dem;
            if (contact_law == CONTACT_HERTZ)
            {
                Real E1 = p1.rdata(realData::E);
                Real E2 = p2.rdata(realData::E);
                Real nu1 = p1.rdata(realData::nu);
                Real nu2 = p2.rdata(realData::nu);
                Real E_eff = E1*E2/((1.0-nu1*nu1)*E2+(1.0-nu2*nu2)*E1);
                Real temp = 2.0*E2*(2.0-nu1)*(1.0+nu1) + 2.0*E1*(2.0-nu2)*(1.0+nu2);
                Real G_eff = E1*E2/temp;
                Real beta = log(DEM::e_n) / sqrt(PI*PI + log(DEM::e_n)*log(DEM::e_n));
                Real temp1 = sqrt(Reff*overlap_n);
                kn_dem = 4.0/3.0*E_eff*temp1;
                kt_dem = 8.0*G_eff*temp1;
                eta_n_dem = 2.0*sqrt(5.0/3.0)*beta*sqrt(meff*E_eff)*sqrt(temp1);
                eta_t_dem = 4.0*sqrt(5.0/3.0)*beta*sqrt(meff*G_eff)*sqrt(temp1);
            }
            else
            {            
                kn_dem = DEM::k_n;
                eta_n_dem = two*sqrt(DEM::k_n*meff)*fabs(log(DEM::e_n));
                eta_n_dem /= sqrt(PI*PI + log(DEM::e_n)*log(DEM::e_n));

                kt_dem = DEM::k_t;
                eta_t_dem = two*sqrt(DEM::k_t*meff)*fabs(log(DEM::e_t));
                eta_t_dem /= sqrt(PI*PI + log(DEM::e_t)*log(DEM::e_t));
            }

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
                // Real overlap_t_tan = dotpdt(overlap_t,tangent);
                // Real vrel_t_tan = dotpdt(vrel_t,tangent);

                ftsp[XDIR] = -(kt_dem*overlap_t[XDIR]
                               + eta_t_dem*vrel_t[XDIR]);
                ftsp[YDIR] = -(kt_dem*overlap_t[YDIR]
                               + eta_t_dem*vrel_t[YDIR]);
                ftsp[ZDIR] = -(kt_dem*overlap_t[ZDIR]
                               + eta_t_dem*vrel_t[ZDIR]);
                // ftsp[XDIR] = -(kt_dem*overlap_t_tan*tangent[XDIR]
                //                + eta_t_dem*vrel_t_tan*tangent[XDIR]);
                // ftsp[YDIR] = -(kt_dem*overlap_t_tan*tangent[YDIR]
                //                + eta_t_dem*vrel_t_tan*tangent[YDIR]);
                // ftsp[ZDIR] = -(kt_dem*overlap_t_tan*tangent[ZDIR]
                //                + eta_t_dem*vrel_t_tan*tangent[ZDIR]);
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


            p1.rdata(realData::fx) += fn[XDIR] + ft[XDIR];
            p1.rdata(realData::fy) += fn[YDIR] + ft[YDIR];
            p1.rdata(realData::fz) += fn[ZDIR] + ft[ZDIR];

            // For debug
            // p1.rdata(realData::overlap_n) = overlap_n;
            // p1.rdata(realData::kn) = kn_dem;
            // p1.rdata(realData::kt) = kt_dem;
            // p1.rdata(realData::eta_n) = eta_n_dem;
            // p1.rdata(realData::eta_t) = eta_t_dem;
            // p1.rdata(realData::overlap_n_y) = overlap_n*normal[YDIR];
            // p1.rdata(realData::fn_y) = fn[YDIR];
            // p1.rdata(realData::vr_n_y) = vrel_trans_norm*normal[YDIR];
            // p1.rdata(realData::overlap_t_y) = overlap_t[YDIR];
            // p1.rdata(realData::ft_y) = ft[YDIR];
            // p1.rdata(realData::vr_t_y) = vrel_t[YDIR];


            //Torque calculation
            //=================
            
           /*           /\
                     r1/   \r2
                      /-x----\
                      d1   d2

            d1+d2=d
            d1=r1*cos(theta)=r1*(d^2+r1^2-r^2)/(2*d*r1) */

            Real dist_c_p1 = 0.5 * (distmag + (p1.rdata(realData::radius)*p1.rdata(realData::radius) 
                            - p2.rdata(realData::radius)*p2.rdata(realData::radius))*distmag_inv);
            
            Real omega_p1[THREEDIM];
            omega_p1[XDIR] = p1.rdata(realData::xangvel);
            omega_p1[YDIR] = p1.rdata(realData::yangvel);
            omega_p1[ZDIR] = p1.rdata(realData::zangvel);
            Real mag_omegap1 = sqrt(dotpdt(omega_p1,omega_p1));
            Real rollfric_torq[THREEDIM]={0.0};
            if (mag_omegap1 > zero)
            {
                Real mag_fn = sqrt(dotpdt(fn,fn));
                Real torq   = -DEM::muR * mag_fn * dist_c_p1;
                rollfric_torq[XDIR]=torq*omega_p1[XDIR]/mag_omegap1;
                rollfric_torq[YDIR]=torq*omega_p1[YDIR]/mag_omegap1;
                rollfric_torq[ZDIR]=torq*omega_p1[ZDIR]/mag_omegap1;
            }

            Real tow[THREEDIM];
            crosspdt(normal, ft, tow);
            //note: tow here needs a distance, so done in final calc

            p1.rdata(realData::taux) += dist_c_p1*tow[XDIR]+rollfric_torq[XDIR];
            p1.rdata(realData::tauy) += dist_c_p1*tow[YDIR]+rollfric_torq[YDIR];
            p1.rdata(realData::tauz) += dist_c_p1*tow[ZDIR]+rollfric_torq[ZDIR];

            Real avgtemp=0.5*(p1.rdata(realData::temperature) + p2.rdata(realData::temperature));
            p1.rdata(realData::temperature)=avgtemp;
        }
    }

});
