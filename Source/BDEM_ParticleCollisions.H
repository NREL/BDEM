amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
{
    ParticleType& p1 = pstruct[i];

    for (const auto& p2 : nbor_data.getNeighbors(i))
    {                
        Real dist_x = p2.pos(0) - p1.pos(0);
        Real dist_y = p2.pos(1) - p1.pos(1);
        Real dist_z = p2.pos(2) - p1.pos(2);

        Real r2 = dist_x*dist_x +
            dist_y*dist_y +
            dist_z*dist_z;

        Real r_lm = p1.rdata(realData::radius) + p2.rdata(realData::radius);

        if ( r2 <= (r_lm - TINYVAL)*(r_lm - TINYVAL) and (p1.id() != p2.id()))
        {
            Real dist_mag     = sqrt(r2);

            AMREX_ASSERT(dist_mag >= TINYVAL);

            Real dist_mag_inv = one/dist_mag;

            Real normal[THREEDIM];
            normal[XDIR] = dist_x * dist_mag_inv;
            normal[YDIR] = dist_y * dist_mag_inv;
            normal[ZDIR] = dist_z * dist_mag_inv;

            Real overlap_n = r_lm - dist_mag;
            Real vrel_trans_norm;
            Real vrel_t[THREEDIM];
	    Real vRoll[THREEDIM];

            cfrelvel(p1, p2, vrel_trans_norm, vrel_t, normal, dist_mag);

            Real m1 = p1.rdata(realData::mass);
            Real m2 = p2.rdata(realData::mass);
            Real meff = m1*m2/(m1+m2);
	    Real rad1 = p1.rdata(realData::radius);
	    Real rad2 = p2.rdata(realData::radius);
	    Real Reff = rad1*rad2/(rad1+rad2);

            Real kn_des = DEM::k_n;
            Real eta_n_des = two*sqrt(DEM::k_n*meff)*fabs(log(DEM::e_n));
            eta_n_des /= sqrt(PI*PI + log(DEM::e_n)*log(DEM::e_n));

            // tangential components not implemented yet
            Real kt_des = DEM::k_t;
            Real eta_t_des = two*sqrt(DEM::k_t*meff)*fabs(log(DEM::e_t));
	    eta_t_des /= sqrt(PI*PI + log(DEM::e_t)*log(DEM::e_t));

            Real fn[THREEDIM];
            Real ft[THREEDIM];
            Real overlap_t[THREEDIM];
            Real mag_overlap_t;

            // calculate the normal contact force
            fn[XDIR] = -(kn_des*overlap_n*normal[XDIR]
                    + eta_n_des*vrel_trans_norm*normal[XDIR]);
            fn[YDIR] = -(kn_des*overlap_n*normal[YDIR]
                    + eta_n_des*vrel_trans_norm*normal[YDIR]);
            fn[ZDIR] = -(kn_des*overlap_n*normal[ZDIR]
                    + eta_n_des*vrel_trans_norm*normal[ZDIR]);

	    
	    Real tangent[THREEDIM];

            // calculate the tangential overlap
            overlap_t[XDIR] = dt*vrel_t[XDIR];
            overlap_t[YDIR] = dt*vrel_t[YDIR];
            overlap_t[ZDIR] = dt*vrel_t[ZDIR];
	    mag_overlap_t = sqrt(dot_product(overlap_t, overlap_t));
	    
	    

	    if (mag_overlap_t > zero)
            {
                Real fnmd = DEM::mu * sqrt(dot_product(fn, fn));
		Real ftsp[THREEDIM];
		Real ftn[THREEDIM];
		Real mag_ftsp;
		Real mag_ftn;
		tangent[XDIR] = overlap_t[XDIR]/mag_overlap_t;
		tangent[YDIR] = overlap_t[YDIR]/mag_overlap_t;
		tangent[ZDIR] = overlap_t[ZDIR]/mag_overlap_t;
		// calculate the tangential contact force
		Real overlap_t_tan = dot_product(overlap_t,tangent);
		Real vrel_t_tan = dot_product(vrel_t,tangent);
		
		ftsp[XDIR] = -(kt_des*overlap_t_tan*tangent[XDIR]
			     + eta_t_des*vrel_t_tan*tangent[XDIR]);
		ftsp[YDIR] = -(kt_des*overlap_t_tan*tangent[YDIR]
			     + eta_t_des*vrel_t_tan*tangent[YDIR]);
		ftsp[ZDIR] = -(kt_des*overlap_t_tan*tangent[ZDIR]
			     + eta_t_des*vrel_t_tan*tangent[ZDIR]);
                ftn[XDIR] = -fnmd * tangent[XDIR];
                ftn[YDIR] = -fnmd * tangent[YDIR];
                ftn[ZDIR] = -fnmd * tangent[ZDIR];
		mag_ftsp = sqrt(dot_product(ftsp,ftsp));
		mag_ftn = sqrt(dot_product(ftn,ftn));
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
	    Real mag_omegap1 = sqrt(dot_product(omega_p1,omega_p1));
	    Real nFr[THREEDIM];
	    Real mag_fn = sqrt(dot_product(fn,fn));
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
		cross_product(normal, Froll, nFr);
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

            Real dist_cl = 0.5 * (dist_mag + (part1_r*part1_r - part2_r*part2_r) * dist_mag_inv);
            dist_cl = dist_mag - dist_cl;

            Real tow[THREEDIM];

            cross_product(normal, ft, tow);

	    p1.rdata(realData::taux) += dist_cl*(tow[XDIR]+nFr[XDIR]);
	    p1.rdata(realData::tauy) += dist_cl*(tow[YDIR]+nFr[YDIR]);
	    p1.rdata(realData::tauz) += dist_cl*(tow[ZDIR]+nFr[ZDIR]);

            Real avgtemp=0.5*(p1.rdata(realData::temperature) + p2.rdata(realData::temperature));
            p1.rdata(realData::temperature)=avgtemp;
        }
    }

});
