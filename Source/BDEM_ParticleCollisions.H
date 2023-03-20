amrex::ParallelFor(np,[=]
                  AMREX_GPU_DEVICE (int i) noexcept
{
    ParticleType& p1 = pstruct[i];

    for(int pc1 = 0; pc1 < p1.idata(intData::num_comp_sphere); pc1++){
        
        // Calculate center of mass of component spheres
        //    - Alternative to recalculating this for each force evaluation would be to store positions in a vector and update with each step
        //    - Adhoc position evaluation like this only an option for simple cases like rod-like glued sphere particle
        //    - For rod-like glued sphere particle, assumption is that spheres extend along the body-fixed x-axis
        Real p1pos_inert[THREEDIM];
        Real p1pos_body[THREEDIM];
        get_inertial_body_fixed_pos(p1, pc1, p1pos_inert, p1pos_body);

        // Calculate the gravitational force and nonzero torque
        p1.rdata(realData::fx) += gravity[XDIR]*p1.rdata(realData::mass)/p1.idata(intData::num_comp_sphere);
        p1.rdata(realData::fy) += gravity[YDIR]*p1.rdata(realData::mass)/p1.idata(intData::num_comp_sphere);
        p1.rdata(realData::fz) += gravity[ZDIR]*p1.rdata(realData::mass)/p1.idata(intData::num_comp_sphere);
  
        Real posvect[THREEDIM];     // Vector from CoM of particle to CoM of component sphere
        posvect[XDIR] = p1pos_inert[XDIR] - p1.pos(0);
        posvect[YDIR] = p1pos_inert[YDIR] - p1.pos(1);
        posvect[ZDIR] = p1pos_inert[ZDIR] - p1.pos(2);

        // Code for glued sphere spin test
        Real f_init[THREEDIM] = {0.0, 0.0, 0.0};
        f_init[init_force_dir] = init_force;
        if(pc1 == init_force_comp){
            p1.rdata(realData::fx) += f_init[XDIR];
            p1.rdata(realData::fy) += f_init[YDIR];
            p1.rdata(realData::fz) += f_init[ZDIR];

            Real torq_init[THREEDIM];
            crosspdt(posvect, f_init, torq_init);
            p1.rdata(realData::taux) += torq_init[XDIR];
            p1.rdata(realData::tauy) += torq_init[YDIR];
            p1.rdata(realData::tauz) += torq_init[ZDIR];
        }

        for (const auto& p2 : nbor_data.getNeighbors(i))
        {                
            for(int pc2 = 0; pc2 < p2.idata(intData::num_comp_sphere); pc2++){
                // Renaming properties for readability
                Real m1 = p1.rdata(realData::mass);
                Real m2 = p2.rdata(realData::mass);
                Real rad1 = p1.rdata(realData::radius);
                Real rad2 = p2.rdata(realData::radius);

                // Calculate positional quantities for p2
                Real p2pos_inert[THREEDIM];
                Real p2pos_body[THREEDIM];
                get_inertial_body_fixed_pos(p2, pc2, p2pos_inert, p2pos_body);

                // Determine if there is any component sphere overlap
                Real dist_x = p2pos_inert[XDIR] - p1pos_inert[XDIR];
                Real dist_y = p2pos_inert[YDIR] - p1pos_inert[YDIR];
                Real dist_z = p2pos_inert[ZDIR] - p1pos_inert[ZDIR];

                Real dist2 = dist_x*dist_x +
                dist_y*dist_y +
                dist_z*dist_z;

                // Note: assuming all component particles have the same radius
                Real rad_sum = rad1 + rad2;
      
                Real distmag     = sqrt(dist2);
                AMREX_ASSERT(distmag >= TINYVAL);
                Real distmag_inv = one/distmag;
            
                bool bridge_exists = false;
                int bridge_index = 0;
                if(liquid_bridging){
                    // Check to see if a liquid bridge exists between the component spheres
                    for(int br=0; br<MAXBRIDGES; br++){
                        if(p1.idata(intData::first_bridge+3*br+0) == p2.id()){
                            if(p1.idata(intData::first_bridge+3*br+1) == pc1 && p1.idata(intData::first_bridge+3*br+2) == pc2){
                                bridge_exists = true;
                                bridge_index = br;
                            }
                        }
                    }
                }
  
                // Initialize variables used for both contact and liquid bridge force calculations
                // Calculate effective quantities
                Real meff = m1*m2/(m1+m2);
                Real Reff = rad1*rad2/(rad1+rad2);
                Real normal[THREEDIM];
                Real tangent[THREEDIM] = {0.0};
                Real overlap_n = rad_sum - distmag;
                Real vrel_trans_norm;
                Real vrel_trans_tan = 0.0;
                Real vrel_n[THREEDIM];
                Real vrel_t[THREEDIM];
                Real contact_loc[THREEDIM];
                Real contact_dist;
                Real contact_vect1[THREEDIM];
                Real contact_vect2[THREEDIM];
                Real vrel[THREEDIM];
                Real angvel1[THREEDIM] = {p1.rdata(realData::xangvel), p1.rdata(realData::yangvel), p1.rdata(realData::zangvel)};
                Real angvel2[THREEDIM] = {p2.rdata(realData::xangvel), p2.rdata(realData::yangvel), p2.rdata(realData::zangvel)};
                Real cpdt1[THREEDIM];
                Real cpdt2[THREEDIM];

                if ( (dist2 <= (rad_sum - TINYVAL)*(rad_sum - TINYVAL) and (p1.id() != p2.id())) || bridge_exists )
                {
                    // Add liquid bridge between particles if there isn't one already, and at least one particle has moisture content
                    if(liquid_bridging && !bridge_exists && (p1.rdata(realData::liquid_volume) + p2.rdata(realData::liquid_volume) > 0.0)){
                        // Find first available bridge index (denoted with -1 value)
                        while(p1.idata(intData::first_bridge+bridge_index*3) != -1){
                            bridge_index++;
                            if(bridge_index == MAXBRIDGES) Abort("\nMaximum number of liquid bridges reached, MAXBRIDGES needs to be increased!\n");
                        }

                        // Calculate the liquid volume belonging to p1 that contributes to the new bridge, and add to bridge total
                        Real Vtot = (p1.rdata(realData::liquid_volume)/(p1.idata(intData::num_comp_sphere)*2.0))*(1.0-sqrt(1.0-(pow(rad2,two)/(pow(rad1,two)+pow(rad2,two)))));
    
                        // Only create new bridge if the particle has sufficient liquid volume to do so
                        if(p1.rdata(realData::total_bridge_volume) + Vtot < p1.rdata(realData::liquid_volume)){
                            p1.rdata(realData::total_bridge_volume) += Vtot;

                            // Fill in relevant bridge indices (p2 id, particle 1 component sphere, particle 2 component sphere)
                            p1.idata(intData::first_bridge+bridge_index*3 + 0) = p2.id();
                            p1.idata(intData::first_bridge+bridge_index*3 + 1) = pc1;
                            p1.idata(intData::first_bridge+bridge_index*3 + 2) = pc2;
                            bridge_exists = true;
                        }
                    }
                    
                    // Normal vector pointing from center of component sphere 2 to component sphere 2
                    normal[XDIR] = dist_x*distmag_inv;
                    normal[YDIR] = dist_y*distmag_inv;
                    normal[ZDIR] = dist_z*distmag_inv;

                    // Calculate the distance from the center of the two component spheres to the point of contact
                    contact_dist = (dist2 + pow(rad1,two) - pow(rad2,two)) / (2.0 * distmag);

                    // Calculate distance from origin to point of contact
                    contact_loc[XDIR] = p1pos_inert[XDIR] + contact_dist*normal[XDIR];
                    contact_loc[YDIR] = p1pos_inert[YDIR] + contact_dist*normal[YDIR];
                    contact_loc[ZDIR] = p1pos_inert[ZDIR] + contact_dist*normal[ZDIR];

                    // Calculate the vector pointing from the center of mass of each composite particle to the point of contact
                    contact_vect1[XDIR] = contact_loc[XDIR] - p1.pos(0);
                    contact_vect1[YDIR] = contact_loc[YDIR] - p1.pos(1);
                    contact_vect1[ZDIR] = contact_loc[ZDIR] - p1.pos(2);

                    contact_vect2[XDIR] = contact_loc[XDIR] - p2.pos(0);
                    contact_vect2[YDIR] = contact_loc[YDIR] - p2.pos(1);
                    contact_vect2[ZDIR] = contact_loc[ZDIR] - p2.pos(2);

                    // Calculate relative velocity of component spheres
                    crosspdt(angvel1, contact_vect1, cpdt1);
                    crosspdt(angvel2, contact_vect2, cpdt2);
  
                    vrel[XDIR] = p1.rdata(realData::xvel) - p2.rdata(realData::xvel) + cpdt1[XDIR] - cpdt2[XDIR];
                    vrel[YDIR] = p1.rdata(realData::yvel) - p2.rdata(realData::yvel) + cpdt1[YDIR] - cpdt2[YDIR];
                    vrel[ZDIR] = p1.rdata(realData::zvel) - p2.rdata(realData::zvel) + cpdt1[ZDIR] - cpdt2[ZDIR];

                    // Calculate normal and tangential components of relative velocity
                    vrel_trans_norm = dotpdt(vrel, normal);
                    vrel_n[XDIR] = vrel_trans_norm*normal[XDIR];
                    vrel_n[YDIR] = vrel_trans_norm*normal[YDIR];
                    vrel_n[ZDIR] = vrel_trans_norm*normal[ZDIR];

                    vrel_t[XDIR] = vrel[XDIR] - vrel_n[XDIR];
                    vrel_t[YDIR] = vrel[YDIR] - vrel_n[YDIR];
                    vrel_t[ZDIR] = vrel[ZDIR] - vrel_n[ZDIR];

                    // Calculate tangential vector here as well for use in liquid bridge model
                    Real mag_vrel_t = sqrt(dotpdt(vrel_t, vrel_t));
                    if(mag_vrel_t > 0.0){
                        tangent[XDIR] = vrel_t[XDIR] / mag_vrel_t;
                        tangent[YDIR] = vrel_t[YDIR] / mag_vrel_t;
                        tangent[ZDIR] = vrel_t[ZDIR] / mag_vrel_t;
                        vrel_trans_tan = dotpdt(vrel, tangent);
                    }
                }
                if ( dist2 <= (rad_sum - TINYVAL)*(rad_sum - TINYVAL) and (p1.id() != p2.id()) )
                {

                    Real kn_dem = DEM::k_n;
                    // TODO: Does eta need to be changed or does it remain the same?
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
                    fn[XDIR] = -(kn_dem*overlap_n*normal[XDIR] + eta_n_dem*vrel_n[XDIR]);
                    fn[YDIR] = -(kn_dem*overlap_n*normal[YDIR] + eta_n_dem*vrel_n[YDIR]);
                    fn[ZDIR] = -(kn_dem*overlap_n*normal[ZDIR] + eta_n_dem*vrel_n[ZDIR]);

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

                    // Total force on composite particle is just a sum across the forces on all component spheres
                    p1.rdata(realData::fx) += fn[XDIR] + ft[XDIR];
                    p1.rdata(realData::fy) += fn[YDIR] + ft[YDIR];
                    p1.rdata(realData::fz) += fn[ZDIR] + ft[ZDIR];

                    //Torque calculation
                    //=================
            
                    // Only particles with a single component sphere should roll
                    Real rollfric_torq[THREEDIM]={0.0};
                    if(p1.idata(intData::num_comp_sphere) == 1){
                        Real dist_c_p1 = 0.5 * (distmag + (p1.rdata(realData::radius)*p1.rdata(realData::radius)
                                        - p2.rdata(realData::radius)*p2.rdata(realData::radius))*distmag_inv);

                        Real omega_p1[THREEDIM];
                        omega_p1[XDIR] = p1.rdata(realData::xangvel);
                        omega_p1[YDIR] = p1.rdata(realData::yangvel);
                        omega_p1[ZDIR] = p1.rdata(realData::zangvel);
                        Real mag_omegap1 = sqrt(dotpdt(omega_p1,omega_p1));
                        if (mag_omegap1 > zero)
                        {
                            Real mag_fn = sqrt(dotpdt(fn,fn));
                            Real torq   = -DEM::muR * mag_fn * dist_c_p1;
                            rollfric_torq[XDIR]=torq*omega_p1[XDIR]/mag_omegap1;
                            rollfric_torq[YDIR]=torq*omega_p1[YDIR]/mag_omegap1;
                            rollfric_torq[ZDIR]=torq*omega_p1[ZDIR]/mag_omegap1;
                        }
                    }
  

                    Real ftot[THREEDIM] = {fn[XDIR]+ft[XDIR], fn[YDIR]+ft[YDIR], fn[ZDIR]+ft[ZDIR]};
                    Real torq[THREEDIM];
                    crosspdt(contact_vect1, ftot, torq);
                    p1.rdata(realData::taux) += torq[XDIR] + rollfric_torq[XDIR];
                    p1.rdata(realData::tauy) += torq[YDIR] + rollfric_torq[YDIR];
                    p1.rdata(realData::tauz) += torq[ZDIR] + rollfric_torq[ZDIR];

                    Real avgtemp=0.5*(p1.rdata(realData::temperature) + p2.rdata(realData::temperature));
                    p1.rdata(realData::temperature)=avgtemp;

                }

                // Calculate the capillary force from the liquid bridge
                if(bridge_exists){
                    // Recalculate bridge volume (not necessary if there's an easy/efficient way to store this) 
                    Real Vp1 = (p1.rdata(realData::liquid_volume)/(p1.idata(intData::num_comp_sphere)*2.0))*(1.0 - sqrt(1.0 - (pow(rad2,two) / (pow(rad1,two) + pow(rad2,two)))));
                    Real Vp2 = (p2.rdata(realData::liquid_volume)/(p2.idata(intData::num_comp_sphere)*2.0))*(1.0 - sqrt(1.0 - (pow(rad1,two) / (pow(rad1,two) + pow(rad2,two)))));
                    Real Vtot = Vp1 + Vp2;

                    // Calculate the separation distance
                    Real Sdist = (distmag - (rad1+rad2) > 0.0) ? distmag - (rad1+rad2):0.0;

                    // Compute rupture distance and remove bridge if necessary
                    Real Drupt = (1.0 + 0.5*DEM::contact_angle)*pow(Vtot,1.0/3.0);
                    if(Sdist > Drupt){
                        p1.idata(intData::first_bridge+bridge_index*3 + 0) = -1;
                        p1.idata(intData::first_bridge+bridge_index*3 + 1) = -1;
                        p1.idata(intData::first_bridge+bridge_index*3 + 2) = -1;
                        p1.rdata(realData::total_bridge_volume) -= Vp1;
                    } else {
                        // Calculate the fitting coefficients from Soulie et al. (2006)
                        Real acoef = -1.1*pow( Vtot/pow(rad2,3.0) ,-0.53);
                        Real bcoef = (-0.148*log(Vtot/pow(rad2,3.0)) - 0.96)*pow(DEM::contact_angle,two) - 0.0082*log(Vtot/pow(rad2,3.0)) + 0.48;
                        Real ccoef = 0.0018*log(Vtot/pow(rad2,3.0)) + 0.078;

                        // Calculate the normal capillary force and normal/tangential viscous forces 
                        Real Fcap = PI*DEM::gamma*sqrt(rad1*rad2)*(ccoef + exp(acoef*Sdist/rad2 + bcoef));
                        
                        Real Fvisc_n = (Sdist > 0.0) ? 6.0*PI*DEM::mu_liq*Reff*vrel_trans_norm*(Reff/Sdist):0.0;
                        Real Fvisc_t = (Sdist > 0.0) ?  6.0*PI*DEM::mu_liq*Reff*vrel_trans_tan*(8.0/15.0 * log(Reff/Sdist) + 0.9588): 0.0;

                        Real bridge_force[THREEDIM];
                        Real torq[THREEDIM];
                        bridge_force[XDIR] = (Fcap + Fvisc_n)*normal[XDIR] + Fvisc_t*tangent[XDIR];
                        bridge_force[YDIR] = (Fcap + Fvisc_n)*normal[YDIR] + Fvisc_t*tangent[YDIR];
                        bridge_force[ZDIR] = (Fcap + Fvisc_n)*normal[ZDIR] + Fvisc_t*tangent[ZDIR];

                        p1.rdata(realData::fx) += bridge_force[XDIR];
                        p1.rdata(realData::fy) += bridge_force[YDIR];
                        p1.rdata(realData::fz) += bridge_force[ZDIR];

                        crosspdt(contact_vect1, bridge_force, torq);
                        p1.rdata(realData::taux) += torq[XDIR];
                        p1.rdata(realData::tauy) += torq[YDIR];
                        p1.rdata(realData::tauz) += torq[ZDIR];
                    }
                }
            }
        }
    }
});
