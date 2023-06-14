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
        if(glued_sphere_particles){
            get_inertial_pos(p1, pc1, p1pos_inert);
        } else {
            p1pos_inert[XDIR] = p1.pos(0);
            p1pos_inert[YDIR] = p1.pos(1);
            p1pos_inert[ZDIR] = p1.pos(2);
        }

        // Calculate the gravitational force and nonzero torque
        p1.rdata(realData::fx) += gravity[XDIR]*p1.rdata(realData::mass)/p1.idata(intData::num_comp_sphere);
        p1.rdata(realData::fy) += gravity[YDIR]*p1.rdata(realData::mass)/p1.idata(intData::num_comp_sphere);
        p1.rdata(realData::fz) += gravity[ZDIR]*p1.rdata(realData::mass)/p1.idata(intData::num_comp_sphere);
  
        Real posvect[THREEDIM];     // Vector from CoM of particle to CoM of component sphere
        posvect[XDIR] = p1pos_inert[XDIR] - p1.pos(0);
        posvect[YDIR] = p1pos_inert[YDIR] - p1.pos(1);
        posvect[ZDIR] = p1pos_inert[ZDIR] - p1.pos(2);

        // TODO: Should code for specific test case be moved into another fcn?
        // Code for glued sphere spin test
        Real f_init[THREEDIM] = {0.0, 0.0, 0.0};
        f_init[init_force_dir] = init_force;
        if(pc1 == init_force_comp && p1.id() == 1){
            p1.rdata(realData::fx) += f_init[XDIR];
            p1.rdata(realData::fy) += f_init[YDIR];
            p1.rdata(realData::fz) += f_init[ZDIR];

            Real torq_init[THREEDIM];
            crosspdt(posvect, f_init, torq_init);
            p1.rdata(realData::taux) += torq_init[XDIR];
            p1.rdata(realData::tauy) += torq_init[YDIR];
            p1.rdata(realData::tauz) += torq_init[ZDIR];
        }

        // Code for cantilever beam test
        // Force is only applied to right-most particle (which generated first with id = 1)
        if((cb_force != 0.0 || cb_torq != 0) && p1.id() == 1){
            Real cb_force_vect[THREEDIM] = {0.0, 0.0, 0.0};
            Real cb_torq_vect[THREEDIM] = {0.0, 0.0, 0.0};
            cb_force_vect[cb_dir] = cb_force;        
            cb_torq_vect[cb_dir] = cb_torq;        
            p1.rdata(realData::fx) += cb_force_vect[XDIR];
            p1.rdata(realData::fy) += cb_force_vect[YDIR];
            p1.rdata(realData::fz) += cb_force_vect[ZDIR];
            p1.rdata(realData::taux) += cb_torq_vect[XDIR];
            p1.rdata(realData::tauy) += cb_torq_vect[YDIR];
            p1.rdata(realData::tauz) += cb_torq_vect[ZDIR];
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
                if(glued_sphere_particles){
                    get_inertial_pos(p2, pc2, p2pos_inert);
                } else {
                    p2pos_inert[XDIR] = p2.pos(0);
                    p2pos_inert[YDIR] = p2.pos(1);
                    p2pos_inert[ZDIR] = p2.pos(2);
                }

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
                        if(p1.idata(intData::first_bridge+3*br+0) == (int)p2.id()){
                            if(p1.idata(intData::first_bridge+3*br+1) == pc1 && p1.idata(intData::first_bridge+3*br+2) == pc2){
                                bridge_exists = true;
                                bridge_index = br;
                            }
                        }
                    }
                }
                
                // Check to see if particles are bonded
                bool bond_exists = false;
                int bond_index = 0;
                if(bonded_sphere_particles){
                    for(int b=0; b<MAXBONDS; b++){
                        if(p1.idata(intData::first_bond+b) == (int)p2.id()){
                            bond_exists = true;
                            bond_index = b;
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

                // Calculate quantities used in contact, liquid bridge, and particle bond force evaluation
                if ( (dist2 <= (rad_sum - TINYVAL)*(rad_sum - TINYVAL) and (p1.id() != p2.id())) || bridge_exists || bond_exists )
                {
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

                    // Add liquid bridge between particles if there isn't one already, particles are not bonded, and at least one particle has moisture content
                    if(liquid_bridging && !bridge_exists && !bond_exists && (p1.rdata(realData::liquid_volume) + p2.rdata(realData::liquid_volume) > 0.0)){
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
                            p1.idata(intData::first_bridge+bridge_index*3 + 0) = (int)p2.id();
                            p1.idata(intData::first_bridge+bridge_index*3 + 1) = pc1;
                            p1.idata(intData::first_bridge+bridge_index*3 + 2) = pc2;
                            bridge_exists = true;
                        }
                    }
                }

                // Evaluate the bond force
                if(bond_exists)
                {
                    // Bond area and moments
                    Real bond_rad = std::min(rad1, rad2) * DEM::bond_radius_factor;
                    Real A_bond = PI * pow(bond_rad,two);
                    Real I_bond = PI * pow(bond_rad,4) / 4.0;
                    Real I_p_bond = I_bond * 2.0;

                    // Stiffnesses for Bernoulli beam model and bond damping calculations
                    Real k_n_bond = DEM::E_bond * A_bond / rad_sum;
                    Real k_t_bond = DEM::G_bond * A_bond / rad_sum;
                    Real k_tor_bond = DEM::G_bond * I_p_bond / rad_sum;
                    Real k_ben_bond = DEM::E_bond * I_bond / rad_sum;

                    /* ----- COMMENTING OUT OLD TBBM CODE, NOT WORKING ------ */
                    // // TBBM Coefficients
                    // Real phi_bond = 20.0*bond_rad*bond_rad*(1.0+DEM::nu_bond) / (3.0*rad_sum*rad_sum);
                    // Real k1 = DEM::E_bond * A_bond / rad_sum;
                    // Real k2 = 12.0*DEM::E_bond*I_bond / (pow(rad_sum,3.0)*(1.0+phi_bond));
                    // Real k3 = 0.5*rad_sum*k2;
                    // Real k4 = DEM::E_bond*I_bond / (rad_sum * (1.0 + DEM::nu_bond));
                    // Real k5 = DEM::E_bond*I_bond / rad_sum * ((4.0+phi_bond)/(1.0+phi_bond));
                    // Real k6 = DEM::E_bond*I_bond / rad_sum * ((2.0-phi_bond)/(1.0+phi_bond));

                    // // Calculate orthonormal basis vectors
                    // Real sdotn;
                    // Real checkpar[THREEDIM];
                    // Real temp1[THREEDIM];
                    // crosspdt(normal, s_vect, checkpar);
                    // if(sqrt(dotpdt(checkpar,checkpar)) > 1.0e-5){
                    //     sdotn = dotpdt(normal, s_vect);
                    //     temp1[XDIR] = s_vect[XDIR] - sdotn*normal[XDIR]; temp1[YDIR] = s_vect[YDIR] - sdotn*normal[YDIR]; temp1[ZDIR] = s_vect[ZDIR] - sdotn*normal[ZDIR];
                    // } else {
                    //     sdotn = dotpdt(normal, s_vect2);
                    //     temp1[XDIR] = s_vect2[XDIR] - sdotn*normal[XDIR]; temp1[YDIR] = s_vect2[YDIR] - sdotn*normal[YDIR]; temp1[ZDIR] = s_vect2[ZDIR] - sdotn*normal[ZDIR];
                    // }
                    // Real temp_mag = sqrt(dotpdt(temp1, temp1));
                    // Real basis1[THREEDIM] = {temp1[XDIR]/temp_mag, temp1[YDIR]/temp_mag, temp1[ZDIR]/temp_mag};
                    // Real temp2[THREEDIM];
                    // crosspdt(normal, basis1, temp2);
                    // Real temp_mag2 = sqrt(dotpdt(temp2, temp2));
                    // Real basis2[THREEDIM] = {temp2[XDIR]/temp_mag2, temp2[YDIR]/temp_mag2, temp2[ZDIR]/temp_mag2};   // Just in case?

                    // // Calculate p1 and p2 incremental displacements
                    // Real inertial_trans_disp1[THREEDIM] = {p1.rdata(realData::xvel)*dt, p1.rdata(realData::yvel)*dt, p1.rdata(realData::zvel)*dt};
                    // Real inertial_trans_disp2[THREEDIM] = {p2.rdata(realData::xvel)*dt, p2.rdata(realData::yvel)*dt, p2.rdata(realData::zvel)*dt};
                    // Real inertial_ang_disp1[THREEDIM] = {p1.rdata(realData::xangvel)*dt, p1.rdata(realData::yangvel)*dt, p1.rdata(realData::zangvel)*dt};
                    // Real inertial_ang_disp2[THREEDIM] = {p2.rdata(realData::xangvel)*dt, p2.rdata(realData::yangvel)*dt, p2.rdata(realData::zangvel)*dt};

                    // // Calculate local bond frame displacements
                    // Real Ux1 = dotpdt(inertial_trans_disp1, normal);
                    // Real Uy1 = dotpdt(inertial_trans_disp1, basis1);
                    // Real Uz1 = dotpdt(inertial_trans_disp1, basis2);
                    // Real tx1 = dotpdt(inertial_ang_disp1, normal);
                    // Real ty1 = dotpdt(inertial_ang_disp1, basis1);
                    // Real tz1 = dotpdt(inertial_ang_disp1, basis2);
                    // Real Ux2 = dotpdt(inertial_trans_disp2, normal);
                    // Real Uy2 = dotpdt(inertial_trans_disp2, basis1);
                    // Real Uz2 = dotpdt(inertial_trans_disp2, basis2);
                    // Real tx2 = dotpdt(inertial_ang_disp2, normal);
                    // Real ty2 = dotpdt(inertial_ang_disp2, basis1);
                    // Real tz2 = dotpdt(inertial_ang_disp2, basis2);

                    // // Calculate local bond frame incremental forces/torques
                    // Real dFx = k1*Ux1 - k1*Ux2;
                    // Real dFy = k2*Uy1 + k3*tz1 - k2*Uy2 + k3*tz2;
                    // Real dFz = k2*Uz1 - k3*ty1 - k2*Uz2 - k3*ty2;
                    // Real dMx = k4*tx1 - k4*tx2;
                    // Real dMy = -k3*Uz1 + k5*ty1 + k3*Uz2 + k6*ty2;
                    // Real dMz = k3*Uy1 + k5*tz1 - k3*Uy2 + k6*tz2;

                    // // Calculate normal and tangential relative angular velocity components
                    // Real angvel_rel[THREEDIM];
                    // angvel_rel[XDIR] = angvel1[XDIR] - angvel2[XDIR];  
                    // angvel_rel[YDIR] = angvel1[YDIR] - angvel2[YDIR];  
                    // angvel_rel[ZDIR] = angvel1[ZDIR] - angvel2[ZDIR];  

                    // Real angvel_rel_norm;
                    // Real angvel_rel_n[THREEDIM];
                    // Real angvel_rel_t[THREEDIM];

                    // angvel_rel_norm = dotpdt(angvel_rel, normal);
                    // angvel_rel_n[XDIR] = angvel_rel_norm*normal[XDIR];
                    // angvel_rel_n[YDIR] = angvel_rel_norm*normal[YDIR];
                    // angvel_rel_n[ZDIR] = angvel_rel_norm*normal[ZDIR];

                    // angvel_rel_t[XDIR] = angvel_rel[XDIR] - angvel_rel_n[XDIR];
                    // angvel_rel_t[YDIR] = angvel_rel[YDIR] - angvel_rel_n[YDIR];
                    // angvel_rel_t[ZDIR] = angvel_rel[ZDIR] - angvel_rel_n[ZDIR];

                    // // Calculate the damping forces
                    // Real F_dn[THREEDIM];
                    // Real F_dt[THREEDIM];
                    // Real M_dn[THREEDIM];
                    // Real M_dt[THREEDIM];

                    // F_dn[XDIR] = DEM::beta_bond * pow(2.0*m1*k_n_bond,0.5) * vrel_n[XDIR];
                    // F_dn[YDIR] = DEM::beta_bond * pow(2.0*m1*k_n_bond,0.5) * vrel_n[YDIR];
                    // F_dn[ZDIR] = DEM::beta_bond * pow(2.0*m1*k_n_bond,0.5) * vrel_n[ZDIR];
                    // F_dt[XDIR] = DEM::beta_bond * pow(2.0*m1*k_t_bond,0.5) * vrel_t[XDIR];
                    // F_dt[YDIR] = DEM::beta_bond * pow(2.0*m1*k_t_bond,0.5) * vrel_t[YDIR];
                    // F_dt[ZDIR] = DEM::beta_bond * pow(2.0*m1*k_t_bond,0.5) * vrel_t[ZDIR];
                    // 
                    // M_dn[XDIR] = DEM::beta_bond * pow(2.0*k_tor_bond/p1.rdata(realData::Iinv),0.5) * angvel_rel_n[XDIR];
                    // M_dn[YDIR] = DEM::beta_bond * pow(2.0*k_tor_bond/p1.rdata(realData::Iinv),0.5) * angvel_rel_n[YDIR];
                    // M_dn[ZDIR] = DEM::beta_bond * pow(2.0*k_tor_bond/p1.rdata(realData::Iinv),0.5) * angvel_rel_n[ZDIR];
                    // M_dt[XDIR] = DEM::beta_bond * pow(2.0*k_ben_bond/p1.rdata(realData::Iinv),0.5) * angvel_rel_t[XDIR];
                    // M_dt[YDIR] = DEM::beta_bond * pow(2.0*k_ben_bond/p1.rdata(realData::Iinv),0.5) * angvel_rel_t[YDIR];
                    // M_dt[ZDIR] = DEM::beta_bond * pow(2.0*k_ben_bond/p1.rdata(realData::Iinv),0.5) * angvel_rel_t[ZDIR];

                    Real radsuminv = 1.0 / (rad1 + rad2);
                    Real wr[THREEDIM] = {p1.rdata(realData::xangvel) - p2.rdata(realData::xangvel), 
                                         p1.rdata(realData::yangvel) - p2.rdata(realData::yangvel), 
                                         p1.rdata(realData::zangvel) - p2.rdata(realData::zangvel)}; 

                    Real wdp = dotpdt(wr, normal);
                    Real wn[THREEDIM] = {wdp*normal[XDIR], wdp*normal[YDIR], wdp*normal[ZDIR]};
                    Real wt[THREEDIM] = {wr[XDIR] - wn[XDIR], wr[YDIR] - wn[YDIR], wr[ZDIR] - wn[ZDIR]};

                    /* ------- Normal bond force ------- */
                    Real nforce[THREEDIM] = {0.0, 0.0, 0.0};
                    Real nforce_damped[THREEDIM] = {0.0, 0.0, 0.0};
                    nforce[XDIR] = -k_n_bond*overlap_n*normal[XDIR];
                    nforce[YDIR] = -k_n_bond*overlap_n*normal[YDIR];
                    nforce[ZDIR] = -k_n_bond*overlap_n*normal[ZDIR];

                    Real minvel = 1e-5 * std::min(rad1,rad2) /dt;
                    Real dvX = 0.01*nforce[XDIR]*dt;
                    Real dvY = 0.01*nforce[YDIR]*dt;
                    Real dvZ = 0.01*nforce[ZDIR]*dt;
                    Real multiplierX = amrex::Math::copysign(1.0,vrel_n[0]);
                    Real multiplierY = amrex::Math::copysign(1.0,vrel_n[1]);
                    Real multiplierZ = amrex::Math::copysign(1.0,vrel_n[2]);
                    nforce_damped[0] = nforce[0] - DEM::beta_bond*fabs(nforce[0])*multiplierX;
                    nforce_damped[1] = nforce[1] - DEM::beta_bond*fabs(nforce[1])*multiplierY;
                    nforce_damped[2] = nforce[2] - DEM::beta_bond*fabs(nforce[2])*multiplierZ;

                    /* ------- Tangent bond force ------- */
                    Real force_tang[THREEDIM] = {p1.rdata(realData::first_bond_v+bond_index*9+0), 
                                                 p1.rdata(realData::first_bond_v+bond_index*9+1), 
                                                 p1.rdata(realData::first_bond_v+bond_index*9+2)};
                    Real tforce_damped[THREEDIM] = {0.0, 0.0, 0.0};
                    Real fdp = dotpdt(force_tang, normal);
                    Real force_norm[THREEDIM] = {fdp*normal[XDIR], fdp*normal[YDIR], fdp*normal[ZDIR]};
                    force_tang[XDIR] -= force_norm[XDIR];
                    force_tang[YDIR] -= force_norm[YDIR];
                    force_tang[ZDIR] -= force_norm[ZDIR];

                    Real dtforce[THREEDIM] = {-vrel_t[XDIR]*k_t_bond*dt, -vrel_t[YDIR]*k_t_bond*dt, -vrel_t[ZDIR]*k_t_bond*dt};

                    force_tang[XDIR] += dtforce[XDIR];
                    force_tang[YDIR] += dtforce[YDIR];
                    force_tang[ZDIR] += dtforce[ZDIR];

                    dvX = 0.01*force_tang[XDIR]*dt;
                    dvY = 0.01*force_tang[YDIR]*dt;
                    dvZ = 0.01*force_tang[ZDIR]*dt;
                    multiplierX = amrex::Math::copysign(1.0,vrel_t[0]);
                    multiplierY = amrex::Math::copysign(1.0,vrel_t[1]);
                    multiplierZ = amrex::Math::copysign(1.0,vrel_t[2]);
                    tforce_damped[0] = force_tang[0] - DEM::beta_bond*fabs(force_tang[0])*multiplierX;
                    tforce_damped[1] = force_tang[1] - DEM::beta_bond*fabs(force_tang[1])*multiplierY;
                    tforce_damped[2] = force_tang[2] - DEM::beta_bond*fabs(force_tang[2])*multiplierZ;

                    p1.rdata(realData::first_bond_v+bond_index*9+0) = tforce_damped[0];
                    p1.rdata(realData::first_bond_v+bond_index*9+1) = tforce_damped[1];
                    p1.rdata(realData::first_bond_v+bond_index*9+2) = tforce_damped[2];

                    /* ------- Normal bond torque ------- */
                    Real torq_norm[THREEDIM] = {p1.rdata(realData::first_bond_v+bond_index*9+3), 
                                                p1.rdata(realData::first_bond_v+bond_index*9+4), 
                                                p1.rdata(realData::first_bond_v+bond_index*9+5)};
                    Real ntorq_damped[THREEDIM] = {0.0, 0.0, 0.0};
                    Real tdp = dotpdt(torq_norm, normal);
                    torq_norm[XDIR] = tdp*normal[XDIR];
                    torq_norm[YDIR] = tdp*normal[YDIR];
                    torq_norm[ZDIR] = tdp*normal[ZDIR];

                    Real dntorque[THREEDIM] = {-wn[XDIR]*k_tor_bond*dt, -wn[YDIR]*k_tor_bond*dt, -wn[ZDIR]*k_tor_bond*dt};
                    torq_norm[XDIR] += dntorque[XDIR];
                    torq_norm[YDIR] += dntorque[YDIR];
                    torq_norm[ZDIR] += dntorque[ZDIR];
                    ntorq_damped[0] = torq_norm[0] - DEM::beta_bond*fabs(torq_norm[0])*amrex::Math::copysign(1.0,wn[0]);
                    ntorq_damped[1] = torq_norm[1] - DEM::beta_bond*fabs(torq_norm[1])*amrex::Math::copysign(1.0,wn[1]);
                    ntorq_damped[2] = torq_norm[2] - DEM::beta_bond*fabs(torq_norm[2])*amrex::Math::copysign(1.0,wn[2]);

                    p1.rdata(realData::first_bond_v+bond_index*9+3) = ntorq_damped[0];
                    p1.rdata(realData::first_bond_v+bond_index*9+4) = ntorq_damped[1];
                    p1.rdata(realData::first_bond_v+bond_index*9+5) = ntorq_damped[2];

                    /* ------- Tangential bond torque ------- */
                    Real torq_tang[THREEDIM] = {p1.rdata(realData::first_bond_v+bond_index*9+6), 
                                                p1.rdata(realData::first_bond_v+bond_index*9+7), 
                                                p1.rdata(realData::first_bond_v+bond_index*9+8)};
                    Real ttorq_damped[THREEDIM] = {0.0, 0.0, 0.0};
                    tdp = dotpdt(torq_tang, normal);
                    torq_tang[XDIR] -= tdp*normal[XDIR];
                    torq_tang[YDIR] -= tdp*normal[YDIR];
                    torq_tang[ZDIR] -= tdp*normal[ZDIR];

                    Real dttorque[THREEDIM] = {-wt[XDIR]*k_ben_bond*dt, -wt[YDIR]*k_ben_bond*dt, -wt[ZDIR]*k_ben_bond*dt};
                    torq_tang[XDIR] += dttorque[XDIR];
                    torq_tang[YDIR] += dttorque[YDIR];
                    torq_tang[ZDIR] += dttorque[ZDIR];
                    ttorq_damped[0] = torq_tang[0] - DEM::beta_bond*fabs(torq_tang[0])*amrex::Math::copysign(1.0,wt[0]);
                    ttorq_damped[1] = torq_tang[1] - DEM::beta_bond*fabs(torq_tang[1])*amrex::Math::copysign(1.0,wt[1]);
                    ttorq_damped[2] = torq_tang[2] - DEM::beta_bond*fabs(torq_tang[2])*amrex::Math::copysign(1.0,wt[2]);

                    p1.rdata(realData::first_bond_v+bond_index*9+6) = ttorq_damped[0];
                    p1.rdata(realData::first_bond_v+bond_index*9+7) = ttorq_damped[1];
                    p1.rdata(realData::first_bond_v+bond_index*9+8) = ttorq_damped[2];

                    // Add on cross product term (should not be included in tau_bond_t update)
                    Real tor[THREEDIM] = {0.0, 0.0, 0.0};
                    crosspdt(tforce_damped, normal, tor); 
                    Real cri = distmag * rad1 * radsuminv;

                    // Total bond forces and torques are added to particle quantities (check sign for crosspdt term?)
                    p1.rdata(realData::fx) += nforce_damped[XDIR] + tforce_damped[XDIR];
                    p1.rdata(realData::fy) += nforce_damped[YDIR] + tforce_damped[YDIR];
                    p1.rdata(realData::fz) += nforce_damped[ZDIR] + tforce_damped[ZDIR];

                    p1.rdata(realData::taux) += ntorq_damped[XDIR] + ttorq_damped[XDIR] - cri*tor[XDIR];
                    p1.rdata(realData::tauy) += ntorq_damped[YDIR] + ttorq_damped[YDIR] - cri*tor[YDIR];
                    p1.rdata(realData::tauz) += ntorq_damped[ZDIR] + ttorq_damped[ZDIR] - cri*tor[ZDIR];
                }

                // Evaluate the contact force if p1 and p2 don't have a bond between them
                // if ( dist2 <= (rad_sum - TINYVAL)*(rad_sum - TINYVAL) and (p1.id() != p2.id()))
                if ( dist2 <= (rad_sum - TINYVAL)*(rad_sum - TINYVAL) and (p1.id() != p2.id()) && !bond_exists)
                {
                    Real kn_dem, kt_dem, eta_n_dem, eta_t_dem;
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
                        eta_n_dem = -2.0*sqrt(5.0/3.0)*beta*sqrt(meff*E_eff)*sqrt(temp1);
                        eta_t_dem = -4.0*sqrt(5.0/3.0)*beta*sqrt(meff*G_eff)*sqrt(temp1);
                    }
                    else
                    {            
                        kn_dem = DEM::k_n;
                        // TODO: Does eta need to be changed or does it remain the same?
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
                        // Real Fcap = PI*DEM::gamma*sqrt(rad1*rad2)*(ccoef + exp(acoef*Sdist/rad2 + bcoef));
                        // Real Fcap = (Sdist > 0.0) ? PI*DEM::gamma*sqrt(rad1*rad2)*(ccoef + exp(acoef*Sdist/rad2 + bcoef)):0.0;
                        Real Fcap = (Sdist > 0.0) ? PI*DEM::gamma*sqrt(rad1*rad2)*(ccoef + exp(acoef*Sdist/rad2 + bcoef)) * 0.5*(tanh(5000.0*(Sdist - Drupt/10.0)) + 1.0):0.0;
                        
                        Real Fvisc_n = (Sdist > 0.0) ? 6.0*PI*DEM::mu_liq*Reff*vrel_trans_norm*(Reff/Sdist):0.0;
                        Real Fvisc_t = (Sdist > 0.0) ?  6.0*PI*DEM::mu_liq*Reff*vrel_trans_tan*(8.0/15.0 * log(Reff/Sdist) + 0.9588): 0.0;

                        Real bridge_force[THREEDIM];
                        Real torq[THREEDIM];
                        bridge_force[XDIR] = (Fcap - Fvisc_n)*normal[XDIR] - Fvisc_t*tangent[XDIR];
                        bridge_force[YDIR] = (Fcap - Fvisc_n)*normal[YDIR] - Fvisc_t*tangent[YDIR];
                        bridge_force[ZDIR] = (Fcap - Fvisc_n)*normal[ZDIR] - Fvisc_t*tangent[ZDIR];

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
