#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <BDEM_ParticleContainer.H>
#include <BDEM_EB.H>
#include <STLtools.H>

//definition
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::k_n      = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::e_n      = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::mu       = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::k_n_wall = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::e_n_wall = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::E_wall  = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::nu_wall  = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::mu_wall  = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::k_t      = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::e_t      = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::muR       = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::k_t_wall = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::e_t_wall = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::muR_wall  = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::tcoll    = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::mu_liq   = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::hminf   = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::CED   = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::contact_angle = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::gamma    = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::wall_gamma    = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::E_bond  = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::G_bond  = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::beta_bond  = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::global_damping  = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::force_damping  = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::angv_damping  = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::nu_bond  = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::bond_radius_factor  = one;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::k_c  = one;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::sigma_max  = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::eps_g  = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::rho_g  = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::mu_g  = zero;
AMREX_GPU_DEVICE_MANAGED int DEM::particles_in_parcel = 1;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::CG_ratio  = 1.;

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    auto start = std::chrono::high_resolution_clock::now();
    
    DEMspecs specs;
    specs.read_dem_specs();
    RealBox real_box;
    for (int n = 0; n < 3; n++)
    {
        real_box.setLo(n, specs.plo[n]);
        real_box.setHi(n, specs.phi[n]);
    }

    IntVect domain_lo(AMREX_D_DECL(0,0,0));
    IntVect domain_hi(AMREX_D_DECL(specs.ncells[XDIR]-1,
                specs.ncells[YDIR]-1,
                specs.ncells[ZDIR]-1));

    const Box domain(domain_lo, domain_hi);

    int coord = 0; //cartesian
    Geometry geom(domain, &real_box, coord, specs.periodic.data());

    BoxArray ba(domain);
    ba.maxSize(specs.max_grid_size);
    DistributionMapping dm(ba);

    EBtools::init_eb(geom,ba,dm); 

    const int ng_cells = specs.ng_cells;
    BDEMParticleContainer bpc(geom, dm, ba, ng_cells,specs.chemptr);
    if(!specs.restartedcase)
    {
        if(specs.init_particles_using_file == 1)
        {
            if(specs.bonded_sphere_particles){
                bpc.InitBondedParticles("particle_input.dat", specs.do_heat_transfer, specs.contact_law, 
                                        specs.cantilever_beam_test, specs.liquid_bridging,
                                        specs.liquid_density, specs.moisture_content, 
                                        specs.moisture_content_stdev,specs.FSP);
            } else {
                bpc.InitParticles("particle_input.dat",specs.do_heat_transfer, 
                                    specs.temp_mean, specs.temp_stdev, specs.contact_law,
                                    specs.liquid_bridging, specs.liquid_density, 
                                    specs.moisture_content, specs.moisture_content_stdev, specs.FSP);
            }
            bpc.InitChemSpecies(specs.species_massfracs.data());
        }
        else
        {
            Print()<<"Doing autogeneration\n";
            bpc.InitParticles (specs.autogen_mincoords.data(),specs.autogen_maxcoords.data(), 
                    specs.autogen_meanvel.data(),  specs.autogen_fluctuation.data(), 
                    specs.autogen_radius, specs.autogen_dens, specs.autogen_E, specs.autogen_nu, 
                    specs.autogen_temp,
                    specs.autogen_species_massfracs.data(),
                    specs.autogen_multi_part_per_cell, specs.autogen_layer_particles, 
                    specs.bonded_sphere_particles,
                    specs.autogen_min_radius, specs.autogen_max_radius,
                    specs.autogen_bp_type, specs.particle_type_list, 
                    specs.use_type_dist, specs.particle_dist_list,
                    specs.liquid_bridging,
                    specs.liquid_density, specs.moisture_content,
                    specs.moisture_content_stdev, specs.FSP);
        }
    }
    else
    {
        bpc.Restart(specs.restartfilename, "particles");
        if(specs.reassign_particle_properties) bpc.reassignParticleProperties(specs.autogen_radius, specs.autogen_dens, specs.autogen_E, specs.autogen_nu);
    }

    if(specs.coloring_striation) bpc.particleColoringStriation(specs.coloring_striation_len, specs.coloring_striation_dir);

    if(specs.modify_species)
    {
        bpc.InitChemSpecies(specs.modified_species_massfracs.data());
    }
    if(specs.stratified_specs)
    {
        bpc.InitChemSpecies(specs.nstrat_domains,specs.strat_mincoords.data(),
                            specs.strat_maxcoords.data(),specs.strat_spec_massfracs.data());
    }
        
    //compute tcoll here over all particles
    DEM::tcoll=bpc.compute_coll_timescale(specs.bonded_sphere_particles);
    Print()<<"tcoll:"<<DEM::tcoll<<"\n";

    Real dt = (specs.constant_dt > 0.0) ? specs.constant_dt:specs.cfl*DEM::tcoll;
    int steps=0;
    Real time=zero;
    Real output_time=zero;
    Real particle_sourcing_time=zero;
    Real output_timeMass=zero;
    Real output_timePrint=zero;
    int output_it=0;
    amrex::Print() << "Time step dt = " << dt << "\n";

    amrex::Print() << "Num particles before eb removal  " << bpc.TotalNumberOfParticles() << "\n";
    //if(EBtools::using_levelset_geometry and !specs.restartedcase)
    if(!specs.restartedcase){
        if(EBtools::using_levelset_geometry)
        {
            bpc.removeParticlesOutsideBoundary(EBtools::lsphi,
                                                EBtools::ebfactory,EBtools::ls_refinement);
            amrex::Print() << "Num particles after eb removal  " << bpc.TotalNumberOfParticles() << "\n";
        }

        for (int stli = 0; stli < specs.stls.size(); stli++)
        {
            if(specs.stls[stli].remove_particles_inside == 1)
            {
                bpc.removeParticlesInsideSTL(specs.outside_point, specs.stls[stli].stlptr);
                amrex::Print() << "Num particles after stl removal " << bpc.TotalNumberOfParticles() << "\n";
            }
        }

        if(specs.remove_eb_overlaps) {
            bpc.removeEBOverlapParticles(EBtools::ebfactory, EBtools::lsphi, EBtools::ls_refinement);
            amrex::Print() << "Num particles after EB overlap removal " << bpc.TotalNumberOfParticles() << "\n";
        }
    }
    if(specs.clip_particles) {
        bpc.clipParticles(specs.clip_particle_dir, specs.clip_particle_val);
        amrex::Print() << "Num particles after clipping " << bpc.TotalNumberOfParticles() << "\n";
    }

    bpc.set_domain_bcs(specs.bclo,specs.bchi);
    bpc.writeParticles(steps+specs.stepoffset, specs.bonded_sphere_particles, specs.pltprefix);

    for (int stli = 0; stli < specs.stls.size(); stli++)
    {
        std::string stlpltfile = amrex::Concatenate( (specs.stls[stli].name + "_").c_str(), steps+specs.stepoffset, 5);
        if(specs.stls[stli].dynamicstl!=0)
        {
            specs.stls[stli].stlptr->write_stl_file(stlpltfile);                
        }
        specs.stls[stli].stlptr->writeVTK(stlpltfile);
    }  

    amrex::Print() << "Num particles after init is " << bpc.TotalNumberOfParticles() << "\n";

    // Calculate the moisture content for each particle
    if(specs.liquid_bridging && specs.recalculate_MC) bpc.computeMoistureContent(specs.moisture_content, specs.moisture_content_stdev, specs.liquid_density, specs.FSP);

    amrex::Real mass_flow_prev;
    amrex::Real mass_flow_next;
    amrex::Real total_speed=1e10;
    amrex::Real t_prev=0.0;
    if(specs.calc_mass_flow) bpc.Calculate_Total_Mass_MaterialPoints(mass_flow_prev, specs.mass_flow_dir, specs.mass_flow_cutoff);

    Real cb_force = 0.0;
    Real cb_torq = 0.0;
    Real cb_intt = specs.cb_load_time + specs.cb_hold_time;
    Real cb_int_stop = cb_intt*specs.cb_interval_num;

    while((steps < specs.maxsteps) and (time < specs.final_time) and (total_speed > specs.avg_speed_stop))
    {
        time += dt;
        output_time += dt;
        output_timeMass += dt;
        output_timePrint += dt;
        particle_sourcing_time += dt;
    
        if(steps>0) specs.init_force = 0.0;
        if(specs.cantilever_beam_test){ 
            if(specs.cb_intervals){
                Real temp_t = (time > cb_intt) ? remainder(time, cb_intt):time;
                if(temp_t > 0 && temp_t <= specs.cb_load_time && time < cb_int_stop) cb_force += dt / (specs.cb_interval_num*specs.cb_load_time) * specs.cb_force_max;
                if(temp_t > 0 && temp_t <= specs.cb_load_time && time < cb_int_stop) cb_torq += dt / (specs.cb_interval_num*specs.cb_load_time) * specs.cb_torq_max;
            } else {
                cb_force = (time < specs.cb_load_time) ? (time/specs.cb_load_time)*specs.cb_force_max:
                                                            (time < specs.cb_unload_time) ? specs.cb_force_max:0.0;
                cb_torq = (time < specs.cb_load_time) ?  (time/specs.cb_load_time)*specs.cb_torq_max:
                                                            (time < specs.cb_unload_time) ? specs.cb_torq_max:0.0;
            }
        }

        if(specs.particle_sourcing==1 && 
            particle_sourcing_time > specs.particle_sourcing_interval
            && time < specs.stop_sourcing_time)

        {
            amrex::Print() << "Doing particle sourcing\n";
            bpc.clearSourcingVolume(specs.particle_sourcing_mincoords.data(), specs.particle_sourcing_maxcoords.data());
            bpc.InitParticles (specs.particle_sourcing_mincoords.data(),specs.particle_sourcing_maxcoords.data(), 
                                specs.particle_sourcing_meanvel.data(),  specs.particle_sourcing_fluctuation.data(), 
                                specs.particle_sourcing_radius, specs.particle_sourcing_dens,
                                specs.particle_sourcing_E, specs.particle_sourcing_nu,
                                specs.particle_sourcing_temp,
                                specs.particle_sourcing_species_massfracs.data(),
                                specs.particle_sourcing_multi_part_per_cell, 
                                specs.particle_sourcing_layer_particles,
                                specs.bonded_sphere_particles,
                                specs.particle_sourcing_min_radius, specs.particle_sourcing_max_radius,
                                specs.particle_sourcing_bp_type, specs.particle_type_list, 
                                specs.use_type_dist, specs.particle_dist_list,
                                specs.liquid_bridging, specs.liquid_density,
                                specs.moisture_content, specs.moisture_content_stdev, specs.FSP);

            if(EBtools::using_levelset_geometry)
            {
                bpc.removeParticlesOutsideBoundary(EBtools::lsphi,
                                                    EBtools::ebfactory,EBtools::ls_refinement);
            }
            amrex::Print() << "Num particles after eb removal  " << bpc.TotalNumberOfParticles() << "\n";
            for (int stli = 0; stli < specs.stls.size(); stli++)
            {
                if(specs.stls[stli].remove_particles_inside == 1)
                {
                    bpc.removeParticlesInsideSTL(specs.outside_point, specs.stls[stli].stlptr);
                    amrex::Print() << "Num particles after stl removal " << bpc.TotalNumberOfParticles() << "\n";
                }
            }

            


            bpc.redist_particles(1,specs.using_softwalls);
            amrex::Print()<<"adding particles\n";
            if(specs.using_softwalls)
            {
                bpc.reassignParticles_softwall(); 
            }
            amrex::Print() << "Num particles after sourcing " << bpc.TotalNumberOfParticles() << "\n";
            particle_sourcing_time=zero;
        }

        if (steps % specs.num_redist == 0)
        {
            if(specs.remove_eb_overlaps) {
                bpc.removeEBOverlapParticles(EBtools::ebfactory, EBtools::lsphi, EBtools::ls_refinement);
                amrex::Print() << "Num particles after EB overlap removal " << bpc.TotalNumberOfParticles() << "\n";
            }

            bpc.redist_particles(1,specs.using_softwalls);
            if(specs.using_softwalls)
            {
                bpc.reassignParticles_softwall(); 
            }
        } 
        else if (steps % specs.num_upd_neighbor == 0)
        {
            bpc.updateNeighbors();
        }


        /*if(specs.stl_geom_present)
            {
            bpc.checkParticlesInsideSTL(specs.outside_point);
            }*/

        if(specs.verlet_scheme){
            BL_PROFILE_VAR("MOVE_PART",movepart);
            specs.verlet_scheme = 1;
            bpc.moveParticles(dt,specs.do_chemistry,specs.minradius_frac,specs.verlet_scheme);
            BL_PROFILE_VAR_STOP(movepart);
        }

        BL_PROFILE_VAR("FORCE_CALC",forceCalc);
        {

            bpc.computeForces(dt,EBtools::ebfactory,EBtools::lsphi,
                                specs.do_heat_transfer,specs.walltemp_vardir,
                                specs.walltemp_polynomial.data(),
                                EBtools::ls_refinement,specs.stl_geom_present, specs.contact_law, steps,
                                specs.gravity, specs.stls, time,
                                specs.bonded_sphere_particles,
                                specs.liquid_bridging, 
                                specs.particle_cohesion,
                                specs.init_force, specs.init_force_dir, specs.init_force_comp,
                                cb_force, cb_torq, specs.cb_dir, specs.drag_model,
                                specs.solve_fibrillation);
        }
        BL_PROFILE_VAR_STOP(forceCalc);
        

        
        BL_PROFILE_VAR("MOVE_PART",movepart);
        if(specs.verlet_scheme) specs.verlet_scheme = 2;
        bpc.moveParticles(dt,specs.do_chemistry,specs.minradius_frac,specs.verlet_scheme);
        BL_PROFILE_VAR_STOP(movepart);

        for (int stli = 0; stli < specs.stls.size(); stli++)
        {

            if(specs.stls[stli].dynamicstl!=0)
            {
                if(specs.stls[stli].dynamicstl==1)
                {
                    specs.stls[stli].stlptr->move_stl(dt,specs.stls[stli].dynamicstl,specs.stls[stli].dynstl_transl_dir.data(),
                                        specs.stls[stli].dynstl_center.data(),specs.stls[stli].dynstl_transl_vel);
                }
                else if(specs.stls[stli].dynamicstl==2)
                {
                    specs.stls[stli].stlptr->move_stl(dt,specs.stls[stli].dynamicstl,specs.stls[stli].dynstl_rot_dir.data(),
                                        specs.stls[stli].dynstl_center.data(),specs.stls[stli].dynstl_rot_vel);

                }
                else if(specs.stls[stli].dynamicstl==3)
                {
                    amrex::Real stl_pos_prv = specs.stls[stli].stl_vib_amp*sin((time+specs.stls[stli].stl_timeoffset)*specs.stls[stli].stl_vib_freq);
                    amrex::Real stl_pos_nxt = specs.stls[stli].stl_vib_amp*sin((time+dt+specs.stls[stli].stl_timeoffset)*specs.stls[stli].stl_vib_freq);
                    amrex::Real stl_vel = (stl_pos_nxt - stl_pos_prv)/dt;
                    specs.stls[stli].stlptr->move_stl(dt,specs.stls[stli].dynamicstl,specs.stls[stli].stl_vib_dir.data(),specs.stls[stli].dynstl_center.data(),stl_vel);
                }
                else
                {
                    amrex::Abort("STL motion not implemented\n");
                }
                if (steps % specs.num_redist == 0)
                {
                    specs.stls[stli].stlptr->boxmap.clear();
                    specs.stls[stli].stlptr->orb_of_triangulation(0,specs.stls[stli].stlptr->num_tri-1,0,specs.stls[stli].stlptr->sorted_indexarray);
                }
            }

        }
        


        if (output_timePrint > specs.screen_output_time)
        {


            Print() << "step:" << steps << "\t" 
                    << "time:" << time << "\t" 
                    << "clock time:" << specs.elapsed_clock_time() << " s\n";
            
            output_timePrint=zero;
        }

        if (specs.mass_flow)
        {
            if (output_timeMass > specs.massflow_output_time)
            {
                bpc.redist_particles(0,specs.using_softwalls);
                if(specs.using_softwalls)
                {
                    bpc.reassignParticles_softwall(); 
                }
                PrintToFile("Total_Mass") << time << "\t" << bpc.TotalNumberOfParticles() <<"\n";
                output_timeMass=zero;
            }
        }

        if (output_time > specs.write_output_time) 
        {
            BL_PROFILE_VAR("OUTPUT_TIME",outputs);
            Print()<<"writing outputs at step,time:"<<steps<<"\t"<<time<<"\n";
            bpc.redist_particles(0,specs.using_softwalls);
            output_it++;
            bpc.writeParticles(output_it+specs.stepoffset, specs.bonded_sphere_particles, specs.pltprefix);
            const std::string& rstfile = specs.pltprefix + amrex::Concatenate("rst", output_it+specs.stepoffset, 5);
            bpc.Checkpoint(rstfile, "particles");
            output_time=zero;

            for (int stli = 0; stli < specs.stls.size(); stli++)
            {
                
                std::string stlpltfile = amrex::Concatenate( (specs.stls[stli].name + "_").c_str(), output_it+specs.stepoffset, 5);
                if(specs.stls[stli].dynamicstl!=0)
                {
                    specs.stls[stli].stlptr->write_stl_file(stlpltfile);
                    Print() << "Writing STL " << specs.stls[stli].name << "\n";
                }
                specs.stls[stli].stlptr->writeVTK(stlpltfile);
                specs.stls[stli].stlptr->printForces();
                specs.stls[stli].stlptr->update_bounding_box();

                
            }

            if(specs.using_softwalls)
            {
                bpc.reassignParticles_softwall(); 
            }

            BL_PROFILE_VAR_STOP(outputs);
            if(specs.calc_mass_flow){
                    Real delta_t = time - t_prev;
                    t_prev = time;
                    bpc.Calculate_Total_Mass_MaterialPoints(mass_flow_next, specs.mass_flow_dir, specs.mass_flow_cutoff);
                    Real flow_out = mass_flow_prev - mass_flow_next;
                    mass_flow_prev = mass_flow_next;
                    PrintToFile(specs.pltprefix + "Mass_Flow.out")<<time<<"\t"<<delta_t<<"\t"<<flow_out<<"\t"<<mass_flow_next<<"\n";
                }
                bpc.Calculate_Total_Speed_MaterialPoints(total_speed);
                total_speed /= bpc.TotalNumberOfParticles();
        }
        steps++;
    }

    bpc.redist_particles(0,specs.using_softwalls);
    bpc.writeParticles(output_it+1+specs.stepoffset, specs.bonded_sphere_particles, specs.pltprefix);
    const std::string& rstfile = specs.pltprefix + amrex::Concatenate("rst", output_it+1+specs.stepoffset, 5);
    bpc.Checkpoint(rstfile, "particles");
    if(specs.using_softwalls)
    {
        bpc.reassignParticles_softwall(); 
    }
    specs.clear();
    

    amrex::Finalize();
}
