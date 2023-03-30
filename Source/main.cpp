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
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::contact_angle = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::gamma    = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::E_bond  = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::G_bond  = zero;
AMREX_GPU_DEVICE_MANAGED amrex::Real DEM::beta_bond  = zero;

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
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

        const int ng_cells = one;
        BDEMParticleContainer bpc(geom, dm, ba, ng_cells,specs.chemptr);
        if(!specs.restartedcase)
        {
            if(specs.init_particles_using_file == 1)
            {
                if(specs.bonded_sphere_particles){
                    bpc.InitBondedParticles("particle_input.dat",specs.do_heat_transfer, specs.cantilever_beam_test);
                } else {
                    bpc.InitParticles("particle_input.dat",specs.do_heat_transfer, specs.glued_sphere_particles, specs.temp_mean, specs.temp_stdev, specs.contact_law);
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
                        specs.autogen_multi_part_per_cell, specs.autogen_max_sphere);
            }
        }
        else
        {
            bpc.Restart(specs.restartfilename, "particles");
        }

        if(specs.modify_species)
        {
            bpc.InitChemSpecies(specs.modified_species_massfracs.data());
        }
        if(specs.stratified_specs)
        {
           bpc.InitChemSpecies(specs.nstrat_domains,specs.strat_mincoords.data(),
                               specs.strat_maxcoords.data(),specs.strat_spec_massfracs.data());
        }

        if(specs.stl_geom_present)
        {
            STLtools::read_stl_file(specs.stlfilename);
        }

        //compute tcoll here over all particles
        DEM::tcoll=bpc.compute_coll_timescale();
        Print()<<"tcoll:"<<DEM::tcoll<<"\n";

        Real dt = specs.cfl*DEM::tcoll;
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
        // FIXME: How do we handle bonded sphere particles that are partially outside boundary?
        if(EBtools::using_levelset_geometry)
        {
            bpc.removeParticlesOutsideBoundary(EBtools::lsphi,
                                               EBtools::ebfactory,EBtools::ls_refinement);
        }
        amrex::Print() << "Num particles after eb removal  " << bpc.TotalNumberOfParticles() << "\n";
        if(specs.stl_geom_present)
        {
            bpc.removeParticlesInsideSTL(specs.outside_point);
        }
        amrex::Print() << "Num particles after stl removal " << bpc.TotalNumberOfParticles() << "\n";

        bpc.set_domain_bcs(specs.bclo,specs.bchi);
        if(specs.visualize_component_spheres && specs.glued_sphere_particles){
            // If using glued sphere model, create new particle container to create particles for each 
            // component sphere (for visualization purposes)
            BDEMParticleContainer bpc_vis(geom, dm, ba, ng_cells,specs.chemptr);
            bpc_vis.createGluedSpheres(bpc);
            bpc_vis.writeParticles(steps+specs.stepoffset, specs.glued_sphere_particles);
        } else {
            bpc.writeParticles(steps+specs.stepoffset, specs.glued_sphere_particles);
        }
        if(specs.dynamicstl!=0)
        {
            std::string stlpltfile = amrex::Concatenate("triplt", steps+specs.stepoffset, 5)+".stl";
            STLtools::write_stl_file(stlpltfile);
        }
        amrex::Print() << "Num particles after init is " << bpc.TotalNumberOfParticles() << "\n";

        // Calculate the moisture content for each particle
        if(specs.liquid_bridging) bpc.computeMoistureContent(specs.moisture_content, specs.contact_angle, specs.liquid_density, specs.FSP);

        while((steps < specs.maxsteps) and (time < specs.final_time))
        {
            time += dt;
            output_time += dt;
            output_timeMass += dt;
            output_timePrint += dt;
            particle_sourcing_time += dt;
        
            Real cb_force = 0.0;
            Real cb_torq = 0.0;
            if(steps>0) specs.init_force = 0.0;
            if(specs.cantilever_beam_test){ 
                cb_force = (time < specs.cb_load_time) ? (time/specs.cb_load_time)*specs.cb_force_max:
                                                         (time < specs.cb_unload_time) ? specs.cb_force_max:0.0;
                cb_torq = (time < specs.cb_load_time) ?  (time/specs.cb_load_time)*specs.cb_torq_max:
                                                         (time < specs.cb_unload_time) ? specs.cb_torq_max:0.0;
            }

            if(specs.particle_sourcing==1 && 
               particle_sourcing_time > specs.particle_sourcing_interval
               && time < specs.stop_sourcing_time)
            {
                bpc.InitParticles (specs.particle_sourcing_mincoords.data(),specs.particle_sourcing_maxcoords.data(), 
                                   specs.particle_sourcing_meanvel.data(),  specs.particle_sourcing_fluctuation.data(), 
                                   specs.particle_sourcing_radius, specs.particle_sourcing_dens,
                                   specs.particle_sourcing_E, specs.particle_sourcing_nu,
                                   specs.particle_sourcing_temp,
                                   specs.particle_sourcing_species_massfracs.data(),
                                   specs.particle_sourcing_multi_part_per_cell, specs.particle_sourcing_max_sphere);

                amrex::Print() << "Num particles before eb removal  " << bpc.TotalNumberOfParticles() << "\n";
                //if(EBtools::using_levelset_geometry and !specs.restartedcase)
                if(EBtools::using_levelset_geometry)
                {
                    bpc.removeParticlesOutsideBoundary(EBtools::lsphi,
                                                       EBtools::ebfactory,EBtools::ls_refinement);
                }
                amrex::Print() << "Num particles after eb removal  " << bpc.TotalNumberOfParticles() << "\n";
                if(specs.stl_geom_present)
                {
                    bpc.removeParticlesInsideSTL(specs.outside_point);
                }
                amrex::Print() << "Num particles after stl removal " << bpc.TotalNumberOfParticles() << "\n";

                bpc.redist_particles(1,specs.using_softwalls);
                amrex::Print()<<"adding particles\n";
                if(specs.using_softwalls)
                {
                    bpc.reassignParticles_softwall(); 
                }
                amrex::Print() << "Num particles after sourcing " << bpc.TotalNumberOfParticles() << "\n";
                particle_sourcing_time=zero;
        
                // FIXME: Just calculating moisture for all particles, but should only be for newly sourced particles
                if(specs.liquid_bridging) bpc.computeMoistureContent(specs.moisture_content, specs.contact_angle, specs.liquid_density, specs.FSP);
            }

            if (steps % specs.num_redist == 0)
            {
                bpc.redist_particles(1,specs.using_softwalls);
                if(specs.using_softwalls)
                {
                    bpc.reassignParticles_softwall(); 
                }
            } 
            else
            {
                bpc.updateNeighbors();
            }
    

            /*if(specs.stl_geom_present)
              {
              bpc.checkParticlesInsideSTL(specs.outside_point);
              }*/

            BL_PROFILE_VAR("FORCE_CALC",forceCalc);
            {
                bpc.computeForces(dt,EBtools::ebfactory,EBtools::lsphi,
                                  specs.do_heat_transfer,specs.walltemp_vardir,
                                  specs.walltemp_polynomial.data(),
                                  EBtools::ls_refinement,specs.stl_geom_present, specs.contact_law, steps,
                                  specs.gravity,
                                  specs.glued_sphere_particles,
                                  specs.bonded_sphere_particles,
                                  specs.liquid_bridging, 
                                  specs.init_force, specs.init_force_dir, specs.init_force_comp,
                                  cb_force, cb_torq, specs.cb_dir);
            }
            BL_PROFILE_VAR_STOP(forceCalc);

            BL_PROFILE_VAR("MOVE_PART",movepart);
            bpc.moveParticles(dt,specs.do_chemistry,specs.minradius_frac,specs.glued_sphere_particles);
            BL_PROFILE_VAR_STOP(movepart);

            if(specs.dynamicstl!=0)
            {
                if(specs.dynamicstl==1)
                {
                    STLtools::move_stl(dt,specs.dynamicstl,specs.dynstl_transl_dir,
                                       specs.dynstl_center.data(),specs.dynstl_transl_vel);
                }
                else
                    if(specs.dynamicstl==2)
                    {
                        STLtools::move_stl(dt,specs.dynamicstl,specs.dynstl_rot_dir,
                                           specs.dynstl_center.data(),specs.dynstl_rot_vel);
                    }
                    else
                    {
                        amrex::Abort("STL motion not implemented\n");
                    }
                if (steps % specs.num_redist == 0)
                {
                    STLtools::boxmap.clear();
                    STLtools::orb_of_triangulation(0,STLtools::num_tri-1,0,STLtools::sorted_indexarray);
                }
            }

            if (output_timePrint > specs.screen_output_time)
            {
                Print()<<"step:"<<steps<<"\t"<<"time:"<<time<<"\n";
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
                if(specs.visualize_component_spheres && specs.glued_sphere_particles){
                    BDEMParticleContainer bpc_vis(geom, dm, ba, ng_cells,specs.chemptr);
                    bpc_vis.createGluedSpheres(bpc);
                    bpc_vis.writeParticles(output_it+specs.stepoffset, specs.glued_sphere_particles);
                } else {
                    bpc.writeParticles(output_it+specs.stepoffset, specs.glued_sphere_particles);
                }
                const std::string& rstfile = amrex::Concatenate("rst", output_it+specs.stepoffset, 5);
                bpc.Checkpoint(rstfile, "particles");
                output_time=zero;
                if(specs.dynamicstl!=0)
                {
                    std::string stlpltfile = amrex::Concatenate("triplt", output_it+specs.stepoffset, 5)+".stl";
                    STLtools::write_stl_file(stlpltfile);
                    STLtools::update_bounding_box();
                }
                if(specs.using_softwalls)
                {
                    bpc.reassignParticles_softwall(); 
                }

                BL_PROFILE_VAR_STOP(outputs);
            }
            steps++;
        }

        bpc.redist_particles(0,specs.using_softwalls);
        if(specs.visualize_component_spheres && specs.glued_sphere_particles){
            BDEMParticleContainer bpc_vis(geom, dm, ba, ng_cells,specs.chemptr);
            bpc_vis.createGluedSpheres(bpc);
            bpc_vis.writeParticles(output_it+1+specs.stepoffset, specs.glued_sphere_particles);
        } else {
            bpc.writeParticles(output_it+1+specs.stepoffset, specs.glued_sphere_particles);
        }
        const std::string& rstfile = amrex::Concatenate("rst", output_it+1+specs.stepoffset, 5);
        bpc.Checkpoint(rstfile, "particles");
        if(specs.using_softwalls)
        {
            bpc.reassignParticles_softwall(); 
        }
        specs.clear();
    }

    amrex::Finalize();
}
