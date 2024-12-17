#ifndef BDEM_SPECS_H_
#define BDEM_SPECS_H_

#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <AllChemistries.H>
#include <STLtools.H>

using namespace amrex;

struct realData
{
  enum 
  { // Particle data stored in p.m_data
    radius,
    radinit,
    xvel,
    yvel,
    zvel,
    fx,
    fy,
    fz,
    xangvel,
    yangvel,
    zangvel,
    taux,
    tauy,
    tauz,
    Iinv,
    volume,
    mass,
    density,
    E,
    nu,
    temperature,
    posx_prvs,
    posy_prvs,
    posz_prvs,
    euler_angle_x,
    euler_angle_y,
    euler_angle_z,
    liquid_volume,   
    total_bridge_volume,
    first_bond_v, // We leave room for 9 real components (ft_x, ft_y, ft_z, taun_x, taun_y, taun_z, taut_x, taut_y, taut_z) * MAXBONDS
    theta_x = first_bond_v + MAXBONDS*9,
    firstspec,
    count
  };
};

struct intData
{
  enum 
  {
    phase = 0,
    near_softwall = 1,
    type_id = 2,
    unique_id = 3,
    first_bridge=4,
    first_bond = first_bridge + MAXBRIDGES,
    count = first_bond + MAXBONDS
  };
};


namespace DEM
{
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real k_n;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real e_n;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real mu;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real k_t;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real e_t;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real muR;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real tcoll;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real k_n_wall;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real e_n_wall;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real E_wall;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real nu_wall;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real mu_wall;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real k_t_wall;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real e_t_wall;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real muR_wall;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real global_damping;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real force_damping;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real angv_damping;

    extern AMREX_GPU_DEVICE_MANAGED amrex::Real mu_liq;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real hminf;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real CED;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real contact_angle;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real gamma;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real wall_gamma;

    extern AMREX_GPU_DEVICE_MANAGED amrex::Real E_bond;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real G_bond;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real beta_bond;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real nu_bond;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real bond_radius_factor;

    extern AMREX_GPU_DEVICE_MANAGED amrex::Real k_c;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real sigma_max;

    extern AMREX_GPU_DEVICE_MANAGED amrex::Real eps_g;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho_g;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real mu_g;

    extern AMREX_GPU_DEVICE_MANAGED int particles_in_parcel;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real CG_ratio;
}

class DEMspecs
{
    public:
        Vector<Real> plo;
        Vector<Real> phi;
        Vector<Real> outside_point;
        Vector<int> ncells;
        Vector<int> periodic;

        Vector<int> bclo;
        Vector<int> bchi;
        int contact_law = 0;
        int ng_cells = one;

        Vector<Real> autogen_mincoords; 
        Vector<Real> autogen_maxcoords; 
        Vector<Real> autogen_meanvel; 
        Vector<Real> autogen_fluctuation; 
        int autogen_multi_part_per_cell=0;
        int autogen_layer_particles=0;
        Real autogen_radius; 
        Real autogen_dens;
        Real autogen_E = 100.0;
        Real autogen_nu = 0.3;
        Real autogen_temp;
        Real autogen_min_radius = -1.0;
        Real autogen_max_radius = -1.0;
        Vector<Real> autogen_species_massfracs;
        int autogen_bp_type = 0;

        bool restartedcase=false;
        bool do_heat_transfer=false;
        bool using_softwalls=false;
        bool remove_eb_overlaps=false;
        bool reassign_particle_properties=false;

        Vector<Real> walltemp_polynomial;
        int walltemp_vardir=0;

        // STL stuff
        


        int particle_sourcing=0;
        Vector<Real> particle_sourcing_mincoords; 
        Vector<Real> particle_sourcing_maxcoords; 
        Vector<Real> particle_sourcing_meanvel; 
        Vector<Real> particle_sourcing_fluctuation; 
        int particle_sourcing_multi_part_per_cell=0;
        int particle_sourcing_layer_particles=0;
        Real particle_sourcing_radius; 
        Real particle_sourcing_dens;
        Real particle_sourcing_E = 100.0;
        Real particle_sourcing_nu = 0.3;
        Real particle_sourcing_temp;
        Real particle_sourcing_interval=1.0; 
        Vector<Real> particle_sourcing_species_massfracs;
        Real stop_sourcing_time; 
        Real particle_sourcing_min_radius = -1.0;
        Real particle_sourcing_max_radius = -1.0;
        int particle_sourcing_bp_type = 0;

        RealVect gravity;
        int init_particles_using_file=1;
        Real cfl=0.02;
        Real constant_dt = 0.0;
        int verlet_scheme = 0;
        Real avg_speed_stop = -1;

        int mass_flow=0;

        int max_grid_size=8;
        Real write_output_time=10.0;
	      Real screen_output_time;
        Real massflow_output_time=1.0;
        int maxsteps=100;
        Real final_time=100.0;
        int num_redist=25;
        int num_upd_neighbor=1;
        int stepoffset=0;
        std::string restartfilename {""};
	std::string pltprefix {""};

        mechanism *chemptr;
        std::string chemistry;
        Vector<Real> species_massfracs;
        int do_chemistry=0;
        Real minradius_frac=0.1;

        int modify_species=0;
        Vector<Real> modified_species_massfracs;

        int stratified_specs=0;
        int nstrat_domains;
        Vector<Real> strat_mincoords; 
        Vector<Real> strat_maxcoords;
        Vector<Real> strat_spec_massfracs; 

        // input options for bonded sphere particles
        int bonded_sphere_particles = 0;
        Real bond_damping_factor = 0.0;
        Vector<int> particle_type_list;  // Generate particles by randomly sampling from user-specified list
        int use_type_dist = 0;
        Vector<Real> particle_dist_list;  // List of percentages for each particle type

        // Options for modeling liquid bridging and viscous liquid forces
        int liquid_bridging = 0;
        Real moisture_content = 0.0;
        Real moisture_content_stdev = 0.0;
        Real FSP = 23.0;     // Fiber saturation point of pine chips by default
        Real contact_angle = 0.0;
        Real surface_tension = 0.0;
        Real liquid_density = 1000.0;   // Assume water density by default
        int recalculate_MC = 0;

        Real init_force = 0.0;
        int init_force_dir = 0;
        int init_force_comp = 0;

        // Cantilever beam test case variables
        int cantilever_beam_test = 0;
        Real cb_force_max = 0.0;
        Real cb_torq_max = 0.0;
        int cb_dir = 1;
        Real cb_load_time = 1.0;
        Real cb_unload_time = 1e20;
        bool cb_intervals = false;
        int cb_interval_num = 0;
        Real cb_hold_time = 0.0;

        // Options for initializing particle temperatures according to normal dist.
        Real temp_mean = -1.0;
        Real temp_stdev = -1.0;

        // Options for clipping particles
        int clip_particles = 0;
        int clip_particle_dir;
        Real clip_particle_val;

        int calc_mass_flow = 0;
        int mass_flow_dir = 1;
        Real mass_flow_cutoff = 0.0;

        // Generic cohesion model variables
        int particle_cohesion = 0;
        Real k_c = 0.2;     // Ratio of normal to cohesive forces
        Real sigma_max = 10.0;    // Maximum cohesive stress

        // Options for creating particle color stripes
        int coloring_striation = 0;
        Real coloring_striation_len = 1.0;
        int coloring_striation_dir = 1;

        // One-way drag modeling
        int drag_model = 0;

        // List of STL files and their settings
        Vector< stl_specs > stls;
        bool stl_geom_present=false;

        void read_dem_specs()
        {
            plo.resize(AMREX_SPACEDIM);
            phi.resize(AMREX_SPACEDIM);
            ncells.resize(AMREX_SPACEDIM);
            periodic.resize(AMREX_SPACEDIM);
            outside_point.resize(AMREX_SPACEDIM);
            
            autogen_mincoords.resize(AMREX_SPACEDIM);
            autogen_maxcoords.resize(AMREX_SPACEDIM);
            autogen_meanvel.resize(AMREX_SPACEDIM);
            autogen_fluctuation.resize(AMREX_SPACEDIM);
            
            particle_sourcing_mincoords.resize(AMREX_SPACEDIM);
            particle_sourcing_maxcoords.resize(AMREX_SPACEDIM);
            particle_sourcing_meanvel.resize(AMREX_SPACEDIM);
            particle_sourcing_fluctuation.resize(AMREX_SPACEDIM);
           
            bclo.resize(AMREX_SPACEDIM);
            bchi.resize(AMREX_SPACEDIM);
            species_massfracs.resize(MAXSPECIES);
            modified_species_massfracs.resize(MAXSPECIES);
            particle_sourcing_species_massfracs.resize(MAXSPECIES,0.0);
            autogen_species_massfracs.resize(MAXSPECIES,0.0);

            walltemp_polynomial.resize(3,0.0);


            ParmParse pp("bdem");
            pp.getarr("prob_lo",plo);
            pp.getarr("prob_hi",phi);
            pp.getarr("ncells", ncells);
            pp.getarr("is_it_periodic",periodic);
            pp.get("max_grid_size", max_grid_size);
            pp.query("ng_cells", ng_cells);

            //outside_point for stl stuff
            Real fac=3.0;
            outside_point[0]=plo[0]+fac*(phi[0]-plo[0]);
            outside_point[1]=plo[1]+fac*(phi[1]-plo[1]);
            outside_point[2]=plo[2]+fac*(phi[2]-plo[2]);
            pp.queryarr("outside_point",outside_point);
            
         
            //by default it is outflow 
            bclo[0]=0;
            bclo[1]=0;
            bclo[2]=0;

            bchi[0]=0;
            bchi[1]=0;
            bchi[2]=0;
            
            pp.queryarr("bc_lower",bclo);
            pp.queryarr("bc_upper",bchi);

            if(bclo[0]==SOFTWALL_BC || bclo[1]==SOFTWALL_BC || bclo[2]==SOFTWALL_BC ||
               bchi[0]==SOFTWALL_BC || bchi[1]==SOFTWALL_BC || bchi[2]==SOFTWALL_BC)
            {
                using_softwalls=true;
            }

            pp.query("remove_eb_overlaps", remove_eb_overlaps);

            // Contact law
            pp.query("contact_law", contact_law);
        
            pp.query("stepoffset",stepoffset);

            pp.get("max_steps", maxsteps);
            pp.get("final_time",final_time);
            pp.query("avg_speed_stop", avg_speed_stop);

	    pp.query("particle_mass_history",mass_flow);
	    if(mass_flow)
	    {
		massflow_output_time=write_output_time/10.0;
		pp.query("write_output_time_massflow", massflow_output_time);
	    }

            write_output_time=final_time/10.0;
            pp.query("write_output_time", write_output_time);

	    screen_output_time=write_output_time/10.0;
	    pp.query("screen_output_time", screen_output_time);

            pp.query("num_redist", num_redist);
            pp.query("num_upd_neighbor", num_upd_neighbor);
            pp.query("cfl", cfl);
            pp.query("constant_dt", constant_dt);
            pp.query("verlet_scheme", verlet_scheme);
            pp.query("init_particles_using_file",init_particles_using_file);
            
            chemistry="nullchem";
            pp.query("chemistry",chemistry);
            if(chemistry == "nullchem")
            {
                nullchem* nullchemptr;
                nullchemptr=new nullchem;
                chemptr=nullchemptr;
            }
            else if(chemistry == "plasticschem")
            {
                plasticschem* plastchemptr;
                plastchemptr=new plasticschem; 
                chemptr=plastchemptr;
            }
            else if(chemistry == "plasticsmixchem")
            {
                plasticsmixchem* plastmixchemptr;
                plastmixchemptr=new plasticsmixchem; 
                chemptr=plastmixchemptr;
            }
            else
            {
                amrex::Print()<<"Chemistry specified:"<<chemistry<<"\n";
                amrex::Abort("Chemistry not implemented yet\n");
            }
            chemptr->define(); 
            species_massfracs[0]=1.0;
            pp.queryarr("species_massfracs",species_massfracs,0,chemptr->nspecies);
            pp.query("do_chemistry",do_chemistry);
            pp.query("minradius_frac",minradius_frac);
        
            pp.query("modify_species",modify_species);
            modified_species_massfracs[0]=1.0;
            if(modify_species)
            {
                pp.queryarr("modified_species_massfracs"
                            ,modified_species_massfracs,0,chemptr->nspecies);
            }
	    
            if(init_particles_using_file==0)
            {
                pp.getarr("mincoords_autogen",autogen_mincoords);
                pp.getarr("maxcoords_autogen",autogen_maxcoords);
                pp.getarr("meanvel_autogen",autogen_meanvel);
                pp.getarr("fluctuation_autogen",autogen_fluctuation);
                pp.get("radius_autogen",autogen_radius);
                pp.get("dens_autogen",autogen_dens);
                if(contact_law != 0){
                    pp.get("E_autogen",autogen_E);
                    pp.get("nu_autogen",autogen_nu);
                }
                autogen_temp=NTP_TEMP;
                pp.query("temp_autogen",autogen_temp);
                autogen_species_massfracs[0]=1.0;
                pp.queryarr("autogen_species_massfracs",species_massfracs,0,chemptr->nspecies);
                pp.query("autogen_multiple_particles_per_cell",autogen_multi_part_per_cell);
                pp.query("autogen_layer_particles",autogen_layer_particles);
                pp.query("min_radius_autogen",autogen_min_radius);
                pp.query("max_radius_autogen",autogen_max_radius);
            }

            pp.query("particle_sourcing",particle_sourcing);
            if(particle_sourcing==1)
            {
                pp.getarr("mincoords_particle_sourcing",particle_sourcing_mincoords);
                pp.getarr("maxcoords_particle_sourcing",particle_sourcing_maxcoords);
                pp.getarr("meanvel_particle_sourcing",particle_sourcing_meanvel);
                pp.getarr("fluctuation_particle_sourcing",particle_sourcing_fluctuation);
                pp.get("radius_particle_sourcing",particle_sourcing_radius);
                pp.get("dens_particle_sourcing",particle_sourcing_dens);
                if(contact_law != 0){
                    pp.get("E_particle_sourcing",particle_sourcing_E);
                    pp.get("nu_particle_sourcing",particle_sourcing_nu);
                }
                particle_sourcing_temp=NTP_TEMP;
                pp.query("temp_particle_sourcing",particle_sourcing_temp);
                pp.get("particle_sourcing_interval",particle_sourcing_interval);
                pp.get("stop_sourcing_time",stop_sourcing_time);
                pp.query("particle_sourcing_multiple_particles_per_cell",particle_sourcing_multi_part_per_cell);
                pp.query("particle_sourcing_layer_particles",particle_sourcing_layer_particles);
                particle_sourcing_species_massfracs[0]=1.0;
                pp.queryarr("particle_sourcing_species_massfracs",species_massfracs,0,chemptr->nspecies);
                pp.query("min_radius_particle_sourcing",particle_sourcing_min_radius);
                pp.query("max_radius_particle_sourcing",particle_sourcing_max_radius);
            }

            pp.get("normal_spring_constant",DEM::k_n);
            pp.get("normal_restitution_coefficient",DEM::e_n);
            pp.get("friction_coefficient",DEM::mu);
	        pp.get("tangential_spring_constant",DEM::k_t);
            pp.get("tangential_restitution_coefficient",DEM::e_t);
            pp.get("rolling_friction_coefficient",DEM::muR);
            
            DEM::k_n_wall = DEM::k_n;
            DEM::e_n_wall = DEM::e_n;
            DEM::E_wall  = 100.0;  // Initial unrealistic value
            DEM::nu_wall  = 0.3;  // Default general nu_wall
            DEM::mu_wall  = DEM::mu;
	        DEM::k_t_wall = DEM::k_t;
            DEM::e_t_wall = DEM::e_t;
            DEM::muR_wall  = DEM::muR;
	    
            pp.query("wall_normal_spring_constant",DEM::k_n_wall);
            pp.query("wall_normal_restitution_coefficient",DEM::e_n_wall);
            pp.query("elastic_modulus_wall",DEM::E_wall);
            pp.query("poissons_ratio_wall",DEM::nu_wall);
            pp.query("wall_friction_coefficient",DEM::mu_wall);
	        pp.query("wall_tangential_spring_constant",DEM::k_t_wall);
            pp.query("wall_tangential_restitution_coefficient",DEM::e_t_wall);
            pp.query("wall_rolling_friction_coefficient",DEM::muR_wall);

            Array<Real,AMREX_SPACEDIM> gravity_input{0.0, 0.0, 0.0};
            pp.get("gravity", gravity_input);
            for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
            {
                    gravity[dir] = gravity_input[dir];
            }
            
            pp.query("restart",restartfilename);
            if(!restartfilename.empty())
            {
                restartedcase=true;
            }

            /* STL related stuff -------------------------------*/
            int number_of_stl_files(0);
            pp.query("number_of_stl_files", number_of_stl_files);

            if ( number_of_stl_files > 0 )
            {
                stl_geom_present = true;
            }

            stls.resize(number_of_stl_files);

            //- Read stl file names
            Vector<std::string> stlnames(number_of_stl_files);
            pp.getarr("stl_names",stlnames);

            //- Now read stl settings independently
            for (int stli = 0; stli < number_of_stl_files; stli++)
            {
                std::string stlname(stlnames[stli]);

                stls[stli].name = stlname;

                stls[stli].stlptr = new STLtools();
                stls[stli].stlptr->name = stls[stli].name;
                pp.query( ( stlname + "_stl_timeoffset" ).c_str(), stls[stli].stl_timeoffset );

                if (stls[stli].stl_timeoffset > 0)
                {
                    stls[stli].stlptr->read_stl_file(stls[stli].name + std::to_string(stls[stli].stl_timeoffset) + ".stl");
                }
                else
                {
                    stls[stli].stlptr->read_stl_file(stls[stli].name  + ".stl");
                }

                stls[stli].stlptr->buildGridData();
                pp.query( ( stlname + "_dynamicstl" ).c_str(), stls[stli].dynamicstl );

                amrex::Print()  << "Found STL named " << stls[stli].name << "\n";

                if ( stls[stli].dynamicstl == 1 )
                {
                    stls[stli].dynstl_transl_dir.resize(3,0.);
                    pp.getarr( ( stlname + "_dynamicstl_transl_dir" ).c_str(), stls[stli].dynstl_transl_dir );
                    pp.get( ( stlname + "_dynamicstl_transl_vel" ).c_str(), stls[stli].dynstl_transl_vel );

                    amrex::Print()  << "STL has dynamic translation with velocity " 
                                    << stls[stli].dynstl_transl_vel 
                                    << ", and direction (" 
                                    << stls[stli].dynstl_transl_dir[0] << " "
                                    << stls[stli].dynstl_transl_dir[1] << " "
                                    << stls[stli].dynstl_transl_dir[2] << ") "
                                    << "\n";
                }                
                else if ( stls[stli].dynamicstl == 2 )
                {
                    stls[stli].dynstl_center.resize(3,0.);
                    stls[stli].dynstl_rot_dir.resize(3,0.);
                    bool use_geo_center;
                    pp.query( ( stlname + "_dynamicstl_geo_rot_center" ).c_str(), use_geo_center );
                    pp.get( ( stlname + "_dynamicstl_rot_vel" ).c_str(), stls[stli].dynstl_rot_vel );
                    pp.getarr ( (stlname + "_dynamicstl_rot_dir" ).c_str(), stls[stli].dynstl_rot_dir );

                    // Normalize the rotation direction
                    amrex::Real mod = normVec(stls[stli].dynstl_rot_dir.data());
                    for (int dir = 0; dir < 3; dir++)
                    {
                        stls[stli].dynstl_rot_dir.data()[dir] /= mod;
                    }
                
                    // Check if the geometric center of the STL should be used as center of rotation
                    if( use_geo_center )
                    {
                        stls[stli].stlptr->getSTLGeoCenter(stls[stli].dynstl_center.data());
                    }
                    else
                    {
                        pp.getarr( ( stlname + "_dynamicstl_rot_center" ).c_str(), stls[stli].dynstl_center );
                    }

                    amrex::Print()  << "STL has dynamic rotation with center ( " 
                                    << stls[stli].dynstl_center[0] << " "
                                    << stls[stli].dynstl_center[1] << " "
                                    << stls[stli].dynstl_center[2] << " )"
                                    << ", velocity " << stls[stli].dynstl_rot_vel 
                                    << ", and direction (" 
                                    << stls[stli].dynstl_rot_dir[0] << " "
                                    << stls[stli].dynstl_rot_dir[1] << " "
                                    << stls[stli].dynstl_rot_dir[2] << " )"<< "\n";
                }
                else if ( stls[stli].dynamicstl == 3 )
                {
                    stls[stli].stl_vib_dir.resize(3,0.);
                    pp.get( ( stlname + "_stl_vib_amp" ).c_str(), stls[stli].stl_vib_amp );
                    pp.get( ( stlname + "_stl_vib_freq" ).c_str(), stls[stli].stl_vib_freq );
                    pp.getarr( ( stlname + "_stl_vib_dir" ).c_str(), stls[stli].stl_vib_dir );
                }

            }
            

            // Reassign particle properties for a restarted case using autogen particle parameters
            pp.query("reassign_particle_properties", reassign_particle_properties);
            pp.query("coloring_striation", coloring_striation);
            if(coloring_striation){
                pp.get("coloring_striation_len", coloring_striation_len);
                pp.get("coloring_striation_dir", coloring_striation_dir);
            }

            pp.query("do_heat_transfer",do_heat_transfer);

            if(do_heat_transfer)
            {
                walltemp_vardir=0;
                pp.query("walltemp_vardir",walltemp_vardir);

                //quadratic polynomial: c0+c1*x+c2*x^2
                walltemp_polynomial[0]=NTP_TEMP;
                walltemp_polynomial[1]=0.0;
                walltemp_polynomial[2]=0.0;
                pp.queryarr("walltemp_polynomial",walltemp_polynomial);
            }

            pp.query("stratified_specs",stratified_specs);
            if(stratified_specs)
            {
                pp.get("nstrat_domains",nstrat_domains);
                strat_mincoords.resize(AMREX_SPACEDIM*nstrat_domains);
                strat_maxcoords.resize(AMREX_SPACEDIM*nstrat_domains);
                strat_spec_massfracs.resize(MAXSPECIES*nstrat_domains);

                for(int i=0;i<nstrat_domains;i++)
                {
                    std::string mincoordstr=amrex::Concatenate("strat_mincoords_",i,2); 
                    std::string maxcoordstr=amrex::Concatenate("strat_maxcoords_",i,2); 
                    std::string specstr=amrex::Concatenate("strat_spec_",i,2);

                    Vector<Real> tmpvec1,tmpvec2,tmpvec3;
                    tmpvec1.resize(AMREX_SPACEDIM);
                    tmpvec2.resize(AMREX_SPACEDIM);
                    tmpvec3.resize(chemptr->nspecies);

                    pp.getarr(mincoordstr.c_str(),tmpvec1,0,AMREX_SPACEDIM);
                    pp.getarr(maxcoordstr.c_str(),tmpvec2,0,AMREX_SPACEDIM);
                    pp.getarr(specstr.c_str(),tmpvec3,0,chemptr->nspecies);
                    for(int d=0;d<AMREX_SPACEDIM;d++)
                    {
                        strat_mincoords[i*AMREX_SPACEDIM+d]=tmpvec1[d];
                        strat_maxcoords[i*AMREX_SPACEDIM+d]=tmpvec2[d];
                    }
                    for(int sp=0;sp<chemptr->nspecies;sp++)
                    {
                        strat_spec_massfracs[i*chemptr->nspecies+sp]=tmpvec3[sp];
                    }
                }
            }

            pp.query("liquid_bridging", liquid_bridging);
            if(liquid_bridging){
                // Read in liquid bridging variables
                pp.query("moisture_content", moisture_content); // Moisture content using Green (or wet) basis
                pp.query("moisture_content_stdev", moisture_content_stdev); // Set to zero if since MC value desired
                pp.query("contact_angle", DEM::contact_angle);
                pp.query("FSP", FSP);
                pp.query("liquid_density", liquid_density);
                pp.query("liquid_viscosity", DEM::mu_liq);
                pp.query("surface_tension", DEM::gamma);
                pp.query("wall_surface_tension", DEM::wall_gamma);
                pp.query("recalculate_MC", recalculate_MC);
                pp.query("hminf", DEM::hminf);
                pp.query("cohesion_energy_density", DEM::CED);

                // Convert contact angle to rad
                DEM::contact_angle *= PI/180.0;

                // Convert MC and FSP from percentages
                FSP /= 100.0; 
                moisture_content /= 100.0;
                moisture_content_stdev /= 100.0;

                // TODO: add more input checks
                if(moisture_content < 0.0 || moisture_content > 100.0) Abort("\nMoisture content should be a percentage between 0 and 100.\n");
            }

            // Spin test variables
            pp.query("init_force", init_force);
            pp.query("init_force_dir", init_force_dir);
            pp.query("init_force_comp", init_force_comp);

            // Normal dist. temperature variables
            pp.query("temp_mean", temp_mean);
            pp.query("temp_stdev", temp_stdev);

            // Bonded sphere model variables
            pp.query("bonded_sphere_particles", bonded_sphere_particles);
            if(bonded_sphere_particles){
                pp.get("young_mod_bond", DEM::E_bond);
                pp.get("poisson_rat_bond", DEM::nu_bond);
                pp.query("bond_damping_factor", DEM::beta_bond);
                pp.query("bond_radius_factor", DEM::bond_radius_factor);
                DEM::G_bond = DEM::E_bond / (2.0 * (1.0 + DEM::nu_bond));

                pp.query("cantilever_beam_test", cantilever_beam_test);
                if(cantilever_beam_test){
                    if(restartedcase || !init_particles_using_file || particle_sourcing) Abort("\nCantilever beam test case not compatible with restart, autogen, or sourcing!\n"); 
                    pp.query("cb_force_max", cb_force_max);
                    pp.query("cb_torq_max", cb_torq_max);
                    pp.get("cb_dir", cb_dir);
                    pp.get("cb_load_time", cb_load_time);
                    pp.query("cb_unload_time", cb_unload_time);
                    pp.query("cb_intervals", cb_intervals);
                    if(cb_intervals){
                      pp.get("cb_interval_num", cb_interval_num);
                      pp.get("cb_hold_time", cb_hold_time);
                    }
                }
            }

            if(bonded_sphere_particles){
                pp.query("autogen_bp_type", autogen_bp_type);
                pp.query("particle_sourcing_bp_type", particle_sourcing_bp_type);
                if (pp.contains("particle_type_list")){
                    int num_types = pp.countval("particle_type_list");
                    particle_type_list.resize(num_types);
                    pp.getarr("particle_type_list", particle_type_list);
                    pp.query("use_type_dist", use_type_dist);
                    if(use_type_dist){
                        int num_types_dist = pp.countval("particle_dist_list");
                        if(num_types != num_types_dist) amrex::Abort("Particle type distribution list length much match particle type list length!\n");
                        particle_dist_list.resize(num_types_dist);
                        pp.getarr("particle_dist_list", particle_dist_list);
                        
                        // Normalize particle distribution list
                        Real dsum = 0.0;
                        for(int i=0; i<num_types; i++) dsum += particle_dist_list[i];
                        if(dsum <= 0.0) amrex::Abort("Invalid distribution values supplied!\n");
                        for(int i=0; i<num_types; i++) particle_dist_list[i] = particle_dist_list[i] / dsum;
                    }
                }
            }

            pp.query("global_damping", DEM::global_damping);
            pp.query("force_damping", DEM::force_damping);
            pp.query("angv_damping", DEM::angv_damping);

            pp.query("clip_particles", clip_particles);
            if(clip_particles){
              pp.get("clip_particle_dir", clip_particle_dir);
              pp.get("clip_particle_val", clip_particle_val);
            }

            pp.query("calc_mass_flow", calc_mass_flow);
            if(calc_mass_flow){
                pp.get("mass_flow_dir", mass_flow_dir);
                pp.get("mass_flow_cutoff", mass_flow_cutoff);
            }

            pp.query("particle_cohesion", particle_cohesion);
            if(particle_cohesion){
                pp.get("cohesion_coefficient", DEM::k_c);
                pp.get("cohesion_max_stress", DEM::sigma_max);
            }

            pp.query("drag_model", drag_model);
            if(drag_model == 1){
                pp.get("gas_vfrac", DEM::eps_g);
                pp.get("gas_visc", DEM::mu_g);
                pp.get("gas_density", DEM::rho_g);
            }

            //- CG options
            // Ratio between the radius of the Coarse Grained particle and
            // the radius of the original particle
            //! Notice that the current implementation follow the work of
            //! Nakamura et al., Coarse-Grained discrete element method
            //! for shear granular flows, Chem. Eng. J. Adv. (2020)  

            pp.query("CG_ratio",DEM::CG_ratio);      
            
            if ( DEM::CG_ratio > 1.0 )
            {
                DEM::particles_in_parcel = std::pow( DEM::CG_ratio, 3.0 );

                amrex::Print()  << "Using Coarse Grained particles with one parcel representing "
                                << DEM::particles_in_parcel << " particles.\n"
                                << "The corresponding radius ratio is "
                                << "( r_CG / r_p ) = " << DEM::CG_ratio << "\n";
                amrex::Print()  << "WARNING: the CG model only works when the stiffness k is specified in the input! \n";

                //! See table 1 of Nakamura et al.

                //- Scale the stiffness
                DEM::k_n        *= ( DEM::CG_ratio * DEM::CG_ratio ) / (  ( 2. * DEM::CG_ratio ) - 1. );
                DEM::k_n_wall    = ( DEM::CG_ratio * DEM::CG_ratio ) / ( ( 1. / DEM::k_n_wall )  + ( ( DEM::CG_ratio - 1. ) / DEM::k_n ) );
                DEM::k_t        *= ( DEM::CG_ratio * DEM::CG_ratio ) / (  ( 2. * DEM::CG_ratio ) - 1. );
                DEM::k_t_wall    = ( DEM::CG_ratio * DEM::CG_ratio ) / ( ( 1. / DEM::k_t_wall )  + ( ( DEM::CG_ratio - 1. ) / DEM::k_t ) );
                
                //- Scale the restitution coefficient
                DEM::e_n        =   std::sqrt( 1. - (1./DEM::CG_ratio) * ( 1. - (DEM::e_n*DEM::e_n) ) );
                DEM::e_n_wall   =   std::sqrt( 1. - (1./DEM::CG_ratio) * ( 1. - (DEM::e_n_wall*DEM::e_n_wall) ) );
                DEM::e_t        =   std::sqrt( 1. - (1./DEM::CG_ratio) * ( 1. - (DEM::e_t*DEM::e_t) ) );
                DEM::e_t_wall   =   std::sqrt( 1. - (1./DEM::CG_ratio) * ( 1. - (DEM::e_t_wall*DEM::e_t_wall) ) );

                //- The sliding friction coefficient cannot be pre-scaled and it is computed during collisions
            }

        }
        void clear()
        {
            plo.clear();
            phi.clear();
            ncells.clear();
            periodic.clear();
            autogen_mincoords.clear(); 
            autogen_maxcoords.clear(); 
            autogen_meanvel.clear(); 
            autogen_fluctuation.clear();
            for (int stli = 0; stli < stls.size(); stli++)
            {
                delete stls[stli].stlptr;
            }
            stls.clear();
             
            delete(chemptr);
        }
};
#endif
