#ifndef BDEM_SPECS_H_
#define BDEM_SPECS_H_

#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <AllChemistries.H>

using namespace amrex;

struct realData
{
  enum 
  { // Particle data stored in p.m_data
    radius = 0,
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
    temperature,
    posx_prvs,
    posy_prvs,
    posz_prvs,
    firstspec,
    count
  };
};

struct intData
{
  enum 
  {
    phase = 0,
    near_softwall,
    count
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
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real mu_wall;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real k_t_wall;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real e_t_wall;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real muR_wall;

    //For test case
    // extern AMREX_GPU_DEVICE_MANAGED amrex::Real fy_test;
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

        Vector<Real> autogen_mincoords; 
        Vector<Real> autogen_maxcoords; 
        Vector<Real> autogen_meanvel; 
        Vector<Real> autogen_fluctuation; 
        int autogen_multi_part_per_cell=0;
        Real autogen_radius; 
        Real autogen_dens;
        Real autogen_temp;
        Vector<Real> autogen_species_massfracs;

        bool restartedcase=false;
        bool stl_geom_present=false;
        bool do_heat_transfer=false;
        bool using_softwalls=false;

        Vector<Real> walltemp_polynomial;
        int walltemp_vardir=0;

        int dynamicstl=0;
        int dynstl_transl_dir=0;
        Real dynstl_transl_vel=0.0;
        int dynstl_rot_dir=0;
        Vector<Real> dynstl_center;
        Real dynstl_rot_vel=0.0;

        int particle_sourcing=0;
        Vector<Real> particle_sourcing_mincoords; 
        Vector<Real> particle_sourcing_maxcoords; 
        Vector<Real> particle_sourcing_meanvel; 
        Vector<Real> particle_sourcing_fluctuation; 
        int particle_sourcing_multi_part_per_cell=0;
        Real particle_sourcing_radius; 
        Real particle_sourcing_dens;
        Real particle_sourcing_temp;
        Real particle_sourcing_interval; 
        Vector<Real> particle_sourcing_species_massfracs;
        Real stop_sourcing_time; 

        RealVect gravity;
        int init_particles_using_file=1;
        Real cfl=0.02;

        int mass_flow=0;

        int max_grid_size=8;
        Real write_output_time=10.0;
	Real screen_output_time;
        Real massflow_output_time=1.0;
        int maxsteps=100;
        Real final_time=100.0;
        int num_redist=25;
        int stepoffset=0;
        std::string restartfilename {""};
        std::string stlfilename {""};

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

        // For 2-particle test case
        Real fy_test = 0.0;

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
            dynstl_center.resize(3,0.0);

            ParmParse pp("bdem");
            pp.getarr("prob_lo",plo);
            pp.getarr("prob_hi",phi);
            pp.getarr("ncells", ncells);
            pp.getarr("is_it_periodic",periodic);
            pp.get("max_grid_size", max_grid_size);

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
        
            pp.query("stepoffset",stepoffset);

            pp.get("max_steps", maxsteps);
            pp.get("final_time",final_time);

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
            pp.query("cfl", cfl);
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
                autogen_temp=NTP_TEMP;
                pp.query("temp_autogen",autogen_temp);
                autogen_species_massfracs[0]=1.0;
                pp.queryarr("autogen_species_massfracs",species_massfracs,0,chemptr->nspecies);
                pp.query("autogen_multiple_particles_per_cell",autogen_multi_part_per_cell);
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
                pp.get("temp_particle_sourcing",particle_sourcing_temp);
                pp.get("particle_sourcing_interval",particle_sourcing_interval);
                pp.get("stop_sourcing_time",stop_sourcing_time);
                pp.query("particle_sourcing_multiple_particles_per_cell",particle_sourcing_multi_part_per_cell);
                particle_sourcing_species_massfracs[0]=1.0;
                pp.queryarr("particle_sourcing_species_massfracs",species_massfracs,0,chemptr->nspecies);
            }

            pp.get("normal_spring_constant",DEM::k_n);
            pp.get("normal_restitution_coefficient",DEM::e_n);
            pp.get("friction_coefficient",DEM::mu);
	    pp.get("tangential_spring_constant",DEM::k_t);
            pp.get("tangential_restitution_coefficient",DEM::e_t);
            pp.get("rolling_friction_coefficient",DEM::muR);
            
            DEM::k_n_wall = DEM::k_n;
            DEM::e_n_wall = DEM::e_n;
            DEM::mu_wall  = DEM::mu;
	    DEM::k_t_wall = DEM::k_t;
            DEM::e_t_wall = DEM::e_t;
            DEM::muR_wall  = DEM::muR;
	    
            pp.query("wall_normal_spring_constant",DEM::k_n_wall);
            pp.query("wall_normal_restitution_coefficient",DEM::e_n_wall);
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
            pp.query("stl_file",stlfilename);
            if(!stlfilename.empty())
            {
                stl_geom_present=true;
            }
            pp.query("dynamicstl",dynamicstl);

            if(!stl_geom_present)
            {
                dynamicstl=0;
            }

            if(dynamicstl>0)
            {
                if(dynamicstl==1)
                {
                    pp.query("dynamicstl_transl_dir",dynstl_transl_dir);
                    pp.query("dynamicstl_transl_vel",dynstl_transl_vel);
                }
                else
                    if(dynamicstl==2)
                    {
                        pp.query("dynamicstl_rot_dir",dynstl_rot_dir);
                        pp.queryarr("dynamicstl_rot_center",dynstl_center);
                        pp.query("dynamicstl_rot_vel",dynstl_rot_vel);
                    }
                    else
                    {
                        amrex::Abort("dynamic stl movement function not implemented yet\n");
                    }
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


            //For 2-particle test case only
            pp.query("force_per_dt", fy_test);
            //2-particle test case related up here
            

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
            delete(chemptr);
        }
};
#endif
