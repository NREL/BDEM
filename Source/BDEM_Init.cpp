#include <BDEM_ParticleContainer.H>
#include <BDEM_Collide_Utils.H>
#include <stl_tools/STLtools.H>
#include <AMReX_Random.H>
#include <BDEM_BondedParticles.H>

void BDEMParticleContainer::InitParticles (const std::string& filename,
                                           bool &do_heat_transfer,
                                           int glued_sphere_particles, 
                                           int glued_sphere_types,
                                           Real temp_mean, Real temp_stdev,
                                           int contact_law, int liquid_bridging,
                                           Real liquid_density, Real MC_avg, 
                                           Real MC_stdev, Real FSP)
{

    // only read the file on the IO proc
    if (ParallelDescriptor::IOProcessor())  
    {
        std::ifstream ifs;
        ifs.open(filename.c_str(), std::ios::in);

        if (!ifs.good())
        {
            amrex::FileOpenFailed(filename);
        }

        int np = -1;
        ifs >> np >> std::ws;

        // Issue an error if nparticles = 0 is specified
        if ( np == -1 )
        {
            Abort("\nCannot read number of particles from particle_input.dat: file is corrupt.\
                    \nPerhaps you forgot to specify the number of particles on the first line??? ");
        }

        // we add all the particles to grid 0 and tile 0 and let
        // Redistribute() put them in the right places.
        const int lev  = 0;
        const int grid = 0;
        const int tile = 0;

        Gpu::HostVector<ParticleType> host_particles;
        
        auto& particle_tile = DefineAndReturnParticleTile(lev,grid,tile);
        int p_types[BP_TYPES] = {BP_NP0, BP_NP1, BP_NP2, BP_NP3, BP_NP4,
                                 BP_NP5, BP_NP6, BP_NP7, BP_NP8, BP_NP9,
                                 BP_NP10, BP_NP11, BP_NP12, BP_NP13, BP_NP14,
                                 BP_NP15, BP_NP16, BP_NP17, BP_NP18, BP_NP19, BP_NP20};

        for (int i = 0; i < np; i++) 
        {
            ParticleType p;
            // Set id and cpu for this particle
            p.id()  = ParticleType::NextID();
            p.cpu() = ParallelDescriptor::MyProc();

            // Read from input file
            // NOTE: Assumed that all glued sphere particle radii are equal for time being
            ifs >> p.idata(intData::phase);
            ifs >> p.pos(0);
            ifs >> p.pos(1);
            ifs >> p.pos(2);
            ifs >> p.rdata(realData::radius);
            ifs >> p.rdata(realData::density);
            if(contact_law == 1){
                ifs >> p.rdata(realData::E);
                ifs >> p.rdata(realData::nu);
            } else {
                p.rdata(realData::E) = 100.0;
                p.rdata(realData::nu) = 0.3;
            }
            ifs >> p.rdata(realData::xvel);
            ifs >> p.rdata(realData::yvel);
            ifs >> p.rdata(realData::zvel);

            if(do_heat_transfer)
            {
                ifs >> p.rdata(realData::temperature);

                // Overwrite temperature in file if a valid mean temperature and temp stdev are specified
                if(temp_mean > 0.0 && temp_stdev > 0.0){
                    p.rdata(realData::temperature) = max(amrex::RandomNormal(temp_mean, temp_stdev),1.0);
                }
            }
            else
            {
               p.rdata(realData::temperature)=NTP_TEMP;
            }
            
           

            if(glued_sphere_particles)
            {
                // Extracting number of component spheres and Euler angles for glued-sphere cylinders
                // TODO: How do we make sure cylinders don't intersect upon initialization?
                if(!glued_sphere_types) ifs >> p.idata(intData::num_comp_sphere);
                ifs >> p.rdata(realData::euler_angle_x);
                ifs >> p.rdata(realData::euler_angle_y);
                ifs >> p.rdata(realData::euler_angle_z);

                // Use Euler angles to calculate initial quaternion components
                Real eax = p.rdata(realData::euler_angle_x);
                Real eay = p.rdata(realData::euler_angle_y);
                Real eaz = p.rdata(realData::euler_angle_z);
                p.rdata(realData::q0) =  cos(eaz/2.0)*cos(eay/2.0)*cos(eax/2.0) + sin(eaz/2.0)*sin(eay/2.0)*sin(eax/2.0);
                p.rdata(realData::q1) = -sin(eaz/2.0)*sin(eay/2.0)*cos(eax/2.0) + cos(eaz/2.0)*cos(eay/2.0)*sin(eax/2.0);
                p.rdata(realData::q2) =  cos(eaz/2.0)*sin(eay/2.0)*cos(eax/2.0) + sin(eaz/2.0)*cos(eay/2.0)*sin(eax/2.0);
                p.rdata(realData::q3) =  sin(eaz/2.0)*cos(eay/2.0)*cos(eax/2.0) - cos(eaz/2.0)*sin(eay/2.0)*sin(eax/2.0);

                // Calculate principal axis components in inertial frame (for visualization)
                Real pa_body[THREEDIM] = {1.0, 0.0, 0.0};
                Real pa_inert[THREEDIM];
                rotate_vector_to_inertial(p, pa_body, pa_inert);
                p.rdata(realData::pax) = pa_inert[XDIR];
                p.rdata(realData::pay) = pa_inert[YDIR];
                p.rdata(realData::paz) = pa_inert[ZDIR];
            } else{
                p.idata(intData::num_comp_sphere) = 1;
                p.rdata(realData::euler_angle_x) = zero;
                p.rdata(realData::euler_angle_y) = zero;
                p.rdata(realData::euler_angle_z) = zero;
                p.rdata(realData::q0) = 1.0;
                p.rdata(realData::q1) = zero;
                p.rdata(realData::q2) = zero;
                p.rdata(realData::q3) = zero;
                p.rdata(realData::pax) = zero;
                p.rdata(realData::pay) = zero;
                p.rdata(realData::paz) = zero;
                p.idata(intData::type_id) = 0;
            }

            //set initial radius
            p.rdata(realData::radinit)=p.rdata(realData::radius);

            // Set other particle properties
            // NOTE: Calculation of particle mass for glued sphere particles assumes no component sphere overlap
            p.rdata(realData::volume)      = fourbythree*PI*pow(p.rdata(realData::radius),three)*p.idata(intData::num_comp_sphere);
            p.rdata(realData::mass)        = p.rdata(realData::density)*p.rdata(realData::volume);
            p.rdata(realData::Iinv)        = 2.5/(p.rdata(realData::mass)*pow(p.rdata(realData::radius),two));
            
            // User input for angular velocity when using glued sphere particles (at least for testing purposes)
            if(glued_sphere_particles){
                ifs >> p.rdata(realData::xangvel);
                ifs >> p.rdata(realData::yangvel);
                ifs >> p.rdata(realData::zangvel);
                if(glued_sphere_types) {
                  ifs >> p.idata(intData::type_id);
                  p.idata(intData::num_comp_sphere) = p_types[p.idata(intData::type_id)];
                }
            } else {
                p.rdata(realData::xangvel)     = zero;
                p.rdata(realData::yangvel)     = zero;
                p.rdata(realData::zangvel)     = zero;
            }

            if(glued_sphere_particles){
                // NOTE: Approximating moments of inertia for glued sphere particle using cylinder formula
                p.rdata(realData::Ixinv) = 2.0/(p.rdata(realData::mass)*(pow(p.rdata(realData::radius),two)) );
                p.rdata(realData::Iyinv) = 12.0/(p.rdata(realData::mass)*(3.0*pow(p.rdata(realData::radius),two)+pow(2.0*p.idata(intData::num_comp_sphere)*p.rdata(realData::radius),two)));
                p.rdata(realData::Izinv) = 12.0/(p.rdata(realData::mass)*(3.0*pow(p.rdata(realData::radius),two)+pow(2.0*p.idata(intData::num_comp_sphere)*p.rdata(realData::radius),two)));

            } else {
                p.rdata(realData::Ixinv) = 2.5/(p.rdata(realData::mass)*pow(p.rdata(realData::radius),two));
                p.rdata(realData::Iyinv) = 2.5/(p.rdata(realData::mass)*pow(p.rdata(realData::radius),two));
                p.rdata(realData::Izinv) = 2.5/(p.rdata(realData::mass)*pow(p.rdata(realData::radius),two));
            }

            p.rdata(realData::fx) = zero;
            p.rdata(realData::fy) = zero;
            p.rdata(realData::fz) = zero;
            p.rdata(realData::taux) = zero;
            p.rdata(realData::tauy) = zero;
            p.rdata(realData::tauz) = zero;
            p.rdata(realData::theta_x) = zero;

            // Set bond components to zero
            for(int b=0; b<MAXBONDS*9; b++){
                p.rdata(realData::first_bond_v+b) = zero;
            }

            // Set bridge indices to -1 to indicate no existing bridges
            for(int br=0; br<MAXBRIDGES; br++){
                p.idata(intData::first_bridge+3*br) = -1;
                p.idata(intData::first_bridge+3*br+1) = -1;
                p.idata(intData::first_bridge+3*br+2) = -1;
            }

            // If using liquid bridging, calculate particle liquid and recalculate particle mass and density
            if(liquid_bridging){
                Real MC = (MC_stdev > 0.0) ? min(max(amrex::RandomNormal(MC_avg, MC_stdev),0.0),0.9):MC_avg;
                p.rdata(realData::liquid_volume) = (MC > FSP) ? (p.rdata(realData::density)*p.rdata(realData::volume)/liquid_density)*(MC - FSP)/(1 - MC):0;
                p.rdata(realData::mass) = p.rdata(realData::density)*p.rdata(realData::volume) * (1.0 + MC/(1.0 - MC));
                p.rdata(realData::density) = p.rdata(realData::mass) / p.rdata(realData::volume);
            } else {
                p.rdata(realData::liquid_volume) = zero;
            }

            // Keep track of how much particle liquid volume is already used to form bridges
            p.rdata(realData::total_bridge_volume) = zero;

            for(int b=0; b<MAXBONDS; b++) p.idata(intData::first_bond + b) = -1;
            
            //FIXME: get chemistry data from inputs file
            for(int sp=0;sp<MAXSPECIES;sp++)
            {
                p.rdata(realData::firstspec+sp)=0.0;
            }

            // Add everything to the data structure
            host_particles.push_back(p);

            if (!ifs.good())
            {
                amrex::Abort("Error initializing particles from Ascii file. \n");
            }
        }
        
        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + host_particles.size();
        particle_tile.resize(new_size);

        Gpu::copy(Gpu::hostToDevice,
                  host_particles.begin(),
                  host_particles.end(),
                  particle_tile.GetArrayOfStructs().begin() + old_size);
    }
    
    Redistribute();
}

void BDEMParticleContainer::InitBondedParticles (const std::string& filename,
                                                 bool &do_heat_transfer,
                                                 int contact_law,
                                                 int cantilever_beam_test,
                                                 int liquid_bridging, Real liquid_density,
                                                 Real MC_avg, Real MC_stdev, Real FSP)
{

    // only read the file on the IO proc
    if (ParallelDescriptor::IOProcessor())  
    {
        std::ifstream ifs;
        ifs.open(filename.c_str(), std::ios::in);

        if (!ifs.good())
        {
            amrex::FileOpenFailed(filename);
        }

        int np = -1;
        ifs >> np >> std::ws;

        // Issue an error if nparticles = 0 is specified
        if ( np == -1 )
        {
            Abort("\nCannot read number of particles from particle_input.dat: file is corrupt.\
                    \nPerhaps you forgot to specify the number of particles on the first line??? ");
        }
        if(np != 1 && cantilever_beam_test)
        {
            Abort("\nOne particle must be specified for the cantilever beam test case!\n");
        } 

        // we add all the particles to grid 0 and tile 0 and let
        // Redistribute() put them in the right places.
        const int lev  = 0;
        const int grid = 0;
        const int tile = 0;

        Gpu::HostVector<ParticleType> host_particles;
        
        auto& particle_tile = DefineAndReturnParticleTile(lev,grid,tile);

        const ParticleBondData bp_data = ParticleBondData();
        int bp_types[BP_TYPES] = {BP_NP0, BP_NP1, BP_NP2, BP_NP3, BP_NP4, 
                                  BP_NP5, BP_NP6, BP_NP7, BP_NP8, BP_NP9, 
                                  BP_NP10, BP_NP11, BP_NP12, BP_NP13, BP_NP14, 
                                  BP_NP15, BP_NP16, BP_NP17, BP_NP18, BP_NP19, BP_NP20};

        Real pc_pos[THREEDIM];                    
        int bp_ids[200];

        for (int i = 0; i < np; i++) 
        {
            // Read in standard particle properties
            int bp_phase;
            Real bp_pos[THREEDIM];
            Real bp_radius, bp_density, bp_temperature, bp_E, bp_nu;
            Real bp_vel[THREEDIM];

            ifs >> bp_phase;
            ifs >> bp_pos[XDIR];
            ifs >> bp_pos[YDIR];
            ifs >> bp_pos[ZDIR];
            ifs >> bp_radius;
            ifs >> bp_density;
            if(contact_law == 1){
                ifs >> bp_E;
                ifs >> bp_nu;
            } else {
                bp_E = 100.0;
                bp_nu = 0.3;
            }
            ifs >> bp_vel[XDIR];
            ifs >> bp_vel[YDIR];
            ifs >> bp_vel[ZDIR];
            if(do_heat_transfer){
                ifs >> bp_temperature;
            } else {
                bp_temperature = NTP_TEMP;
            }
            amrex::Real eax, eay, eaz;
            ifs >> eax;
            ifs >> eay;
            ifs >> eaz;

            // Last entry should be the bonded particle type
            int bp_type;
            ifs >> bp_type;

            amrex::Real bp_q[4] ={cos(eaz/2.0)*cos(eay/2.0)*cos(eax/2.0) + sin(eaz/2.0)*sin(eay/2.0)*sin(eax/2.0),
                                  -sin(eaz/2.0)*sin(eay/2.0)*cos(eax/2.0) + cos(eaz/2.0)*cos(eay/2.0)*sin(eax/2.0),
                                   cos(eaz/2.0)*sin(eay/2.0)*cos(eax/2.0) + sin(eaz/2.0)*cos(eay/2.0)*sin(eax/2.0),
                                   sin(eaz/2.0)*cos(eay/2.0)*cos(eax/2.0) - cos(eaz/2.0)*sin(eay/2.0)*sin(eax/2.0)};

            // Use the same moisture content percentage for each component particle
            Real MC = 0.0;
            if(liquid_bridging) MC = (MC_stdev > 0.0) ? min(max(amrex::RandomNormal(MC_avg, MC_stdev),0.0),0.9):MC_avg;

            // TODO: checks: appropriate particle for cantilever, valid ID selected, max bonds is sufficient
            for(int j = 0; j<bp_types[bp_type]; j++) bp_ids[j] = ParticleType::NextID();
            for(int j = 0; j<bp_types[bp_type]; j++){
                ParticleType p;
                p.id() = bp_ids[j];
                if(cantilever_beam_test && j == bp_types[bp_type]-1) bp_phase = -1;    // Left-most particle is held inert
                get_bonded_particle_pos(bp_type, j, bp_radius, bp_pos, bp_q, pc_pos);
                bp_init(p, bp_data, bp_phase, pc_pos, bp_radius, bp_density, bp_vel, bp_temperature, j, bp_type, bp_ids, liquid_density, MC, FSP, bp_E, bp_nu);
                host_particles.push_back(p);
            } 
            if (!ifs.good())
            {
                amrex::Abort("Error initializing particles from Ascii file. \n");
            }
        }
        
        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + host_particles.size();
        particle_tile.resize(new_size);

        Gpu::copy(Gpu::hostToDevice,
                  host_particles.begin(),
                  host_particles.end(),
                  particle_tile.GetArrayOfStructs().begin() + old_size);
    }
    
    Redistribute();
}

void BDEMParticleContainer::InitChemSpecies(amrex::Real Yis[MAXSPECIES])
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const auto plo = Geom(lev).ProbLoArray();
    const auto phi = Geom(lev).ProbHiArray();
    const auto dx = Geom(lev).CellSizeArray();
    auto& plev  = GetParticles(lev);

    amrex::Real Yi_captured[MAXSPECIES];

    for(int i=0;i<MAXSPECIES;i++)
    {
        Yi_captured[i]=Yis[i];
    }

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numParticles();
        ParticleType* pstruct = aos().dataPtr();

        // now we move the particles
        amrex::ParallelFor(np,[=]
                AMREX_GPU_DEVICE (int i) noexcept
        {
                ParticleType& p = pstruct[i];

                for(int sp=0;sp<MAXSPECIES;sp++)
                {
                    //p.rdata(realData::firstspec+i)=Yis[i];
                    p.rdata(realData::firstspec+sp)=Yi_captured[sp];
                }
        });
    }
}

void BDEMParticleContainer::InitChemSpecies(int ndomains, Real *mincoords,
        Real *maxcoords,Real *spec_massfracs)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const auto plo = Geom(lev).ProbLoArray();
    const auto phi = Geom(lev).ProbHiArray();
    const auto dx = Geom(lev).CellSizeArray();
    auto& plev  = GetParticles(lev);

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numParticles();
        ParticleType* pstruct = aos().dataPtr();

        // now we move the particles
        amrex::ParallelFor(np,[=]
                AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];

            Real x=p.pos(0);
            Real y=p.pos(1);
            Real z=p.pos(2);

            for(int d=0;d<ndomains;d++)
            {
                if(x>mincoords[d*AMREX_SPACEDIM+0] && x<=maxcoords[d*AMREX_SPACEDIM+0] &&
                        y>mincoords[d*AMREX_SPACEDIM+1] && y<=maxcoords[d*AMREX_SPACEDIM+1] &&
                        z>mincoords[d*AMREX_SPACEDIM+2] && z<=maxcoords[d*AMREX_SPACEDIM+2])
                {
                    for(int sp=0;sp<MAXSPECIES;sp++)
                    {
                        p.rdata(realData::firstspec+sp)=spec_massfracs[d*MAXSPECIES+sp];
                    }

                }

            }

        });
    }
}

void BDEMParticleContainer::removeParticlesOutsideBoundary(const MultiFab *lsmfab,
        const EBFArrayBoxFactory *ebfactory,
        const int ls_refinement,
        const int glued_sphere_particles)
{
    const int lev = 0;
    auto& plev  = GetParticles(lev);
    const FabArray<EBCellFlagFab>* flags = &(ebfactory->getMultiEBCellFlagFab());
    const Geometry& geom = Geom(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numParticles();

        ParticleType* pstruct = aos().dataPtr();

        const Box & bx = mfi.tilebox();

        if ((*flags)[mfi].getType(bx) == FabType::covered)
        {
            amrex::ParallelFor(np,[=]
            AMREX_GPU_DEVICE (int i) noexcept
            {
                ParticleType& p = pstruct[i];
                p.id()=-1;   
            });
        }
        else
        {
            const auto phiarr = lsmfab->array(mfi);
            amrex::ParallelFor(np,[=]
                  AMREX_GPU_DEVICE (int i) noexcept
            {
                ParticleType& p = pstruct[i];
                Real rp = p.rdata(realData::radius);
                for(int pc=0; pc<p.idata(intData::num_comp_sphere); pc++){
                    Real ppos_inert[THREEDIM];
                    if(glued_sphere_particles){
                        get_inertial_pos(p, pc, ppos_inert); 
                    } else {
                        ppos_inert[XDIR] = p.pos(0);
                        ppos_inert[YDIR] = p.pos(1);
                        ppos_inert[ZDIR] = p.pos(2);
                    }
                    Real ls_value = get_levelset_value(ppos_inert, ls_refinement, phiarr, plo, dx);
                    if(ls_value < 0.0)
                    {
                        p.id()=-1;   
                    }
                }
            });
        }
    }
    Redistribute();
}

void BDEMParticleContainer::removeParticlesInsideSTL(Vector<Real> outside_point, const int glued_sphere_particles)
{
    const int lev = 0;
    auto& plev  = GetParticles(lev);
    const Geometry& geom = Geom(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    Array<Real,AMREX_SPACEDIM> po_arr{outside_point[0],outside_point[1],outside_point[2]};

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numParticles();

        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];
            for(int pc = 0; pc<p.idata(intData::num_comp_sphere); pc++){
  
                Real ploc[THREEDIM];
                if(glued_sphere_particles){
                    get_inertial_pos(p, pc, ploc); 
                } else {
                    ploc[XDIR] = p.pos(0);
                    ploc[XDIR] = p.pos(1);
                    ploc[XDIR] = p.pos(2);
                }
                Real ploc_t[3]={0.0};
                Real t1[3],t2[3],t3[3];
                Real outp[]={po_arr[0],po_arr[1],po_arr[2]};
                int num_intersects=0;

                for(int dim=0;dim<3;dim++)
                {
                    for(int j=0;j<3;j++)
                    {
                        ploc_t[dim] += STLtools::eigdirs[3*dim+j]*ploc[j];
                    }
                }

                if( (ploc_t[0]>STLtools::bbox_lo[0]) && 
                   (ploc_t[0]<STLtools::bbox_hi[0]) &&
                   (ploc_t[1]>STLtools::bbox_lo[1]) &&
                   (ploc_t[1]<STLtools::bbox_hi[1]) &&
                   (ploc_t[2]>STLtools::bbox_lo[2]) &&
                   (ploc_t[2]<STLtools::bbox_hi[2]) )
                {
                    for(int tr=0;tr<STLtools::num_tri;tr++)
                    {
                        t1[0]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+0];
                        t1[1]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+1];
                        t1[2]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+2];

                        t2[0]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+3];
                        t2[1]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+4];
                        t2[2]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+5];

                        t3[0]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+6];
                        t3[1]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+7];
                        t3[2]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+8];

                        num_intersects += (1-STLtools::lineseg_tri_intersect(outp,ploc,t1,t2,t3));
                    }
                    if(num_intersects%2 == 1)
                    {
                        p.id()=-1;   
                    }
                }
            }
        });
    }

    Redistribute();
}

void BDEMParticleContainer::checkParticlesInsideSTL(Vector<Real> outside_point, const int glued_sphere_particles)
{
    const int lev = 0;
    auto& plev  = GetParticles(lev);
    const Geometry& geom = Geom(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    Array<Real,AMREX_SPACEDIM> po_arr{outside_point[0],outside_point[1],outside_point[2]};

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numParticles();

        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];
            for(int pc = 0; pc<p.idata(intData::num_comp_sphere); pc++){
                Real ploc[THREEDIM];
                if(glued_sphere_particles){
                    get_inertial_pos(p, pc, ploc); 
                } else {
                    ploc[XDIR] = p.pos(0);
                    ploc[XDIR] = p.pos(1);
                    ploc[XDIR] = p.pos(2);
                }
                Real ploc_t[3]={0.0};
                Real t1[3],t2[3],t3[3];
                Real outp[]={po_arr[0],po_arr[1],po_arr[2]};
                int num_intersects=0;

                for(int dim=0;dim<3;dim++)
                {
                    for(int j=0;j<3;j++)
                    {
                        ploc_t[dim] += STLtools::eigdirs[3*dim+j]*ploc[j];
                    }
                }

                if( (ploc_t[0]>STLtools::bbox_lo[0]) && 
                   (ploc_t[0]<STLtools::bbox_hi[0]) &&
                   (ploc_t[1]>STLtools::bbox_lo[1]) &&
                   (ploc_t[1]<STLtools::bbox_hi[1]) &&
                   (ploc_t[2]>STLtools::bbox_lo[2]) &&
                   (ploc_t[2]<STLtools::bbox_hi[2]) )
                {
                    for(int tr=0;tr<STLtools::num_tri;tr++)
                    {
                        t1[0]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+0];
                        t1[1]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+1];
                        t1[2]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+2];

                        t2[0]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+3];
                        t2[1]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+4];
                        t2[2]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+5];

                        t3[0]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+6];
                        t3[1]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+7];
                        t3[2]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+8];

                        num_intersects += (1-STLtools::lineseg_tri_intersect(outp,ploc,t1,t2,t3));
                    }
                    if(num_intersects%2 == 1)
                    {
#ifndef AMREX_USE_GPU
                        amrex::Print()<<"particle inside stl:"<<ploc[0]<<"\t"<<ploc[1]<<"\t"<<ploc[2]<<"\n";
#else
                        amrex::Abort("particle inside stl");
#endif

                    }
                }
            }
        });
    }

}


void BDEMParticleContainer::InitParticles (Real mincoords[THREEDIM],Real maxcoords[THREEDIM], 
                                           Real meanvel[THREEDIM], Real fluctuation[THREEDIM], Real rad, Real dens,
                                           Real E, Real nu, Real temp,
                                           Real spec[MAXSPECIES],
                                           int do_multi_part_per_cell, 
                                           int layer_particles,
                                           int glued_sphere_particles,
                                           int glued_sphere_types,
                                           int bonded_sphere_particles,
                                           int min_sphere, int max_sphere,
                                           Real min_rad, Real max_rad, 
                                           int p_type, Vector<int> type_list,
                                           int liquid_bridging,
                                           Real liquid_density, Real MC_avg, 
                                           Real MC_stdev, Real FSP)
{
    int lev = 0;
    Real x,y,z,x0,y0,z0;

    Real dx = Geom(lev).CellSize(0);
    Real dy = Geom(lev).CellSize(1);
    Real dz = Geom(lev).CellSize(2);
    const Real* plo = Geom(lev).ProbLo();

    std::mt19937 mt(0451);
    std::uniform_real_distribution<double> dist(0.4, 0.6);
    std::uniform_real_distribution<double> dist1(0.15, 0.35);
    std::uniform_real_distribution<double> dist2(0.65, 0.85);

    const ParticleBondData p_data = ParticleBondData();
    int p_types[BP_TYPES] = {BP_NP0, BP_NP1, BP_NP2, BP_NP3, BP_NP4, 
                             BP_NP5, BP_NP6, BP_NP7, BP_NP8, BP_NP9, 
                             BP_NP10, BP_NP11, BP_NP12, BP_NP13, BP_NP14, 
                             BP_NP15, BP_NP16, BP_NP17, BP_NP18, BP_NP19, BP_NP20};

    Real pc_pos[THREEDIM];                    
    Real bp_pos[THREEDIM];
    int bp_ids[200];
    int bp_phase = 0;

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) 
    {

        const Box& tile_box = mfi.tilebox();
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id, tile_id)];

        Gpu::HostVector<ParticleType> host_particles;

        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) 
        {
            if(do_multi_part_per_cell == 0 || layer_particles > 0)
            {
                int layers = (layer_particles > 1) ? layer_particles:1;

                for(int pl=0; pl<layers; pl++){

                    // Particles layered in the y-direction
                    x = plo[XDIR] + (iv[XDIR] + dist(mt))*dx;
                    y = (layers > 1) ? plo[YDIR] + (iv[YDIR] + ((1.0+pl)/(1.0+layers)))*dy:plo[YDIR] + (iv[YDIR] + dist(mt))*dy;
                    z = plo[ZDIR] + (iv[ZDIR] + dist(mt))*dz;

                    if(x>=mincoords[XDIR] && x<=maxcoords[XDIR] &&
                       y>=mincoords[YDIR] && y<=maxcoords[YDIR] &&
                       z>=mincoords[ZDIR] && z<=maxcoords[ZDIR])
                    {
                        Real MC = 0.0;
                        if(liquid_bridging) MC = (MC_stdev > 0.0) ? min(max(amrex::RandomNormal(MC_avg, MC_stdev),0.0),0.9):MC_avg;
                        if(bonded_sphere_particles){
                            int type = (p_type == -1) ? ceil(amrex::Random()*BP_TYPES) -1:p_type;
                            if(type_list.size() > 0){
                                type = type_list[ceil(amrex::Random()*type_list.size()) - 1];
                                if ( type < 0 || type > BP_TYPES-1) Abort("\nInvalid particle type specified .\n");
                            }
                            bp_pos[XDIR] = x; bp_pos[YDIR] = y; bp_pos[ZDIR] = z;
                            // Real eax = amrex::Random()*PI/2.0;
                            // Real eay = amrex::Random()*PI/2.0;
                            // Real eaz = amrex::Random()*PI/2.0;
                            Real eax = PI/2.0;
                            Real eay = amrex::Random()*PI/20.0;
                            Real eaz = amrex::Random()*PI/20.0;
                            amrex::Real quats[4] ={cos(eaz/2.0)*cos(eay/2.0)*cos(eax/2.0) + sin(eaz/2.0)*sin(eay/2.0)*sin(eax/2.0),
                                                  -sin(eaz/2.0)*sin(eay/2.0)*cos(eax/2.0) + cos(eaz/2.0)*cos(eay/2.0)*sin(eax/2.0),
                                                   cos(eaz/2.0)*sin(eay/2.0)*cos(eax/2.0) + sin(eaz/2.0)*cos(eay/2.0)*sin(eax/2.0),
                                                   sin(eaz/2.0)*cos(eay/2.0)*cos(eax/2.0) - cos(eaz/2.0)*sin(eay/2.0)*sin(eax/2.0)};
                            Real bp_vel[THREEDIM] = {meanvel[XDIR] + fluctuation[XDIR]*(amrex::Random()-half),
                                                     meanvel[YDIR] + fluctuation[XDIR]*(amrex::Random()-half),
                                                     meanvel[ZDIR] + fluctuation[XDIR]*(amrex::Random()-half)};
                            for(int j = 0; j<p_types[type]; j++) bp_ids[j] = ParticleType::NextID();
                            for(int j = 0; j<p_types[type]; j++){
                                ParticleType p;
                                p.id() = bp_ids[j];
                                get_bonded_particle_pos(type, j, rad, bp_pos, quats, pc_pos);
                                bp_init(p, p_data, bp_phase, pc_pos, rad, dens, bp_vel, temp, j, type, bp_ids, liquid_density, MC, FSP, E, nu);
                                host_particles.push_back(p);
                            } 
                        } else if(glued_sphere_particles){
                            int type = (!glued_sphere_types) ? -2:(p_type == -1) ? ceil(amrex::Random()*BP_TYPES) -1:p_type;
                            if(type_list.size() > 0){
                                type = type_list[ceil(amrex::Random()*type_list.size()) - 1];
                                if ( type < 0 || type > BP_TYPES-1) Abort("\nInvalid particle type specified .\n");
                            }
                            int ncs = (glued_sphere_types) ? p_types[type]:min_sphere + floor(amrex::Random()*(max_sphere+1 - min_sphere));
                            ParticleType p = generate_particle(x,y,z,
                                                               meanvel[XDIR] + fluctuation[XDIR]*(amrex::Random()-half),
                                                               meanvel[YDIR] + fluctuation[YDIR]*(amrex::Random()-half),
                                                               meanvel[ZDIR] + fluctuation[ZDIR]*(amrex::Random()-half),
                                                               dens, rad, E, nu, 
                                                               temp, spec, liquid_density, MC, FSP, ncs, type);
                            host_particles.push_back(p);
                        } else {
                            // Calculate radius using random distribution if min and max values are specified
                            Real p_rad = (min_rad > 0.0 && max_rad > 0.0 && max_rad > min_rad) ? min_rad + (max_rad - min_rad)*amrex::Random():rad;

                            ParticleType p = generate_particle(x,y,z,
                                                               meanvel[XDIR] + fluctuation[XDIR]*(amrex::Random()-half),
                                                               meanvel[YDIR] + fluctuation[YDIR]*(amrex::Random()-half),
                                                               meanvel[ZDIR] + fluctuation[ZDIR]*(amrex::Random()-half),
                                                               dens, p_rad, E, nu, 
                                                               temp, spec, liquid_density, MC, FSP);
                            host_particles.push_back(p);
                        }
                    }
                }
            }
            else
            {
                x0 = plo[XDIR]+iv[XDIR]*dx;
                y0 = plo[YDIR]+iv[YDIR]*dy;
                z0 = plo[ZDIR]+iv[ZDIR]*dz;

                for(int k=0;k<2;k++)
                {
                    for(int j=0;j<2;j++)
                    {
                        for(int i=0;i<2;i++)
                        {
                            //x = x0 + (i+dist(mt))*half*dx;
                            //y = y0 + (j+dist(mt))*half*dy;
                            //z = z0 + (k+dist(mt))*half*dz;
                            // x = x0 + (i+half)*half*dx;
                            // y = y0 + (j+half)*half*dy;
                            // z = z0 + (k+half)*half*dz;
                            x = (i == 0) ? x0 + dist1(mt)*dx : x0 + dist2(mt)*dx;
                            y = (j == 0) ? y0 + dist1(mt)*dy : y0 + dist2(mt)*dy;
                            z = (k == 0) ? z0 + dist1(mt)*dz : z0 + dist2(mt)*dz;

                            if(x>=mincoords[XDIR] and x<=maxcoords[XDIR] and 
                               y>=mincoords[YDIR] and y<=maxcoords[YDIR] and
                               z>=mincoords[ZDIR] and z<=maxcoords[ZDIR])
                            {
                                Real MC = 0.0;
                                if(liquid_bridging) MC = (MC_stdev > 0.0) ? min(max(amrex::RandomNormal(MC_avg, MC_stdev),0.0),0.9):MC_avg;
                                if(bonded_sphere_particles){
                                    int type = (p_type == -1) ? ceil(amrex::Random()*BP_TYPES) -1:p_type;
                                    if(type_list.size() > 0){
                                        type = type_list[ceil(amrex::Random()*type_list.size()) - 1];
                                        if ( type < 0 || type > BP_TYPES-1) Abort("\nInvalid particle type specified .\n");
                                    }
                                    bp_pos[XDIR] = x; bp_pos[YDIR] = y; bp_pos[ZDIR] = z;
                                    Real eax = amrex::Random()*PI/2.0;
                                    Real eay = amrex::Random()*PI/2.0;
                                    Real eaz = amrex::Random()*PI/2.0;
                                    amrex::Real quats[4] ={cos(eaz/2.0)*cos(eay/2.0)*cos(eax/2.0) + sin(eaz/2.0)*sin(eay/2.0)*sin(eax/2.0),
                                                          -sin(eaz/2.0)*sin(eay/2.0)*cos(eax/2.0) + cos(eaz/2.0)*cos(eay/2.0)*sin(eax/2.0),
                                                           cos(eaz/2.0)*sin(eay/2.0)*cos(eax/2.0) + sin(eaz/2.0)*cos(eay/2.0)*sin(eax/2.0),
                                                           sin(eaz/2.0)*cos(eay/2.0)*cos(eax/2.0) - cos(eaz/2.0)*sin(eay/2.0)*sin(eax/2.0)};
                                    Real bp_vel[THREEDIM] = {meanvel[XDIR] + fluctuation[XDIR]*(amrex::Random()-half),
                                                             meanvel[YDIR] + fluctuation[XDIR]*(amrex::Random()-half),
                                                             meanvel[ZDIR] + fluctuation[XDIR]*(amrex::Random()-half)};
                                    for(int pi = 0; pi<p_types[type]; pi++) bp_ids[pi] = ParticleType::NextID();
                                    for(int pi = 0; pi<p_types[type]; pi++){
                                        ParticleType p;
                                        p.id() = bp_ids[pi];
                                        get_bonded_particle_pos(type, pi, rad, bp_pos, quats, pc_pos);
                                        bp_init(p, p_data, bp_phase, pc_pos, rad, dens, bp_vel, temp, pi, type, bp_ids, liquid_density, MC, FSP, E, nu);
                                        host_particles.push_back(p);
                                    } 
                                } else if(glued_sphere_particles){
                                    int type = (!glued_sphere_types) ? -2:(p_type == -1) ? ceil(amrex::Random()*BP_TYPES) -1:p_type;
                                    if(type_list.size() > 0){
                                        type = type_list[ceil(amrex::Random()*type_list.size()) - 1];
                                        if ( type < 0 || type > BP_TYPES-1) Abort("\nInvalid particle type specified .\n");
                                    }
                                    int ncs = (glued_sphere_types) ? p_types[type]:min_sphere + floor(amrex::Random()*(max_sphere+1 - min_sphere));
                                    ParticleType p = generate_particle(x,y,z,
                                                                       meanvel[XDIR] + fluctuation[XDIR]*(amrex::Random()-half),
                                                                       meanvel[YDIR] + fluctuation[YDIR]*(amrex::Random()-half),
                                                                       meanvel[ZDIR] + fluctuation[ZDIR]*(amrex::Random()-half),
                                                                       dens, rad, E, nu, 
                                                                       temp, spec, liquid_density, MC, FSP, ncs, type);
                                    host_particles.push_back(p);
                                } else {
                                    // Calculate radius using random distribution if min and max values are specified
                                    Real p_rad = (min_rad > 0.0 && max_rad > 0.0 && max_rad > min_rad) ? min_rad + (max_rad - min_rad)*amrex::Random():rad;

                                    ParticleType p = generate_particle(x,y,z,
                                                                       meanvel[XDIR] + fluctuation[XDIR]*(amrex::Random()-half),
                                                                       meanvel[YDIR] + fluctuation[YDIR]*(amrex::Random()-half),
                                                                       meanvel[ZDIR] + fluctuation[ZDIR]*(amrex::Random()-half),
                                                                       dens, p_rad, E, nu, 
                                                                       temp, spec, liquid_density, MC, FSP);
                                    host_particles.push_back(p);
                                }
                            }
                        }
                    } 
                }
            }
        }

        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + host_particles.size();
        particle_tile.resize(new_size);

        Gpu::copy(Gpu::hostToDevice,
                  host_particles.begin(),
                  host_particles.end(),
                  particle_tile.GetArrayOfStructs().begin() + old_size);

    }

    // We shouldn't need this if the particles are tiled with one tile per grid, but otherwise
    // we do need this to move particles from tile 0 to the correct tile.
    Redistribute();
}

BDEMParticleContainer::ParticleType BDEMParticleContainer::generate_particle(Real x,Real y,Real z,
                                                                             Real velx, Real vely, Real velz,
                                                                             Real dens, Real rad, Real E, Real nu,
                                                                             Real temp, Real spec[MAXSPECIES], 
                                                                             Real liquid_density, Real MC, Real FSP,
                                                                             int num_sphere, int p_type)
{
    ParticleType p;
    p.id()  = ParticleType::NextID();
    p.cpu() = ParallelDescriptor::MyProc();                

    p.pos(XDIR) = x;
    p.pos(YDIR) = y;
    p.pos(ZDIR) = z;

    p.idata(intData::phase) = 0;
    p.rdata(realData::radius) = rad;
    p.rdata(realData::radinit) = rad;
    p.rdata(realData::E) = E;
    p.rdata(realData::nu) = nu;

    p.rdata(realData::density) = dens;
    p.rdata(realData::xvel) = velx;
    p.rdata(realData::yvel) = vely;
    p.rdata(realData::zvel) = velz;
    p.rdata(realData::temperature)=temp;

    p.rdata(realData::fx) = zero;
    p.rdata(realData::fy) = zero;
    p.rdata(realData::fz) = zero;
    p.rdata(realData::taux) = zero;
    p.rdata(realData::tauy) = zero;
    p.rdata(realData::tauz) = zero;
    p.rdata(realData::theta_x) = zero;

    p.idata(intData::num_comp_sphere) = num_sphere;
    p.rdata(realData::euler_angle_x) = PI/2.0;                 
    p.rdata(realData::euler_angle_y) = amrex::Random()*PI/20.0;
    p.rdata(realData::euler_angle_z) = amrex::Random()*PI/20.0;

    Real eax = p.rdata(realData::euler_angle_x);
    Real eay = p.rdata(realData::euler_angle_y);
    Real eaz = p.rdata(realData::euler_angle_z);
    p.rdata(realData::q0) =  cos(eaz/2.0)*cos(eay/2.0)*cos(eax/2.0) + sin(eaz/2.0)*sin(eay/2.0)*sin(eax/2.0);
    p.rdata(realData::q1) = -sin(eaz/2.0)*sin(eay/2.0)*cos(eax/2.0) + cos(eaz/2.0)*cos(eay/2.0)*sin(eax/2.0);
    p.rdata(realData::q2) =  cos(eaz/2.0)*sin(eay/2.0)*cos(eax/2.0) + sin(eaz/2.0)*cos(eay/2.0)*sin(eax/2.0);
    p.rdata(realData::q3) =  sin(eaz/2.0)*cos(eay/2.0)*cos(eax/2.0) - cos(eaz/2.0)*sin(eay/2.0)*sin(eax/2.0);

    Real pa_body[THREEDIM] = {1.0, 0.0, 0.0};
    Real pa_inert[THREEDIM];
    rotate_vector_to_inertial(p, pa_body, pa_inert);
    p.rdata(realData::pax) = pa_inert[XDIR];
    p.rdata(realData::pay) = pa_inert[YDIR];
    p.rdata(realData::paz) = pa_inert[ZDIR];
    p.idata(intData::type_id) = p_type;

    // Set other particle properties
    p.rdata(realData::volume)      = fourbythree*PI*pow(p.rdata(realData::radius),three)*p.idata(intData::num_comp_sphere);
    p.rdata(realData::mass)        = p.rdata(realData::density)*p.rdata(realData::volume);
    p.rdata(realData::Iinv)        = 2.5/(p.rdata(realData::mass)*pow(p.rdata(realData::radius),two));
    p.rdata(realData::xangvel)     = zero;
    p.rdata(realData::yangvel)     = zero;
    p.rdata(realData::zangvel)     = zero;
    p.rdata(realData::Ixinv) = 2.0/(p.rdata(realData::mass)*(pow(p.rdata(realData::radius),two)) );
    p.rdata(realData::Iyinv) = 12.0/(p.rdata(realData::mass)*(3.0*pow(p.rdata(realData::radius),two)+pow(2.0*p.idata(intData::num_comp_sphere)*p.rdata(realData::radius),two)));
    p.rdata(realData::Izinv) = 12.0/(p.rdata(realData::mass)*(3.0*pow(p.rdata(realData::radius),two)+pow(2.0*p.idata(intData::num_comp_sphere)*p.rdata(realData::radius),two)));

    // Set bond components to zero
    for(int b=0; b<MAXBONDS*9; b++){ 
        p.rdata(realData::first_bond_v+b) = zero;
    }

    for(int br=0; br<MAXBRIDGES; br++){
        p.idata(intData::first_bridge+3*br) = -1;
        p.idata(intData::first_bridge+3*br+1) = -1;
        p.idata(intData::first_bridge+3*br+2) = -1;
    }

    // If nonzero MC passed in, calculate liquid volume
    if(MC > 0.0){
        p.rdata(realData::liquid_volume) = (MC > FSP) ? (p.rdata(realData::density)*p.rdata(realData::volume)/liquid_density)*(MC - FSP)/(1 - MC):0;
        p.rdata(realData::mass) = p.rdata(realData::density)*p.rdata(realData::volume) * (1.0 + MC/(1.0 - MC));
        p.rdata(realData::density) = p.rdata(realData::mass) / p.rdata(realData::volume);
    } else {
        p.rdata(realData::liquid_volume) = zero;
    }
    p.rdata(realData::total_bridge_volume) = zero;

    for(int b=0; b<MAXBONDS; b++) p.idata(intData::first_bond + b) = -1;

    for(int i=0;i<MAXSPECIES;i++)
    {
        p.rdata(realData::firstspec+i)=spec[i];
    }

    return(p);
}
