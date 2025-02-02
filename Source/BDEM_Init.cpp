#include <BDEM_ParticleContainer.H>
#include <BDEM_Collide_Utils.H>
#include <stl_tools/STLtools.H>
#include <AMReX_Random.H>
#include <BDEM_BondedParticles.H>

void BDEMParticleContainer::InitParticles (const std::string& filename,
                                           bool &do_heat_transfer,
                                           Real temp_mean, Real temp_stdev,
                                           int contact_law, int liquid_bridging,
                                           Real liquid_density, Real MC_avg, 
                                           Real MC_stdev, Real FSP, const int solve_fibrils)
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
                                 BP_NP15, BP_NP16, BP_NP17, BP_NP18, BP_NP19, 
                                 BP_NP20, BP_NP21, BP_NP22, BP_NP23, BP_NP24, BP_NP25,
                                 BP_NP26, BP_NP27, BP_NP28};

        for (int i = 0; i < np; i++) 
        {
            ParticleType p;
            // Set id and cpu for this particle
            p.id()  = ParticleType::NextID();
            p.cpu() = ParallelDescriptor::MyProc();
            // Unique id for each particle across all procs created using id, cpu, and Cantor pairing function
            p.idata(intData::unique_id) = (int)((p.id() + p.cpu())*(p.id() + p.cpu() + 1) + p.cpu());

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

            if(solve_fibrils == 1)
            {
                ifs >> p.rdata(realData::fraction_of_fibrils);
            }
            else
            {
                p.rdata(realData::fraction_of_fibrils) = 0.;
            }
            
            p.rdata(realData::euler_angle_x) = zero;
            p.rdata(realData::euler_angle_y) = zero;
            p.rdata(realData::euler_angle_z) = zero;
            p.idata(intData::type_id) = 0;

            //set initial radius
            p.rdata(realData::radinit)=p.rdata(realData::radius);

            // Set other particle properties
            // NOTE: Calculation of particle mass for glued sphere particles assumes no component sphere overlap
            p.rdata(realData::volume)      = fourbythree*PI*pow(p.rdata(realData::radius),three);
            p.rdata(realData::mass)        = p.rdata(realData::density)*p.rdata(realData::volume);
            p.rdata(realData::Iinv)        = 2.5/(p.rdata(realData::mass)*pow(p.rdata(realData::radius),two));
            
            // User input for angular velocity when using glued sphere particles (at least for testing purposes)
            p.rdata(realData::xangvel)     = zero;
            p.rdata(realData::yangvel)     = zero;
            p.rdata(realData::zangvel)     = zero;

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
            for(int br=0; br<MAXBRIDGES; br++) p.idata(intData::first_bridge+br) = -1;

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
                                  BP_NP15, BP_NP16, BP_NP17, BP_NP18, BP_NP19, 
                                  BP_NP20, BP_NP21, BP_NP22, BP_NP23, BP_NP24, BP_NP25,
                                  BP_NP26, BP_NP27, BP_NP28};

        Real pc_pos[THREEDIM];     
        // NOTE: size of arrays should be increased if larger particles are introduced        
        int bp_ids[40];
        int bp_unique_ids[40];

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

            for(int j = 0; j<bp_types[bp_type]; j++){
                // // Unique id for each particle across all procs created using id, cpu, and Cantor pairing function
                bp_ids[j] = ParticleType::NextID();
                int cpu_id = ParallelDescriptor::MyProc();
                bp_unique_ids[j] = (int)((bp_ids[j] + cpu_id)*(bp_ids[j] + cpu_id + 1) + cpu_id);
            }
            for(int j = 0; j<bp_types[bp_type]; j++){
                ParticleType p;
                p.id() = bp_ids[j];
                p.idata(intData::unique_id) = bp_unique_ids[j];
                if(cantilever_beam_test && j == bp_types[bp_type]-1) bp_phase = -1;    // Left-most particle is held inert
                get_bonded_particle_pos(bp_type, j, bp_radius, bp_pos, bp_q, pc_pos);
                bp_init(p, bp_data, bp_phase, pc_pos, bp_radius, bp_density, bp_vel, bp_temperature, j, bp_type, bp_unique_ids, liquid_density, MC, FSP, bp_E, bp_nu);
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
        const int ls_refinement)
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
                Real ppos_inert[THREEDIM];
                ppos_inert[XDIR] = p.pos(0);
                ppos_inert[YDIR] = p.pos(1);
                ppos_inert[ZDIR] = p.pos(2);
                Real ls_value = get_levelset_value(ppos_inert, ls_refinement, phiarr, plo, dx);
                // if(ls_value < 0.0)
                if(ls_value < p.rdata(realData::radius))
                {
                    p.id()=-1;   
                }
            });
        }
    }
    Redistribute();
}

void BDEMParticleContainer::removeParticlesInsideSTL(Vector<Real> outside_point, STLtools* stlptr)
{
    const int lev = 0;
    auto& plev  = GetParticles(lev);
    const Geometry& geom = Geom(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    Array<Real,AMREX_SPACEDIM> po_arr{outside_point[0],outside_point[1],outside_point[2]};

//    #define USE_INTERSECT 

    amrex::Print()  << "STL BBs: " 
                    << stlptr->bbox_lo[0] << " " << stlptr->bbox_hi[0] << " "
                    << stlptr->bbox_lo[1] << " " << stlptr->bbox_hi[1] << " "
                    << stlptr->bbox_lo[2] << " " << stlptr->bbox_hi[2] << "\n";

    Real x_lo_bb = stlptr->bbox_lo[0];
    Real x_hi_bb = stlptr->bbox_hi[0];
    Real y_lo_bb = stlptr->bbox_lo[1];
    Real y_hi_bb = stlptr->bbox_hi[1];
    Real z_lo_bb = stlptr->bbox_lo[2];
    Real z_hi_bb = stlptr->bbox_hi[2];
    Real* tri_pts = stlptr->tri_pts;
    int num_tri = stlptr->num_tri;

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
  
            Real ploc[THREEDIM];
            ploc[XDIR] = p.pos(0);
            ploc[YDIR] = p.pos(1);
            ploc[ZDIR] = p.pos(2);
            Real ploc_t[3]={0.0};
            Real t1[3],t2[3],t3[3];
            Real outp[]={po_arr[0],po_arr[1],po_arr[2]};
            int num_intersects=0;
            int min_tri_id=0;

            if( (p.pos(0) + p.rdata(realData::radius) > x_lo_bb ) && 
                (p.pos(0) - p.rdata(realData::radius) < x_hi_bb ) &&
                (p.pos(1) + p.rdata(realData::radius) > y_lo_bb ) &&
                (p.pos(1)- p.rdata(realData::radius) < y_hi_bb ) &&
                (p.pos(2)+ p.rdata(realData::radius) > z_lo_bb ) &&
                (p.pos(2)- p.rdata(realData::radius) < z_hi_bb )
            )
            {
                Real mindist=BIGVAL;
                brutesearch_with_intersections(
                    num_tri,
                    ploc,
                    mindist,
                    min_tri_id,
                    outp,
                    num_intersects,
                    tri_pts
                );           


                // Remove particles with intersections, but also particles that are 
                // too close
                if( mindist < p.rdata(realData::radius) || num_intersects%2 == 1 )
                {
                    p.id()=-1;
                }
            }
        });
    }

    Redistribute();
}

void BDEMParticleContainer::checkParticlesInsideSTL(Vector<Real> outside_point, STLtools* stlptr)
{
    const int lev = 0;
    auto& plev  = GetParticles(lev);
    const Geometry& geom = Geom(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    Array<Real,AMREX_SPACEDIM> po_arr{outside_point[0],outside_point[1],outside_point[2]};

    amrex::Print()  << "STL BBs: " 
                    << stlptr->bbox_lo[0] << " " << stlptr->bbox_hi[0] << " "
                    << stlptr->bbox_lo[1] << " " << stlptr->bbox_hi[1] << " "
                    << stlptr->bbox_lo[2] << " " << stlptr->bbox_hi[2] << "\n";

    Real x_lo_bb = stlptr->bbox_lo[0];
    Real x_hi_bb = stlptr->bbox_hi[0];
    Real y_lo_bb = stlptr->bbox_lo[1];
    Real y_hi_bb = stlptr->bbox_hi[1];
    Real z_lo_bb = stlptr->bbox_lo[2];
    Real z_hi_bb = stlptr->bbox_hi[2];
    Real* tri_pts = stlptr->tri_pts;
    int num_tri = stlptr->num_tri;

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
  
            Real ploc[THREEDIM];
            ploc[XDIR] = p.pos(0);
            ploc[YDIR] = p.pos(1);
            ploc[ZDIR] = p.pos(2);
            Real ploc_t[3]={0.0};
            Real t1[3],t2[3],t3[3];
            Real outp[]={po_arr[0],po_arr[1],po_arr[2]};
            int num_intersects=0;
            int min_tri_id=0;

            if( (p.pos(0) + p.rdata(realData::radius) > x_lo_bb ) && 
                (p.pos(0) - p.rdata(realData::radius) < x_hi_bb ) &&
                (p.pos(1) + p.rdata(realData::radius) > y_lo_bb ) &&
                (p.pos(1)- p.rdata(realData::radius) < y_hi_bb ) &&
                (p.pos(2)+ p.rdata(realData::radius) > z_lo_bb ) &&
                (p.pos(2)- p.rdata(realData::radius) < z_hi_bb )
            )
            {
                Real mindist=BIGVAL;
                brutesearch_with_intersections(
                    num_tri,
                    ploc,
                    mindist,
                    min_tri_id,
                    outp,
                    num_intersects,
                    tri_pts
                );  
                         
                if(num_intersects%2 == 1)
                {
#ifndef AMREX_USE_GPU
                    amrex::Print()<<"particle inside stl:"<<ploc[0]<<"\t"<<ploc[1]<<"\t"<<ploc[2]<<"\n";
#else
                    amrex::Abort("particle inside stl");
#endif

                }
            }
        });
    }

}

void BDEMParticleContainer::reassignParticleProperties(Real reinit_rad, Real reinit_dens, Real reinit_E, Real reinit_nu)
{
    const int lev = 0;
    auto& plev  = GetParticles(lev);
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

        amrex::ParallelFor(np,[=]
                AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];

            // Assign new particle parameters
            p.rdata(realData::radius) = reinit_rad;
            p.rdata(realData::density) = reinit_dens;
            p.rdata(realData::E) = reinit_E;
            p.rdata(realData::nu) = reinit_nu;

            // Recalculate impacted quantities
            p.rdata(realData::radinit)=p.rdata(realData::radius);
            p.rdata(realData::volume)      = fourbythree*PI*pow(p.rdata(realData::radius),three);
            p.rdata(realData::mass)        = p.rdata(realData::density)*p.rdata(realData::volume);
            p.rdata(realData::Iinv)        = 2.5/(p.rdata(realData::mass)*pow(p.rdata(realData::radius),two));
        });
    }
}

void BDEMParticleContainer::particleColoringStriation(Real striation_len, int striation_dir)
{
    const int lev = 0;
    auto& plev  = GetParticles(lev);
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

        amrex::ParallelFor(np,[=]
                AMREX_GPU_DEVICE (int i) noexcept
        {   
            ParticleType& p = pstruct[i];

            amrex::Real rem = fmod(p.pos(striation_dir), 2.0*striation_len);
            
            if(rem < striation_len) {
                p.idata(intData::phase) = 1;
            } else {
                p.idata(intData::phase) = 2;
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
                                           int bonded_sphere_particles,
                                           Real min_rad, Real max_rad, 
                                           int p_type, Vector<int> type_list,
                                           int use_type_dist, Vector<Real> dist_list,
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
                             BP_NP15, BP_NP16, BP_NP17, BP_NP18, BP_NP19, 
                             BP_NP20, BP_NP21, BP_NP22, BP_NP23, BP_NP24, BP_NP25,
                             BP_NP26, BP_NP27, BP_NP28};

    Real pc_pos[THREEDIM];                    
    Real bp_pos[THREEDIM];
    int bp_ids[40];
    int bp_unique_ids[40];
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
                                if(use_type_dist){
                                    Real rnum = amrex::Random();
                                    Real dsum = 0.0;
                                    bool found_type = false;
                                    int type_idx = 0;
                                    while(!found_type) {
                                      dsum += dist_list[type_idx];
                                      if(rnum < dsum){
                                        found_type = true;
                                      } else{
                                        type_idx++;    
                                        if(type_idx >= type_list.size()) amrex::Abort("Something went wrong with type dist calculations!\n");
                                      }
                                    }
                                    type = type_list[type_idx];
                                } else {
                                    type = type_list[ceil(amrex::Random()*type_list.size()) - 1];
                                }
                                if ( type < 0 || type > BP_TYPES-1) Abort("\nInvalid particle type specified .\n");
                            }
                            bp_pos[XDIR] = x; bp_pos[YDIR] = y; bp_pos[ZDIR] = z;
                            // Real eax = amrex::Random()*PI/2.0;
                            // Real eay = amrex::Random()*PI/2.0;
                            // Real eaz = amrex::Random()*PI/2.0;
                            Real eax = PI/2.0;
                            // Real eay = amrex::Random()*PI/20.0;
                            Real eay = amrex::Random()*PI/2.0;
                            Real eaz = amrex::Random()*PI/20.0;
                            amrex::Real quats[4] ={cos(eaz/2.0)*cos(eay/2.0)*cos(eax/2.0) + sin(eaz/2.0)*sin(eay/2.0)*sin(eax/2.0),
                                                  -sin(eaz/2.0)*sin(eay/2.0)*cos(eax/2.0) + cos(eaz/2.0)*cos(eay/2.0)*sin(eax/2.0),
                                                   cos(eaz/2.0)*sin(eay/2.0)*cos(eax/2.0) + sin(eaz/2.0)*cos(eay/2.0)*sin(eax/2.0),
                                                   sin(eaz/2.0)*cos(eay/2.0)*cos(eax/2.0) - cos(eaz/2.0)*sin(eay/2.0)*sin(eax/2.0)};
                            Real bp_vel[THREEDIM] = {meanvel[XDIR] + fluctuation[XDIR]*(amrex::Random()-half),
                                                     meanvel[YDIR] + fluctuation[XDIR]*(amrex::Random()-half),
                                                     meanvel[ZDIR] + fluctuation[XDIR]*(amrex::Random()-half)};
                            for(int j = 0; j<p_types[type]; j++){
                                // Unique id for each particle across all procs created using id, cpu, and Cantor pairing function
                                bp_ids[j] = ParticleType::NextID();
                                int cpu_id = ParallelDescriptor::MyProc();
                                bp_unique_ids[j] = (int)((bp_ids[j] + cpu_id)*(bp_ids[j] + cpu_id + 1) + cpu_id);
                            }
                            for(int j = 0; j<p_types[type]; j++){
                                ParticleType p;
                                p.id() = bp_ids[j];
                                p.idata(intData::unique_id) = bp_unique_ids[j];
                                get_bonded_particle_pos(type, j, rad, bp_pos, quats, pc_pos);
                                bp_init(p, p_data, bp_phase, pc_pos, rad, dens, bp_vel, temp, j, type, bp_unique_ids, liquid_density, MC, FSP, E, nu);
                                host_particles.push_back(p);
                            } 
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
                                        if(use_type_dist){
                                            Real rnum = amrex::Random();
                                            Real dsum = 0.0;
                                            bool found_type = false;
                                            int type_idx = 0;
                                            while(!found_type) {
                                              dsum += dist_list[type_idx];
                                              if(rnum < dsum){
                                                found_type = true;
                                              } else{
                                                type_idx++;    
                                                if(type_idx >= type_list.size()) amrex::Abort("Something went wrong with type dist calculations!\n");
                                              }
                                            }
                                            type = type_list[type_idx];
                                        } else {
                                            type = type_list[ceil(amrex::Random()*type_list.size()) - 1];
                                        }
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
                                    for(int pi = 0; pi<p_types[type]; pi++){
                                        // Unique id for each particle across all procs created using id, cpu, and Cantor pairing function
                                        bp_ids[pi] = ParticleType::NextID();
                                        int cpu_id = ParallelDescriptor::MyProc();
                                        bp_unique_ids[pi] = (int)((bp_ids[pi] + cpu_id)*(bp_ids[pi] + cpu_id + 1) + cpu_id);
                                    }
                                    for(int pi = 0; pi<p_types[type]; pi++){
                                        ParticleType p;
                                        p.id() = bp_ids[pi];
                                        p.idata(intData::unique_id) = bp_unique_ids[pi];
                                        get_bonded_particle_pos(type, pi, rad, bp_pos, quats, pc_pos);
                                        bp_init(p, p_data, bp_phase, pc_pos, rad, dens, bp_vel, temp, pi, type, bp_unique_ids, liquid_density, MC, FSP, E, nu);
                                        host_particles.push_back(p);
                                    } 
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
                                                                             int p_type)
{
    ParticleType p;
    p.id()  = ParticleType::NextID();
    p.cpu() = ParallelDescriptor::MyProc();                
    p.idata(intData::unique_id) = (int)((p.id() + p.cpu())*(p.id() + p.cpu() + 1) + p.cpu());

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

    p.rdata(realData::euler_angle_x) = PI/2.0;                 
    p.rdata(realData::euler_angle_y) = amrex::Random()*PI/20.0;
    p.rdata(realData::euler_angle_z) = amrex::Random()*PI/20.0;

    Real eax = p.rdata(realData::euler_angle_x);
    Real eay = p.rdata(realData::euler_angle_y);
    Real eaz = p.rdata(realData::euler_angle_z);

    p.idata(intData::type_id) = p_type;

    // Set other particle properties
    p.rdata(realData::volume)      = fourbythree*PI*pow(p.rdata(realData::radius),three);
    p.rdata(realData::mass)        = p.rdata(realData::density)*p.rdata(realData::volume);
    p.rdata(realData::Iinv)        = 2.5/(p.rdata(realData::mass)*pow(p.rdata(realData::radius),two));
    p.rdata(realData::xangvel)     = zero;
    p.rdata(realData::yangvel)     = zero;
    p.rdata(realData::zangvel)     = zero;

    // Set bond components to zero
    for(int b=0; b<MAXBONDS*9; b++){ 
        p.rdata(realData::first_bond_v+b) = zero;
    }

    for(int br=0; br<MAXBRIDGES; br++) p.idata(intData::first_bridge+br) = -1;

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

    p.rdata(realData::fraction_of_fibrils) = 0.;

    return(p);
}
