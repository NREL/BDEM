#include <BDEM_EB.H>

namespace EBtools
{
    EBFArrayBoxFactory* ebfactory=NULL;
    MultiFab* lsphi=NULL;
    int ls_refinement=1;
    bool using_levelset_geometry=false;

    void make_cylinder_levelset(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm)
    {
        int ls_ref = ls_refinement;

        // Define nGrow of level-set and EB
        int nghost = 1;

        bool inside        = true;
        amrex::Real radius = 0.0002;
        amrex::Real height = -1.;
        int direction      = 0;
        Vector<amrex::Real> centervec(3);
        Vector<amrex::Real> blovec(3);
        Vector<amrex::Real> bhivec(3);

        amrex::ParmParse pp("cylinder");
        pp.query("dem_inside", inside);
        pp.get("radius",     radius);
        pp.get("height",     height);
        pp.get("direction",  direction);
        pp.getarr("center",    centervec,  0, 3);
        Array<amrex::Real,3> center={centervec[0], centervec[1], centervec[2]};

        /****************************************************************************
         *                                                                          *
         * Build standard EB Factories                                              *
         *                                                                          *
         ***************************************************************************/

        amrex::Print() << " " << std::endl;
        amrex::Print() << " Internal Flow: " << inside << std::endl;
        amrex::Print() << " Radius:    " << radius    << std::endl;
        amrex::Print() << " Height:    " << height    << std::endl;
        amrex::Print() << " Direction: " << direction << std::endl;
        amrex::Print() << " Center:    " << center[0] << ", "
        << center[1] << ", "
        << center[2] << std::endl;

        //Define EB
        EB2::CylinderIF cylinder(radius, height, direction, center, inside);
        auto cyl_gshop = EB2::makeShop(cylinder);

        //make domain finer for levelset
        Box dom_ls = geom.Domain();
        dom_ls.refine(ls_ref);
        Geometry geom_ls(dom_ls);

        int required_coarsening_level = 0;
        int max_coarsening_level=10;
        if (ls_refinement > 1) 
        {
            int tmp = ls_refinement;
            while (tmp >>= 1) ++required_coarsening_level;
        }

        // Build EB
        EB2::Build(cyl_gshop, geom_ls, required_coarsening_level, max_coarsening_level);

        const EB2::IndexSpace & ebis   = EB2::IndexSpace::top();
        const EB2::Level &      eblev  = ebis.getLevel(geom);
        //create lslev
        const EB2::Level & lslev = ebis.getLevel(geom_ls);

        //build factory
        ebfactory = new EBFArrayBoxFactory(eblev, geom, ba, dm,
                                           {nghost, nghost,
                                               nghost}, EBSupport::full);

        //create nodal multifab with level-set refinement
        BoxArray ls_ba = amrex::convert(ba, IntVect::TheNodeVector());
        ls_ba.refine(ls_ref);
        lsphi = new MultiFab;
        lsphi->define(ls_ba, dm, 1, nghost);

        //call signed distance
        amrex::FillSignedDistance (*lsphi,lslev,*ebfactory,ls_ref);

    }

    void make_hopper_levelset(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm)
    {
        int ls_ref = ls_refinement;
        // Define nGrow of level-set and EB
        int nghost = 1;

        amrex::Real exit_radius     = 0.0002;
        amrex::Real cyl_radius      = 0.0002;
        amrex::Real frustum_height  = 0.0002;
        Vector<amrex::Real> centervec(3);

        amrex::ParmParse pp("hopper");
        pp.get("exit_radius",     exit_radius);
        pp.get("cyl_radius",      cyl_radius);
        pp.get("frustum_height",  frustum_height);

        const auto plo = geom.ProbLoArray();
        const auto phi = geom.ProbHiArray();


        Array<amrex::Real,3> frustum_point ={exit_radius,0.0,0.0};
        Array<amrex::Real,3> frustum_normal={frustum_height,(exit_radius-cyl_radius),0.0};

        EB2::PlaneIF frustum(frustum_point, frustum_normal);

        Array<amrex::Real,3> cyl_point={cyl_radius,frustum_height,0.0};
        Array<amrex::Real,3> cyl_normal={1.0,0.0,0.0};
        EB2::PlaneIF cyl(cyl_point, cyl_normal);

        Array<Real,3> center={0.5*(plo[0]+phi[0]),0.5*(plo[1]+phi[1]),0.0};
        auto hopper = EB2::translate(EB2::lathe(EB2::makeUnion(cyl, frustum)),center);

        amrex::Print() << " " << std::endl;
        amrex::Print() << " Exit Radius:        " << exit_radius    << std::endl;
        amrex::Print() << " Cylinder Radius:    " << cyl_radius     << std::endl;
        amrex::Print() << " Frustum height      " << frustum_height << std::endl;

        //Define EB
        auto hopper_gshop = EB2::makeShop(hopper);

        //make domain finer for levelset
        Box dom_ls = geom.Domain();
        dom_ls.refine(ls_ref);
        Geometry geom_ls(dom_ls);

        int required_coarsening_level = 0;
        int max_coarsening_level=10;
        if (ls_refinement > 1) 
        {
            int tmp = ls_refinement;
            while (tmp >>= 1) ++required_coarsening_level;
        }

        // Build EB
        EB2::Build(hopper_gshop, geom_ls, 
                   required_coarsening_level, max_coarsening_level);

        const EB2::IndexSpace & ebis   = EB2::IndexSpace::top();
        const EB2::Level &      eblev  = ebis.getLevel(geom);
        //create lslev
        const EB2::Level & lslev = ebis.getLevel(geom_ls);

        //build factory
        ebfactory = new EBFArrayBoxFactory(eblev, geom, ba, dm,
                                           {nghost, nghost,
                                               nghost}, EBSupport::full);

        //create nodal multifab with level-set refinement
        BoxArray ls_ba = amrex::convert(ba, IntVect::TheNodeVector());
        ls_ba.refine(ls_ref);
        lsphi = new MultiFab;
        lsphi->define(ls_ba, dm, 1, nghost);

        amrex::FillSignedDistance (*lsphi,lslev,*ebfactory,ls_ref);

    }
    void make_wedge_hopper_levelset(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm)
    {
        int ls_ref = ls_refinement;
        // Define nGrow of level-set and EB
        int nghost = 1;

        amrex::Real exit_size     = 0.0002;
        amrex::Real bin_size      = 0.0002;
        amrex::Real funnel_height  = 0.0002;
        Vector<amrex::Real> centervec(3);

        amrex::ParmParse pp("wedge_hopper");
        pp.get("exit_size",     exit_size);
        pp.get("bin_size",      bin_size);
        pp.get("funnel_height",  funnel_height);

        const auto plo = geom.ProbLoArray();
        const auto phi = geom.ProbHiArray();


        Array<amrex::Real,3> funnel_point1 ={0.5*exit_size,0.0,0.0};
        Array<amrex::Real,3> funnel_normal1={funnel_height,0.5*(exit_size-bin_size),0.0};
        EB2::PlaneIF funnel1(funnel_point1, funnel_normal1);
        Array<amrex::Real,3> bin_point1={0.5*bin_size,funnel_height,0.0};
        Array<amrex::Real,3> bin_normal1={1.0,0.0,0.0};
        EB2::PlaneIF bin1(bin_point1, bin_normal1);

        Array<amrex::Real,3> funnel_point2 ={-0.5*exit_size,0.0,0.0};
        Array<amrex::Real,3> funnel_normal2={-funnel_height,0.5*(exit_size-bin_size),0.0};
        EB2::PlaneIF funnel2(funnel_point2, funnel_normal2);
        Array<amrex::Real,3> bin_point2={-0.5*bin_size,funnel_height,0.0};
        Array<amrex::Real,3> bin_normal2={-1.0,0.0,0.0};
        EB2::PlaneIF bin2(bin_point2, bin_normal2);

        Array<Real,3> center={0.5*(plo[0]+phi[0]),plo[1],0.5*(plo[2]+phi[2])};
        auto hopper = EB2::translate(EB2::makeUnion(funnel1,bin1,funnel2,bin2),center);

        //Define EB
        auto hopper_gshop = EB2::makeShop(hopper);

        //make domain finer for levelset
        Box dom_ls = geom.Domain();
        dom_ls.refine(ls_ref);
        Geometry geom_ls(dom_ls);

        int required_coarsening_level = 0;
        int max_coarsening_level=10;
        if (ls_refinement > 1) 
        {
            int tmp = ls_refinement;
            while (tmp >>= 1) ++required_coarsening_level;
        }

        // Build EB
        EB2::Build(hopper_gshop, geom_ls, 
                   required_coarsening_level, max_coarsening_level);

        const EB2::IndexSpace & ebis   = EB2::IndexSpace::top();
        const EB2::Level &      eblev  = ebis.getLevel(geom);
        //create lslev
        const EB2::Level & lslev = ebis.getLevel(geom_ls);

        //build factory
        ebfactory = new EBFArrayBoxFactory(eblev, geom, ba, dm,
                                           {nghost, nghost,
                                               nghost}, EBSupport::full);

        //create nodal multifab with level-set refinement
        BoxArray ls_ba = amrex::convert(ba, IntVect::TheNodeVector());
        ls_ba.refine(ls_ref);
        lsphi = new MultiFab;
        lsphi->define(ls_ba, dm, 1, nghost);

        amrex::FillSignedDistance (*lsphi,lslev,*ebfactory,ls_ref);

    }
    void make_silo_levelset(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm)
    {
        int ls_ref = ls_refinement;
        // Define nGrow of level-set and EB
        int nghost = 1;

        amrex::Real radius = 0.0002;
        amrex::Real height = -1.;
        amrex::Real channel_width  = -1.0; //along (dir+2)%3 (x)
        amrex::Real channel_depth  = -1.0; //along (dir+1)%3 (z)
        amrex::Real channel_offset  =  0.0; //along (dir+2)%3 (x)
        int direction      = 0;
        Vector<amrex::Real> centervec(3);

        amrex::ParmParse pp("silo");
        pp.get("radius",     radius);
        pp.get("height",     height);
        pp.get("direction",  direction);
        pp.getarr("center",    centervec,  0, 3);
        pp.get("channel_width",   channel_width);
        pp.get("channel_depth",   channel_depth);
        pp.get("channel_offset",   channel_offset);

        const auto plo = geom.ProbLoArray();
        const auto phi = geom.ProbHiArray();

        int dir1=(direction+1)%3;
        int dir2=(direction+2)%3;

        Array<amrex::Real,3> center={centervec[0], centervec[1], centervec[2]};

        Array<amrex::Real,3> blo;
        Array<amrex::Real,3> bhi;

        blo[direction] = plo[direction]-(phi[direction]-plo[direction]); //y
        blo[dir2]      = center[dir2]-channel_offset-0.5*channel_width;   //x
        blo[dir1]      = center[dir1]-0.5*channel_depth;                 //z

        bhi[direction] = center[direction];                              //y
        bhi[dir2]      = center[dir2]-channel_offset+0.5*channel_width;   //x
        bhi[dir1]      = center[dir1]+0.5*channel_depth;                 //z

        amrex::Print() << " " << std::endl;
        amrex::Print() << " Radius:    " << radius    << std::endl;
        amrex::Print() << " Height:    " << height    << std::endl;
        amrex::Print() << " Direction: " << direction << std::endl;
        amrex::Print() << " Center:    " << center[0] << ", "
        << center[1] << ", "
        << center[2] << std::endl;
        amrex::Print() << " blo:    " << blo[0] << ", "
        << blo[1] << ", "
        << blo[2] << std::endl;
        amrex::Print() << " bhi:    " << bhi[0] << ", "
        << bhi[1] << ", "
        << bhi[2] << std::endl;

        //Define EB
        EB2::CylinderIF cylinder(radius, height, direction, center, false);
        EB2::BoxIF channel(blo, bhi, false);
        auto slotcyl = EB2::makeComplement(EB2::makeUnion(cylinder, channel));

        auto slotcyl_gshop = EB2::makeShop(slotcyl);

        //make domain finer for levelset
        Box dom_ls = geom.Domain();
        dom_ls.refine(ls_ref);
        Geometry geom_ls(dom_ls);

        int required_coarsening_level = 0;
        int max_coarsening_level=10;
        if (ls_refinement > 1) 
        {
            int tmp = ls_refinement;
            while (tmp >>= 1) ++required_coarsening_level;
        }

        // Build EB
        EB2::Build(slotcyl_gshop, geom_ls, 
                   required_coarsening_level, max_coarsening_level);

        const EB2::IndexSpace & ebis   = EB2::IndexSpace::top();
        const EB2::Level &      eblev  = ebis.getLevel(geom);
        //create lslev
        const EB2::Level & lslev = ebis.getLevel(geom_ls);

        //build factory
        ebfactory = new EBFArrayBoxFactory(eblev, geom, ba, dm,
                                           {nghost, nghost,
                                               nghost}, EBSupport::full);

        //create nodal multifab with level-set refinement
        BoxArray ls_ba = amrex::convert(ba, IntVect::TheNodeVector());
        ls_ba.refine(ls_ref);
        lsphi = new MultiFab;
        lsphi->define(ls_ba, dm, 1, nghost);

        amrex::FillSignedDistance (*lsphi,lslev,*ebfactory,ls_ref);
    }

    void init_eb(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm)
    {
        std::string geom_kind="all_regular"; 
        amrex::ParmParse pp("bdem");

        pp.query("kind_of_geometry",geom_kind);
        pp.query("refine_level_set",ls_refinement);

        if(geom_kind!="all_regular")
        {
            if(geom_kind=="cylinder")
            {
                using_levelset_geometry=true;
                make_cylinder_levelset(geom,ba,dm);
            }
            else if(geom_kind=="silo")
            {
                using_levelset_geometry=true;
                make_silo_levelset(geom,ba,dm);
            }
            else if(geom_kind=="hopper")
            {
                using_levelset_geometry=true;
                make_hopper_levelset(geom,ba,dm);
            }
            else if(geom_kind=="wedge_hopper")
            {
                using_levelset_geometry=true;
                make_wedge_hopper_levelset(geom,ba,dm);
            }
            else if(geom_kind=="eb2")
            {
                using_levelset_geometry=true;
                int required_coarsening_level = 0;
                int max_coarsening_level=50;
                int nghost = 1;
                if (ls_refinement > 1) 
                {
                    int tmp = ls_refinement;
                    while (tmp >>= 1) ++required_coarsening_level;
                }
                Box dom_ls = geom.Domain();
                dom_ls.refine(ls_refinement);
                Geometry geom_ls(dom_ls);
                amrex::EB2::Build(geom_ls, required_coarsening_level, max_coarsening_level);

                const EB2::IndexSpace & ebis   = EB2::IndexSpace::top();
                const EB2::Level &      eblev  = ebis.getLevel(geom);
                //create lslev
                const EB2::Level & lslev = ebis.getLevel(geom_ls);

                //build factory
                ebfactory = new EBFArrayBoxFactory(eblev, geom, ba, dm,
                                                   {nghost, nghost,nghost}, 
                                                   EBSupport::full);

                //create nodal multifab with level-set refinement
                BoxArray ls_ba = amrex::convert(ba, IntVect::TheNodeVector());
                ls_ba.refine(ls_refinement);
                lsphi = new MultiFab;
                lsphi->define(ls_ba, dm, 1, nghost);

                //call signed distance
                amrex::FillSignedDistance (*lsphi,lslev,*ebfactory,ls_refinement);
            }
            else
            {
                amrex::Abort("embedded geometry not implemented yet\n");   
            }
        }

        if(using_levelset_geometry)
        {
            lsphi->FillBoundary(geom.periodicity()); 
            const std::string& pltfile = "ebplt";
            Box dom_ls = geom.Domain();
            dom_ls.refine(ls_refinement);
            Geometry geom_ls(dom_ls);
            BoxArray plot_ba=ba;
            plot_ba.refine(ls_refinement); 
            MultiFab plotmf(plot_ba,dm,lsphi->nComp(),0);
            amrex::average_node_to_cellcenter(plotmf, 0, *lsphi, 0, lsphi->nComp());
            WriteSingleLevelPlotfile(pltfile, plotmf, {"phi"}, geom_ls, 0.0, 0);
        }
    }
}
