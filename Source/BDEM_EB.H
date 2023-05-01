#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <BDEM_CheckPair.H>
#include <BDEM_ParticleContainer.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Polynomial.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Translation.H>

#include <AMReX_EB_utils.H>

namespace EBtools
{
    extern EBFArrayBoxFactory* ebfactory;
    extern MultiFab *lsphi;
    extern int ls_refinement;
    extern bool using_levelset_geometry;
    void init_eb(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm);
    void make_cylinder_levelset(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm);
    void make_silo_levelset(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm);
    void make_hopper_levelset(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm);
    void init_eb(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm);
    void make_wedge_hopper_levelset(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm);
    void make_ring_levelset(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm);
}

