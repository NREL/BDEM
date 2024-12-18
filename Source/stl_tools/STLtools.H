#ifndef STLTOOLS_H_
#define STLTOOLS_H_

#include <STLUtils.H>

using namespace amrex;

    //- Stuff for STL grid
struct stlgrid
{
    int size[3];
    int numberOfCells;
    Real bbmax[3];
    Real bbmin[3];
    Real delta[3];
};


struct STLtools
{
    std::string name;

    AMREX_GPU_DEVICE_MANAGED Real *tri_pts=NULL;
    AMREX_GPU_DEVICE_MANAGED Real *tri_normals=NULL;
    AMREX_GPU_DEVICE_MANAGED Real *tri_mat=NULL;
    AMREX_GPU_DEVICE_MANAGED int *domndir=NULL;
    AMREX_GPU_DEVICE_MANAGED Real *rot_tri_pts=NULL;

    // std::vector<Real> tri_pts_0;
    
    AMREX_GPU_DEVICE_MANAGED int *sorted_indexarray=NULL;

    Gpu::ManagedVector<Real>* tri_pts_vec=NULL; 
    Gpu::ManagedVector<Real>* tri_normals_vec=NULL; 
    Gpu::ManagedVector<Real>* tri_mat_vec=NULL; 
    Gpu::ManagedVector<int>* domndir_vec=NULL;
    Gpu::ManagedVector<Real>* rot_tri_pts_vec=NULL;
    Gpu::ManagedVector<int>* sorted_indexarray_vec=NULL;

    //- Stress fields
    Gpu::ManagedVector<Real>* pressure_vec      = NULL;
    Gpu::ManagedVector<Real>* shear_stress_vec  = NULL;
    AMREX_GPU_DEVICE_MANAGED Real *pressure     = NULL;
    AMREX_GPU_DEVICE_MANAGED Real *shear_stress = NULL;      

    AMREX_GPU_DEVICE_MANAGED int num_tri=0;
    AMREX_GPU_DEVICE_MANAGED int ndata_per_tri=9;
    AMREX_GPU_DEVICE_MANAGED int ndata_per_normal=3;
    int nlines_per_facet=7;

    Real eigdirs[9]={0.0};
    Real bbox_lo[3]={0.0};
    Real bbox_hi[3]={0.0};
    Real expfac=2.0;

    AMREX_GPU_DEVICE_MANAGED stlgrid grid;

    //- This array has the id of the tris in each cell
    AMREX_GPU_DEVICE_MANAGED int* tris_in_grid;
    Gpu::ManagedVector<int>* tris_in_grid_vec = NULL;

    //- This array has as many elements as cells in the grid
    //  and tells where the cell data is starting in tris_in_grid
    AMREX_GPU_DEVICE_MANAGED int* cell_start;
    Gpu::ManagedVector<int>* cell_start_vec = NULL;
    AMREX_GPU_DEVICE_MANAGED int* tris_per_cell;
    Gpu::ManagedVector<int>* tris_per_cell_vec = NULL;

    
    std::map<std::pair<int,int>,Array<Real,6>> boxmap;



    void read_stl_file(std::string fname);
    void write_stl_file(std::string fname);
    void move_stl(Real timestep,int movetype,amrex::Real movedir[3],amrex::Real movecenter[3],Real movevel);
    void set_tri_mat_and_domndir(Real p[9],Real mat[9],Real rot_pts[9],int &dndir);
    void searchtriangulation(int startindex,int endindex,int *indexarray,
    Real p[3],Real &mindistsq,int &idmin,int &distcalcs) const;
    void boundingbox(int startindex,int endindex,int *indexarray,Real bx[6]);
    void update_bounding_box();
    Real getcentroid(int tri_id,int dir);
    void getSTLGeoCenter(Real center[3]);
    Real getNormalComponent(int id,Real p[3]);
    AMREX_GPU_HOST_DEVICE void brutesearch(int startindex,int endindex,int *indexarray,
                     Real p[3],Real &mindist,int &idmin) const;
    void orb_of_triangulation(int startindex,int endindex,int dir,int *indexarray);
    void sort_triangle_ids(int start,int end,int dir,int *indexarray);
    AMREX_GPU_HOST_DEVICE void brutesearch_with_intersections(int startindex,int endindex,int *indexarray,
                     Real p[3],Real &mindist,int &idmin, Real outp[3], int &num_intersects) const;
    AMREX_GPU_HOST_DEVICE void brutesearch(int startindex,int endindex,int *indexarray,
                     Real p[3],Real &mindist,int &idmin, Real dist[3]) const;

    // Uses spatial binning to get the closest triangle in the neighborhood
    AMREX_GPU_HOST_DEVICE void get_closest_local_tri
    (
        Real p[3],
        Real &mindist, 
        int& idmin,
        const Real time,
        const int movetype,
        const Real u[3],
        const Real center[3],
        const Real movevel
    ) const;
    
    void writeVTK(std::string fname) const;

    void printForces() const;

    // Builds space binning addressing of tris
    void buildGridData();
};

struct stl_specs
{
    std::string name = "";
    STLtools* stlptr = NULL;
    int dynamicstl=0;
    Vector<Real> dynstl_transl_dir;
    Real dynstl_transl_vel=0.0;
    Vector<Real> dynstl_rot_dir;
    Vector<Real> dynstl_center;
    Real dynstl_rot_vel=0.0;
    int sum_stl_forces=0;
    Real stl_vib_freq = 0.0;
    Real stl_vib_amp = 0.0;
    Vector<Real> stl_vib_dir; 
    Real stl_timeoffset = 0.0;
};


#endif
