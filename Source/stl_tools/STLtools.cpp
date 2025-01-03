#include<STLtools.H>

void STLtools::read_stl_file(std::string fname)
{
    std::string tmpline,tmp1,tmp2;
    num_tri = 0;

    Real p[9],mat[9],rot_pts[9];
    int dndir;

    std::ifstream infile(fname.c_str());

    Print()<<"STL file name:"<<fname<<"\n";
    while(std::getline(infile,tmpline))
    {

        if(tmpline.find("endfacet")!=std::string::npos)
        {
            num_tri++;
        }
    }

    Print()<<"number of triangles:"<<num_tri<<"\n";

    tri_pts_vec=new Gpu::ManagedVector<Real>;
    tri_normals_vec=new Gpu::ManagedVector<Real>;

    tri_pts_vec->resize(num_tri*ndata_per_tri);
    tri_normals_vec->resize(num_tri*ndata_per_normal);

    tri_mat_vec=new Gpu::ManagedVector<Real>;
    domndir_vec=new Gpu::ManagedVector<int>;
    rot_tri_pts_vec=new Gpu::ManagedVector<Real>;

    tri_mat_vec->resize(num_tri*9);
    rot_tri_pts_vec->resize(num_tri*9);
    domndir_vec->resize(num_tri);

    pressure_vec = new Gpu::ManagedVector<Real>;
    shear_stress_vec = new Gpu::ManagedVector<Real>;
    
    pressure_vec->resize(num_tri);
    shear_stress_vec->resize(num_tri*3);
    
    // Reopen file. It makes things easier
    infile.close();
    infile.open(fname.c_str());

    int tri_id = -1;
    int vertex_id = 0;

    while(std::getline(infile,tmpline))
    {

        std::istringstream iss(tmpline);
        std::string word;
        iss >> word;

        // Check if normal or vertex
        if ( word == "facet" )
        {
            tri_id++; // This is the first entry for a triangle
            iss >> word; // Ignore "normal"
            iss >> (*tri_normals_vec)[tri_id*ndata_per_normal+0];
            iss >> (*tri_normals_vec)[tri_id*ndata_per_normal+1];
            iss >> (*tri_normals_vec)[tri_id*ndata_per_normal+2];
        }
        else if ( word == "vertex" )
        {

            iss >> (*tri_pts_vec)[tri_id*ndata_per_tri + 0 + 3*vertex_id];
            iss >> (*tri_pts_vec)[tri_id*ndata_per_tri + 1 + 3*vertex_id];
            iss >> (*tri_pts_vec)[tri_id*ndata_per_tri + 2 + 3*vertex_id];

            vertex_id++;

            if (vertex_id > 2)
            {
                vertex_id = 0;
            }                    
        }

        // No other check is required
    }    

    if ( tri_id != num_tri - 1 )
    {
        amrex::Abort("Something went wrong when reading the STL file\n");
    }


    for(int i=0;i<num_tri;i++)
    {
        for(int n=0;n<9;n++)
        {
            p[n]=(*tri_pts_vec)[i*ndata_per_tri+n];
        }

        set_tri_mat_and_domndir(p,mat,rot_pts,dndir);

        for(int n=0;n<9;n++)
        {
            (*tri_mat_vec)[i*9+n]=mat[n];
            (*rot_tri_pts_vec)[i*9+n]=rot_pts[n];
        }
        (*domndir_vec)[i]=dndir;
    }

    tri_pts     = tri_pts_vec->dataPtr();
    tri_normals = tri_normals_vec->dataPtr();

    domndir = domndir_vec->dataPtr();
    tri_mat  = tri_mat_vec->dataPtr();

    rot_tri_pts=rot_tri_pts_vec->dataPtr();

    pressure = pressure_vec->dataPtr();
    shear_stress = shear_stress_vec->dataPtr();

    update_bounding_box();

    // tri_pts_0.resize(num_tri*ndata_per_tri);
    // for (int i = 0; i < tri_pts_0.size(); i++)
    // {
    //     tri_pts_0[i] = tri_pts[i];
    // }

    /*for(int i=0;i<num_tri;i++)
        {
        Print()<<"Normals:"
        <<tri_normals[i*ndata_per_normal+0]<<"\t"
        <<tri_normals[i*ndata_per_normal+1]<<"\t"
        <<tri_normals[i*ndata_per_normal+2]<<"\n";

        for(int j=0;j<3;j++)
        {
        Print()<<"point "<<j<<" :"
        <<tri_pts[i*ndata_per_tri+3*j+0]<<"\t"
        <<tri_pts[i*ndata_per_tri+3*j+1]<<"\t"
        <<tri_pts[i*ndata_per_tri+3*j+2]<<"\n";
        }
        }*/
    
    sorted_indexarray_vec=new Gpu::ManagedVector<int>;
    sorted_indexarray_vec->resize(num_tri);
    boxmap.clear();
    for(int i=0;i<STLtools::num_tri;i++)
    {
        (*sorted_indexarray_vec)[i]=i;
    }
    sorted_indexarray=sorted_indexarray_vec->dataPtr();
    
    orb_of_triangulation(0,num_tri-1,0,sorted_indexarray);
}

void STLtools::buildGridData(int gsize[3])
{
    grid.bbmax[0] = bbox_hi[0];
    grid.bbmax[1] = bbox_hi[1];
    grid.bbmax[2] = bbox_hi[2];
    grid.bbmin[0] = bbox_lo[0];
    grid.bbmin[1] = bbox_lo[1];
    grid.bbmin[2] = bbox_lo[2];
    grid.size[0] = gsize[0];
    grid.size[1] = gsize[1];
    grid.size[2] = gsize[2];
    grid.delta[0] = (bbox_hi[0] - bbox_lo[0])/grid.size[0];
    grid.delta[1] = (bbox_hi[1] - bbox_lo[1])/grid.size[1];
    grid.delta[2] = (bbox_hi[2] - bbox_lo[2])/grid.size[2];

    //- Enlarge the bounding box
    grid.bbmax[0] += 2*grid.delta[0];
    grid.bbmax[1] += 2*grid.delta[1];
    grid.bbmax[2] += 2*grid.delta[2];
    grid.bbmin[0] -= 2*grid.delta[0];
    grid.bbmin[1] -= 2*grid.delta[1];
    grid.bbmin[2] -= 2*grid.delta[2];

    grid.size[0] += 4;
    grid.size[1] += 4;
    grid.size[2] += 4;


    grid.numberOfCells = grid.size[0]*grid.size[1]*grid.size[2];

    cell_start_vec = new Gpu::ManagedVector<int>;
    cell_start_vec->resize(grid.numberOfCells);
    cell_start = cell_start_vec->dataPtr();

    tris_per_cell_vec = new Gpu::ManagedVector<int>;
    tris_per_cell_vec->resize(grid.numberOfCells);
    tris_per_cell = tris_per_cell_vec->dataPtr();

    //- Now build the tree on the CPU on one core
    std::vector<std::vector<int>> tritree(grid.numberOfCells);

    amrex::Print() << "Building stl tree for " << name;

    for (int tri = 0; tri < num_tri; tri++)
    {
        Real A[3];
        Real B[3];
        Real C[3];

        for (int coord = 0; coord < 3; coord++)
        {
            A[coord] = tri_pts[ndata_per_tri*tri+coord];
            B[coord] = tri_pts[ndata_per_tri*tri+3+coord];
            C[coord] = tri_pts[ndata_per_tri*tri+6+coord];
        }
        
        for (int cell = 0; cell < grid.numberOfCells; cell++)
        {
            Real cellCentroidDist[3];
            int cellIJK[3];

            ijkFromID(cell,grid.size,cellIJK);

            Real cellCentroidCoord[3] =
            {
                (0.5 + cellIJK[0])*grid.delta[0] + grid.bbmin[0],
                (0.5 + cellIJK[1])*grid.delta[1] + grid.bbmin[1],
                (0.5 + cellIJK[2])*grid.delta[2] + grid.bbmin[2],
            };

            int closestType;

            Real distMag = point_tri_distance(
                cellCentroidCoord,
                A,
                B,
                C,
                closestType,
                cellCentroidDist
            );

            if  (
                    std::fabs(cellCentroidDist[0]) < 0.5*grid.delta[0] + 1e-16
                    &&
                    std::fabs(cellCentroidDist[1]) < 0.5*grid.delta[1] + 1e-16
                    &&
                    std::fabs(cellCentroidDist[2]) < 0.5*grid.delta[2] + 1e-16
                )
            {
                tritree[cell].push_back(tri);
            }
        }    
    }

    // Flatten tree
    int nelems = 0;
    for (int cell = 0; cell < grid.numberOfCells; cell++)
    {
        int treesize = tritree[cell].size();
        cell_start[cell] = -1;
        tris_per_cell[cell] = treesize;

        if(treesize > 0)
        {
            cell_start[cell] = nelems;
        }

        nelems += treesize;
    }

    tris_in_grid_vec = new Gpu::ManagedVector<int>;
    tris_in_grid_vec->resize(nelems);
    tris_in_grid = tris_in_grid_vec->dataPtr();
    int local_id = 0;
    for (int cell = 0; cell < grid.numberOfCells; cell++)
    {
        for (int i = 0; i < tritree[cell].size(); i++)
        {
            tris_in_grid[local_id] = tritree[cell][i];
            local_id++;
        }
    }   

    if(local_id != nelems)
    {
        amrex::Abort("\nError when building binned data structure for stl");
    }
    // for (int i = 0; i < grid.numberOfCells; i++)
    // {
    //     Print() << "Cell " << i << " has: ";
    //     for(int j=0; j< tris_per_cell[i]; j++)
    //     {
    //         Print() << " " << tris_in_grid[cell_start[i] + j];
    //     }
    //     Print() << "\n";
    // }
    
    

    #define WRITE_STL_GRID
    #ifdef WRITE_STL_GRID
    //- Write grid
    std::ofstream gridFile((name + "_grid.vtk").c_str()); 
    gridFile << "# vtk DataFile Version 2.0\n";
    gridFile << "Simple Structured Points Example\n";
    gridFile << "ASCII\n";
    gridFile << "DATASET STRUCTURED_POINTS\n";   

    // Write dimensions
    gridFile << "DIMENSIONS " << grid.size[0] + 1 << " " << grid.size[1] + 1 << " " << grid.size[2] + 1 << "\n";

    // Write origin
    gridFile << "ORIGIN " << grid.bbmin[0] << " " << grid.bbmin[1] << " " << grid.bbmin[2] << "\n";
    
    // Write spacing
    gridFile << "SPACING " << grid.delta[0] << " " << grid.delta[1] << " " << grid.delta[2] << "\n";

    // Write point data header
    gridFile << "CELL_DATA " << grid.numberOfCells << "\n";
    
    // Define a scalar field
    gridFile << "SCALARS firstTri int 1\n";
    gridFile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < grid.numberOfCells; i++)
    {
        if ( cell_start[i] > -1)
        { 
            gridFile << tris_in_grid[cell_start[i]] << "\n";
        }
        else 
        {
            gridFile << -1 << "\n";;
        }
    }
    gridFile << "SCALARS nTri int 1\n";
    gridFile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < grid.numberOfCells; i++)
    {
        gridFile << tris_per_cell[i] << "\n";
    }

    gridFile.close();    

    #endif



    amrex::Print() << "Done building stl tree for " << name;
}

void STLtools::getSTLGeoCenter(Real center[3])
{
    center[0] = 0.;
    center[1] = 0.;
    center[2] = 0.;
    
    for(int i=0;i<STLtools::num_tri;i++)
    {
        for (int trinode = 0; trinode < 3; trinode++)
        {
            for (int dim = 0; dim < 3; dim++)
            {
                center[dim] += (*tri_pts_vec)[i*ndata_per_tri + dim + 3*trinode];
            }
        }                
    }        

    for (int dim = 0; dim < 3; dim++)
    {
        center[dim] /= ( num_tri * 3. );
    }

}

void STLtools::boundingbox(int startindex,int endindex,int *indexarray,Real bx[6])
{
    // NOTE: Returns in min/max xyz coordinates across all stl elements?
    bx[0]=1e50;
    bx[1]=1e50;
    bx[2]=1e50;

    bx[3]=-1e50;
    bx[4]=-1e50;
    bx[5]=-1e50;

    for(int i=startindex;i<=endindex;i++)
    {
        int tri_id=indexarray[i];
        for(int tri=0;tri<3;tri++)
        {
            bx[0]=std::min(bx[0],tri_pts[ndata_per_tri*tri_id+tri*3+0]);
            bx[3]=std::max(bx[3],tri_pts[ndata_per_tri*tri_id+tri*3+0]);
            
            bx[1]=std::min(bx[1],tri_pts[ndata_per_tri*tri_id+tri*3+1]);
            bx[4]=std::max(bx[4],tri_pts[ndata_per_tri*tri_id+tri*3+1]);
            
            bx[2]=std::min(bx[2],tri_pts[ndata_per_tri*tri_id+tri*3+2]);
            bx[5]=std::max(bx[5],tri_pts[ndata_per_tri*tri_id+tri*3+2]);
        }
    }
}

Real STLtools::getcentroid(int tri_id,int dir)
{
    Real bx[6];
    bx[0]=1e50;
    bx[1]=1e50;
    bx[2]=1e50;

    bx[3]=-1e50;
    bx[4]=-1e50;
    bx[5]=-1e50;

    for(int tri=0;tri<3;tri++)
    {
        bx[0]=std::min(bx[0],tri_pts[ndata_per_tri*tri_id+tri*3+0]);
        bx[3]=std::max(bx[3],tri_pts[ndata_per_tri*tri_id+tri*3+0]);

        bx[1]=std::min(bx[1],tri_pts[ndata_per_tri*tri_id+tri*3+1]);
        bx[4]=std::max(bx[4],tri_pts[ndata_per_tri*tri_id+tri*3+1]);

        bx[2]=std::min(bx[2],tri_pts[ndata_per_tri*tri_id+tri*3+2]);
        bx[5]=std::max(bx[5],tri_pts[ndata_per_tri*tri_id+tri*3+2]);
    }

    return(bx[dir]);

    /*Real cent;

        cent=0.333333*(tri_pts[tri_id*ndata_per_tri+dir+0] +
        tri_pts[tri_id*ndata_per_tri+dir+3] +
        tri_pts[tri_id*ndata_per_tri+dir+6]);

        return(cent);*/
}


void STLtools::sort_triangle_ids(int start,int end,int dir,int *indexarray)
{
    // NOTE, as coded sorts triangle IDs based on min/max xyz extents?
    
    //selection sort
    for(int it=0;it<(end-start);it++)
    {
        Real cent1=getcentroid(indexarray[start+it+1],dir);
        int i;
        for(i=start;i<=(start+it);i++)
        {
            Real cent2=getcentroid(indexarray[i],dir);
            if(cent1<=cent2)
            {
                break;
            }
        }
        for(int j=start+it+1;j>=i+1;j--)
        {
            int temp=indexarray[j];
            indexarray[j]=indexarray[j-1];
            indexarray[j-1]=temp;
        }
    }
}

void STLtools::orb_of_triangulation(int startindex,int endindex,
                            int dir,int *indexarray)
{

    int nelm=endindex-startindex+1;
    Array<Real,6> bx;

    sort_triangle_ids(startindex,endindex,dir,indexarray);
    boundingbox(startindex,endindex,indexarray,bx.data());

    std::pair<int,int> indexkey = std::make_pair(startindex, endindex);
    boxmap.insert(std::make_pair(indexkey,bx));

    if(nelm>0)
    {
        int mid=amrex::Math::floor(0.5*Real(startindex+endindex));
        orb_of_triangulation(startindex,mid-1,(dir+1)%6,indexarray);
        orb_of_triangulation(mid+1,endindex,(dir+1)%6,indexarray);
    }
}

AMREX_GPU_HOST_DEVICE void STLtools::brutesearch(int startindex,int endindex,int *indexarray,
                    Real p[3],Real &mindist,int &idmin) const
{
    mindist=1e50;
    idmin=-1;
    Real t1[3],t2[3],t3[3], distTri[3];
    for(int id=startindex;id<=endindex;id++)
    {
        int tri=indexarray[id];
        for(int dim=0;dim<3;dim++)
        {
            t1[dim]=tri_pts[ndata_per_tri*tri+dim+0];
            t2[dim]=tri_pts[ndata_per_tri*tri+dim+3];
            t3[dim]=tri_pts[ndata_per_tri*tri+dim+6];
        }

        int closestType = -1;
        Real ptdist=point_tri_distance(p,t1,t2,t3,closestType,distTri);
        if(ptdist<mindist)
        {
            mindist=ptdist;
            idmin=id;
        }
    }
}

AMREX_GPU_HOST_DEVICE void STLtools::brutesearch(int startindex,int endindex,int *indexarray,
                    Real p[3],Real &mindist,int &idmin, Real dist[3]) const
{
    mindist=1e50;
    idmin=-1;
    Real t1[3],t2[3],t3[3], distTri[3];
    for(int id=startindex;id<=endindex;id++)
    {
        int tri=indexarray[id];
        for(int dim=0;dim<3;dim++)
        {
            t1[dim]=tri_pts[ndata_per_tri*tri+dim+0];
            t2[dim]=tri_pts[ndata_per_tri*tri+dim+3];
            t3[dim]=tri_pts[ndata_per_tri*tri+dim+6];
        }

        int closestType = -1;
        Real ptdist=point_tri_distance(p,t1,t2,t3,closestType,distTri);
        //Real ptdist=point_tri_distance(p,t1,t2,t3,distTri);
        if(ptdist<mindist)
        {
            mindist=ptdist;
            idmin=id;
            copyVector(dist,distTri)
        }
    }
}

Real STLtools::getNormalComponent(int id,Real p[3])
{
    Real cent[3]={0.0};
    Real n[3]={0.0};

    n[0]=tri_normals[ndata_per_normal*id+0];
    n[1]=tri_normals[ndata_per_normal*id+1];
    n[2]=tri_normals[ndata_per_normal*id+2];

    for(int t=0;t<3;t++)
    {
        cent[0] += tri_pts[ndata_per_tri*id+t*3+0];
        cent[1] += tri_pts[ndata_per_tri*id+t*3+1];
        cent[2] += tri_pts[ndata_per_tri*id+t*3+2];
    }

    cent[0]*=0.333333;
    cent[1]*=0.333333;
    cent[2]*=0.333333;



    Real vec_dot_n=(p[0]-cent[0])*n[0]+(p[1]-cent[1])*n[1]+(p[2]-cent[2])*n[2];
    Real vecmag=std::sqrt( Distance2(p,cent) );

    Real returnval=vec_dot_n/vecmag;

    if(vecmag==0.0)
    {
        //worst rare case when p is the centroid
        returnval=1.0;
    }

    return(returnval);
}

Real boxdistsq(const Real box[6],Real p[3])
{
    Real distsq;
    Real diff[3];

    for(int i=0;i<3;i++)
    {
        diff[i]=std::max(box[i]-p[i],0.0);
        diff[i]=std::max(diff[i],p[i]-box[i+3]);
    }
    return(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);
}

void STLtools::searchtriangulation(int startindex,int endindex,int *indexarray,
                            Real p[3],Real &mindistsq,int &idmin,int &ndistcalcs) const
{
    Real bx0[6],bx1[6],bx2[6];
    Real d0sq,d1sq,d2sq,dmin2;
    Real dfacesq,dface;
    Real dmin1sq=1e50;
    Real dmin2sq=1e50;
    int mid=0;
    Real TOL=1e-10;
    Real t1[3],t2[3],t3[3];

    if(endindex>=startindex)
    {
        int localidmin;
        mid=amrex::Math::floor(0.5*(Real(startindex)+Real(endindex)));
        brutesearch(mid,mid,indexarray,p,dfacesq,localidmin);
        dfacesq*=dfacesq;
        ndistcalcs++;
        
        if(dfacesq<mindistsq)
        {
            mindistsq=dfacesq;
            idmin=localidmin;
        }

        d1sq=boxdistsq(boxmap.at(std::make_pair(startindex,mid-1)).data(),p);
        d2sq=boxdistsq(boxmap.at(std::make_pair(mid+1,endindex)).data(),p);

        if(d1sq < mindistsq+TOL)
        {
            searchtriangulation(startindex,mid-1,indexarray,p,mindistsq,idmin,ndistcalcs);
        }
        if(d2sq < mindistsq+TOL)
        {
            searchtriangulation(mid+1,endindex,indexarray,p,mindistsq,idmin,ndistcalcs);
        }
    }
}
void STLtools::set_tri_mat_and_domndir(Real p[9],Real mat[9],Real rot_pts[9],int &dndir)
{
    // Seems that mat is (essentially) a 3x3 matrix comprised of the normal vector components,
    // and two orthogonal matrix components
    // TODO: What is rot_pts?
    amrex::Real v1[3],v2[3],n[3];
    amrex::Real nhat[3],v1hat[3],v2hat[3];
    amrex::Real pr[3][3];
    int dir1,dir2,dir3;


    for(int i=0;i<3;i++)
    {
        v1[i]=p[0*3+i]-p[1*3+i];
        v2[i]=p[0*3+i]-p[2*3+i];
    }
    CrossProd(v1,v2,n);
    getunitvec(n,nhat);

    dndir=0;
    amrex::Real maxval=nhat[0];
    for(int i=0;i<3;i++)
    {
        if(nhat[i]>maxval)
        {
            maxval=nhat[i];
            dndir=i;
        }
    }
    getunitvec(v1,v1hat);
    CrossProd(nhat,v1hat,v2hat);

    dir1=dndir;
    dir2=(dndir+1)%3;
    dir3=(dndir+2)%3;

    mat[dir1*3+0]=nhat[0];
    mat[dir1*3+1]=nhat[1];
    mat[dir1*3+2]=nhat[2];

    mat[dir2*3+0]=v1hat[0];
    mat[dir2*3+1]=v1hat[1];
    mat[dir2*3+2]=v1hat[2];

    mat[dir3*3+0]=v2hat[0];
    mat[dir3*3+1]=v2hat[1];
    mat[dir3*3+2]=v2hat[2];

    for(int pt=0;pt<3;pt++)
    {
        for(int i=0;i<3;i++)
        {
            rot_pts[3*pt+i]=0.0;
            for(int j=0;j<3;j++)
            {
                rot_pts[3*pt+i] += mat[i*3+j]*p[3*pt+j];
            }
        }
    }
    
}
void STLtools::update_bounding_box()
{
    Real c_of_mass[3]={0.0};
    Real centr[3],centr_t[3];
    Real inertia_mat[9]={0.0};
    Real eigen_vectors[9];

    for(int i=0;i<num_tri;i++)
    {
        for(int dim=0;dim<3;dim++)
        {
            c_of_mass[dim] += 0.3333333*( tri_pts[ndata_per_tri*i+dim+0]
                                        + tri_pts[ndata_per_tri*i+dim+3]
                                        + tri_pts[ndata_per_tri*i+dim+6]);
        }
    }
    c_of_mass[0] /= num_tri;
    c_of_mass[1] /= num_tri;
    c_of_mass[2] /= num_tri;

    for(int i=0;i<num_tri;i++)
    {
        for(int dim=0;dim<3;dim++)
        {
            centr[dim] = 0.3333333*( tri_pts[ndata_per_tri*i+dim+0]
                                    + tri_pts[ndata_per_tri*i+dim+3]
                                    + tri_pts[ndata_per_tri*i+dim+6]);
        }

        // NOTE: why do we need inertia?
        inertia_mat[0] += (centr[0]-c_of_mass[0])*(centr[0]-c_of_mass[0]);
        inertia_mat[1] += (centr[0]-c_of_mass[0])*(centr[1]-c_of_mass[1]);
        inertia_mat[2] += (centr[0]-c_of_mass[0])*(centr[2]-c_of_mass[2]);
        inertia_mat[4] += (centr[1]-c_of_mass[1])*(centr[1]-c_of_mass[1]);
        inertia_mat[5] += (centr[1]-c_of_mass[1])*(centr[2]-c_of_mass[2]);
        inertia_mat[8] += (centr[2]-c_of_mass[2])*(centr[2]-c_of_mass[2]);
    }

    inertia_mat[3]=inertia_mat[1];
    inertia_mat[6]=inertia_mat[2];
    inertia_mat[7]=inertia_mat[5];


    // Print()<<"inertia matrix==================\n";
    // for(int i=0;i<3;i++)
    // {
    //     for(int j=0;j<3;j++)
    //     {
    //         Print()<<inertia_mat[3*j+i]<<"\t";
    //     }
    //     Print()<<"\n";
    // }
    // Print()<<"==================\n";

    geteigvectors(inertia_mat,eigen_vectors);

    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            eigdirs[3*i+j]=eigen_vectors[3*j+i]; 
        }
    }

    bbox_lo[0]=1e50;
    bbox_lo[1]=1e50;
    bbox_lo[2]=1e50;
    
    bbox_hi[0]=-1e50;
    bbox_hi[1]=-1e50;
    bbox_hi[2]=-1e50;
    
    //! This does not make any sense to me 
    //for(int i=0;i<num_tri;i++)
    // {
    //     for(int dim=0;dim<3;dim++)
    //     {
    //         centr[dim] = 0.3333333*( tri_pts[ndata_per_tri*i+dim+0]
    //                                + tri_pts[ndata_per_tri*i+dim+3]
    //                                + tri_pts[ndata_per_tri*i+dim+6]);
    //     }

    //     centr_t[0]=0.0;
    //     centr_t[1]=0.0;
    //     centr_t[2]=0.0;
    //     for(int dim=0;dim<3;dim++)
    //     {
    //         for(int j=0;j<3;j++)
    //         {
    //             centr_t[dim] += eigdirs[3*dim+j]*centr[j];
    //         }
    //     }

    //     for(int dim=0;dim<3;dim++)
    //     {
    //         if(centr_t[dim] < bbox_lo[dim])
    //         {
    //             bbox_lo[dim]=centr_t[dim];
    //         }
    //         if(centr_t[dim] > bbox_hi[dim])
    //         {
    //             bbox_hi[dim] = centr_t[dim];
    //         }
    //     }

    // }

    //! Let's compute the bb the old-fashioned way
    // Loop over the triangle and take the max/min of the nodal values
    for(int i=0;i<num_tri;i++)
    {
        for (int n = 0; n < 3; n++)
        {
            for(int dim=0;dim<3;dim++)
            {
                bbox_lo[dim] = amrex::min( tri_pts[ndata_per_tri*i+dim+n*3], bbox_lo[dim] );
                bbox_hi[dim] = amrex::max( tri_pts[ndata_per_tri*i+dim+n*3], bbox_hi[dim] );
            }
        }
        

    }

    // //expand the bounding box
    // Real len[3],bbox_cent[3];
    // for(int dim=0;dim<3;dim++)
    // {
    //     len[dim] = bbox_hi[dim]-bbox_lo[dim];
    //     bbox_cent[dim] = 0.5*(bbox_hi[dim]+bbox_lo[dim]);
    // }
    
    // for(int dim=0;dim<3;dim++)
    // {
    //     bbox_lo[dim] = bbox_cent[dim]-0.5*expfac*len[dim];
    //     bbox_hi[dim] = bbox_cent[dim]+0.5*expfac*len[dim];
    // }

    // Print()<<"inertia matrix==================\n";
    // for(int i=0;i<3;i++)
    // {
    //     for(int j=0;j<3;j++)
    //     {
    //         Print()<<inertia_mat[3*j+i]<<"\t";
    //     }
    //     Print()<<"\n";
    // }
    // Print()<<"eigen vectors==================\n";
    // for(int i=0;i<3;i++)
    // {
    //     for(int j=0;j<3;j++)
    //     {
    //         Print()<<eigen_vectors[3*j+i]<<"\t";
    //     }
    //     Print()<<"\n";
    // }
    // Print()<<"==================\n";


}

void STLtools::write_stl_file(std::string fname)
{

    std::string filename = fname + ".stl";
    
    if( fileExists(filename) ) return; // Never overwrite
    
    Real t1[3],t2[3],t3[3],n[3];
    std::ofstream outfile(filename.c_str());

    outfile<<"solid bdemsolid\n";

    for(int i=0;i<num_tri;i++)
    {
        t1[0]=tri_pts[i*ndata_per_tri+0];
        t1[1]=tri_pts[i*ndata_per_tri+1];
        t1[2]=tri_pts[i*ndata_per_tri+2];
        
        t2[0]=tri_pts[i*ndata_per_tri+3];
        t2[1]=tri_pts[i*ndata_per_tri+4];
        t2[2]=tri_pts[i*ndata_per_tri+5];
        
        t3[0]=tri_pts[i*ndata_per_tri+6];
        t3[1]=tri_pts[i*ndata_per_tri+7];
        t3[2]=tri_pts[i*ndata_per_tri+8];

        tri_n(t1,t2,t3,n);

        outfile<<"facet normal"<<n[0]<<"\t"<<n[1]<<"\t"<<n[2]<<"\n";
        outfile<<"outer loop\n";
        outfile<<"vertex "<<t1[0]<<"\t"<<t1[1]<<"\t"<<t1[2]<<"\n";
        outfile<<"vertex "<<t2[0]<<"\t"<<t2[1]<<"\t"<<t2[2]<<"\n";
        outfile<<"vertex "<<t3[0]<<"\t"<<t3[1]<<"\t"<<t3[2]<<"\n";
        outfile<<"endloop\n";
        outfile<<"endfacet\n";
    }

    outfile<<"endsolid bdemsolid";
    outfile.close();
}

void STLtools::writeVTK(std::string fname) const
{
    // We do not have a list of points in the vtk file.
    // Points will be duplicated.

    std::string filename = fname + ".vtk";
    
    if( fileExists(filename) ) return; // Never overwrite

    Real t1[3],t2[3],t3[3],n[3];
    std::ofstream vtkFile(filename.c_str());

    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Triangles with scalar and vector fields\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET UNSTRUCTURED_GRID\n";

    // Write points with duplicates
    vtkFile << "POINTS " << num_tri*3 << " float\n";
    for(int i=0;i<num_tri;i++)
    {
        vtkFile << tri_pts[i*ndata_per_tri+0] << " ";
        vtkFile << tri_pts[i*ndata_per_tri+1] << " ";
        vtkFile << tri_pts[i*ndata_per_tri+2] << "\n";
        vtkFile << tri_pts[i*ndata_per_tri+3] << " ";
        vtkFile << tri_pts[i*ndata_per_tri+4] << " ";
        vtkFile << tri_pts[i*ndata_per_tri+5] << "\n";            
        vtkFile << tri_pts[i*ndata_per_tri+6] << " ";
        vtkFile << tri_pts[i*ndata_per_tri+7] << " ";
        vtkFile << tri_pts[i*ndata_per_tri+8] << "\n";
    }

    // Write triangles
    vtkFile << "\nCELLS " << num_tri << " " << num_tri * 4 << "\n";
    for(int i=0;i<num_tri;i++)
    {
        vtkFile << "3 " << 3*i << " " << 3*i + 1 << " " << 3*i + 2 << "\n";
    }

    // Write cell types (5 corresponds to VTK_TRIANGLE)
    vtkFile << "\nCELL_TYPES " << num_tri << "\n";
    for (int i = 0; i < num_tri; i++) 
    {
        vtkFile << "5\n";
    }

    // Write scalar cell data
    vtkFile << "\nCELL_DATA " << num_tri << "\n";
    vtkFile << "SCALARS pressure float 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < num_tri; i++) 
    {
        vtkFile << pressure[i] << "\n";
    }

    // Write vector cell data
    vtkFile << "\nVECTORS shear_stress float\n";
    for (int i = 0; i < num_tri; i++) 
    {
        vtkFile << shear_stress[3*i] << " " << shear_stress[3*i + 1] << " " << shear_stress[3*i + 2] <<  "\n";
    }

    vtkFile.close();

}

void STLtools::move_stl(Real timestep,int movetype,amrex::Real movedir[3],amrex::Real movecenter[3],Real movevel)
{

    Real coord[3],newcoord[3],norm[3];
    Real P1[3],P2[3],P3[3];

    for(int i=0;i<num_tri;i++)
    {
        for(int tr=0;tr<3;tr++)
        {
            coord[0]=tri_pts[i*ndata_per_tri+tr*3+0];
            coord[1]=tri_pts[i*ndata_per_tri+tr*3+1];
            coord[2]=tri_pts[i*ndata_per_tri+tr*3+2];

            move_this_point(coord,newcoord,timestep,movetype,movedir,movecenter,movevel);
            
            tri_pts[i*ndata_per_tri+tr*3+0]=newcoord[0];
            tri_pts[i*ndata_per_tri+tr*3+1]=newcoord[1];
            tri_pts[i*ndata_per_tri+tr*3+2]=newcoord[2];

        }

        for(int dim=0;dim<3;dim++)
        {
            P1[dim]=tri_pts[i*ndata_per_tri+0+dim];
            P2[dim]=tri_pts[i*ndata_per_tri+3+dim];
            P3[dim]=tri_pts[i*ndata_per_tri+6+dim];
        }
        
        tri_n(P1,P2,P3,norm);

        tri_normals[i*ndata_per_normal+0]=norm[0];
        tri_normals[i*ndata_per_normal+1]=norm[1];
        tri_normals[i*ndata_per_normal+2]=norm[2];
    }

  //  move_this_point(grid.bbmin,newcoord,timestep,movetype,movedir,movecenter,movevel);
  //  move_this_point(grid.bbmax,newcoord,timestep,movetype,movedir,movecenter,movevel);
}

void STLtools::printForces() const
{
    Real ptot = 0.;
    Real sstot[3] = {0., 0., 0.};

    for (int i = 0; i < num_tri; i++)
    {
        ptot += pressure[i];
        sstot[0] += shear_stress[3*i];
        sstot[1] += shear_stress[3*i + 1];
        sstot[2] += shear_stress[3*i + 2];
    }
    
    Print() << "pressure = " << ptot << "   shear stress = ( " << sstot[0] << " " << sstot[1] << " " << sstot[2] << " )\n";

}



