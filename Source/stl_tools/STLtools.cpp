#include<STLtools.H>
#include<MoveUtils.H>

namespace STLtools
{
    AMREX_GPU_DEVICE_MANAGED Real *tri_pts=NULL;
    AMREX_GPU_DEVICE_MANAGED Real *tri_normals=NULL;
    AMREX_GPU_DEVICE_MANAGED Real *tri_mat=NULL;
    AMREX_GPU_DEVICE_MANAGED int *domndir=NULL;
    AMREX_GPU_DEVICE_MANAGED Real *rot_tri_pts=NULL;
    
    AMREX_GPU_DEVICE_MANAGED int *sorted_indexarray=NULL;

    Gpu::ManagedVector<Real>* tri_pts_vec=NULL; 
    Gpu::ManagedVector<Real>* tri_normals_vec=NULL; 
    Gpu::ManagedVector<Real>* tri_mat_vec=NULL; 
    Gpu::ManagedVector<int>* domndir_vec=NULL;
    Gpu::ManagedVector<Real>* rot_tri_pts_vec=NULL;
    Gpu::ManagedVector<int>* sorted_indexarray_vec=NULL;

    AMREX_GPU_DEVICE_MANAGED int num_tri=0;
    AMREX_GPU_DEVICE_MANAGED int ndata_per_tri=9;
    AMREX_GPU_DEVICE_MANAGED int ndata_per_normal=3;
    int nlines_per_facet=7;

    Real eigdirs[9]={0.0};
    Real bbox_lo[3]={0.0};
    Real bbox_hi[3]={0.0};
    Real expfac=2.0;
    
    std::map<std::pair<int,int>,Array<Real,6>> boxmap;

    void read_stl_file(std::string fname)
    {
        std::string tmpline,tmp1,tmp2;
        int nlines=0;

        Real p[9],mat[9],rot_pts[9];
        int dndir;

        std::ifstream infile(fname.c_str());
        Print()<<"STL file name:"<<fname<<"\n";

        std::getline(infile,tmpline); //solid <solidname>
        while(!infile.eof())
        {
            std::getline(infile,tmpline);
            if(tmpline.find("endsolid")!=std::string::npos)
            {
                break;
            }
            nlines++;
        }

        if(nlines%nlines_per_facet!=0)
        {
            amrex::Abort("may be there are blank lines in the STL file\n");
        }

        num_tri=nlines/nlines_per_facet;
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

        infile.seekg(0);
        std::getline(infile,tmpline); //solid <solidname>

        for(int i=0;i<num_tri;i++)
        {
            std::getline(infile,tmpline);  //facet normal
            std::istringstream fcnormal(tmpline);
            fcnormal>>tmp1>>tmp2
                >>(*tri_normals_vec)[i*ndata_per_normal+0]
                >>(*tri_normals_vec)[i*ndata_per_normal+1]
                >>(*tri_normals_vec)[i*ndata_per_normal+2];

            std::getline(infile,tmpline); // outer loop

            std::getline(infile,tmpline); //vertex 1
            std::istringstream vertex1(tmpline);
            vertex1>>tmp1
                >>(*tri_pts_vec)[i*ndata_per_tri+0]
                >>(*tri_pts_vec)[i*ndata_per_tri+1]
                >>(*tri_pts_vec)[i*ndata_per_tri+2];

            std::getline(infile,tmpline); //vertex 2
            std::istringstream vertex2(tmpline);
            vertex2>>tmp1 
                >>(*tri_pts_vec)[i*ndata_per_tri+3]
                >>(*tri_pts_vec)[i*ndata_per_tri+4]
                >>(*tri_pts_vec)[i*ndata_per_tri+5];

            std::getline(infile,tmpline); //vertex 3
            std::istringstream vertex3(tmpline);
            vertex3>>tmp1 //vertex
                >>(*tri_pts_vec)[i*ndata_per_tri+6]
                >>(*tri_pts_vec)[i*ndata_per_tri+7]
                >>(*tri_pts_vec)[i*ndata_per_tri+8];

            std::getline(infile,tmpline); //end loop
            std::getline(infile,tmpline); //end facet

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

        update_bounding_box();

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
   
   void boundingbox(int startindex,int endindex,int *indexarray,Real bx[6])
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

    Real getcentroid(int tri_id,int dir)
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


    void sort_triangle_ids(int start,int end,int dir,int *indexarray)
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

    void orb_of_triangulation(int startindex,int endindex,
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

    AMREX_GPU_HOST_DEVICE void brutesearch(int startindex,int endindex,int *indexarray,
                     Real p[3],Real &mindist,int &idmin)
    {
        mindist=1e50;
        idmin=-1;
        Real t1[3],t2[3],t3[3];
        for(int id=startindex;id<=endindex;id++)
        {
            int tri=indexarray[id];
            for(int dim=0;dim<3;dim++)
            {
                t1[dim]=tri_pts[ndata_per_tri*tri+dim+0];
                t2[dim]=tri_pts[ndata_per_tri*tri+dim+3];
                t3[dim]=tri_pts[ndata_per_tri*tri+dim+6];
            }
            Real ptdist=point_tri_distance(p,t1,t2,t3);
            if(ptdist<mindist)
            {
                mindist=ptdist;
                idmin=id;
            }
        }
    }

    Real getNormalComponent(int id,Real p[3])
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

    Real boxdistsq(Real box[6],Real p[3])
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

    void searchtriangulation(int startindex,int endindex,int *indexarray,
                             Real p[3],Real &mindistsq,int &idmin,int &ndistcalcs)
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
    void set_tri_mat_and_domndir(Real p[9],Real mat[9],Real rot_pts[9],int &dndir)
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
    void update_bounding_box()
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
        
        Print()<<"inertia matrix==================\n";
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                Print()<<inertia_mat[3*j+i]<<"\t";
            }
            Print()<<"\n";
        }
        Print()<<"==================\n";

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
        for(int i=0;i<num_tri;i++)
        {
            for(int dim=0;dim<3;dim++)
            {
                centr[dim] = 0.3333333*( tri_pts[ndata_per_tri*i+dim+0]
                                       + tri_pts[ndata_per_tri*i+dim+3]
                                       + tri_pts[ndata_per_tri*i+dim+6]);
            }

            centr_t[0]=0.0;
            centr_t[1]=0.0;
            centr_t[2]=0.0;
            for(int dim=0;dim<3;dim++)
            {
                for(int j=0;j<3;j++)
                {
                    centr_t[dim] += eigdirs[3*dim+j]*centr[j];
                }
            }

            for(int dim=0;dim<3;dim++)
            {
                if(centr_t[dim] < bbox_lo[dim])
                {
                    bbox_lo[dim]=centr_t[dim];
                }
                if(centr_t[dim] > bbox_hi[dim])
                {
                    bbox_hi[dim] = centr_t[dim];
                }
            }

        }

        //expand the bounding box
        Real len[3],bbox_cent[3];
        for(int dim=0;dim<3;dim++)
        {
            len[dim] = bbox_hi[dim]-bbox_lo[dim];
            bbox_cent[dim] = 0.5*(bbox_hi[dim]+bbox_lo[dim]);
        }
        
        for(int dim=0;dim<3;dim++)
        {
            bbox_lo[dim] = bbox_cent[dim]-0.5*expfac*len[dim];
            bbox_hi[dim] = bbox_cent[dim]+0.5*expfac*len[dim];
        }

        Print()<<"inertia matrix==================\n";
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                Print()<<inertia_mat[3*j+i]<<"\t";
            }
            Print()<<"\n";
        }
        Print()<<"eigen vectors==================\n";
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                Print()<<eigen_vectors[3*j+i]<<"\t";
            }
            Print()<<"\n";
        }
        Print()<<"==================\n";


    }
    
    void write_stl_file(std::string fname)
    {
        Real t1[3],t2[3],t3[3],n[3];
        std::ofstream outfile(fname.c_str());

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
    }

    void move_stl(Real timestep,int movetype,int movedir,amrex::Real movecenter[3],Real movevel)
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
    }
}
