#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EB_utils.H>
#include <STLtools.H>

using namespace amrex;

void test_distance()
{
    amrex::Real t1[3],t2[3],t3[3],mat[9],rot_pt[9],p[9];
    amrex::Real t0[3];
    int dndir;

    std::cout<<"enter tri point1:\n";
    std::cin>>t1[0]>>t1[1]>>t1[2];

    std::cout<<"enter tri point2:\n";
    std::cin>>t2[0]>>t2[1]>>t2[2];

    std::cout<<"enter tri point3:\n";
    std::cin>>t3[0]>>t3[1]>>t3[2];

    std::cout<<"enter search point:\n";
    std::cin>>t0[0]>>t0[1]>>t0[2];

    for(int i=0;i<3;i++)
    {
        p[i+0]=t1[i];
        p[i+3]=t2[i];
        p[i+6]=t3[i];
    }

    STLtools::set_tri_mat_and_domndir(p,mat,rot_pt,dndir);

    amrex::Print()<<"dndir:"<<dndir<<"\n";

    const Real strt1 = amrex::second();
    Real dist1,dist2;
    for(int iter=0;iter<1000000;iter++)
    { 
        dist1=STLtools::point_tri_distance(t0,rot_pt,
                dndir,mat);
    }
    Real time1 = amrex::second() - strt1;

    const Real strt2 = amrex::second();
    for(int iter=0;iter<1000000;iter++)
    {
        dist2=STLtools::point_tri_distance(t0,t1,t2,t3);
    }
    Real time2 = amrex::second() - strt2;

    amrex::Print()<<"distance:"<<dist1<<"\t"<<dist2<<"\t"<<time1<<"\t"<<time2<<"\n";
}

void test_intersection()
{
    //check intersection and tri_n
    Real t1[3],t2[3],t3[3],v1[3],v2[3],ip[3],n[3];

    t1[0]=0.0; t1[1]=0.0; t1[2]=0.0;
    t2[0]=0.0; t2[1]=1.0; t2[2]=0.0;
    t3[0]=0.0; t3[0]=0.0; t3[2]=1.0;

    v1[0]=0.01; v1[1]=0.0; v1[2]=0.0;
    v2[0]=-0.02; v2[1]=0.0; v2[2]=0.0;

    STLtools::find_intersection_point(v1,v2,t1,t2,t3,ip);

    Print()<<"ip:"<<ip[0]<<"\t"<<ip[1]<<"\t"<<ip[2]<<"\n";

    STLtools::tri_n(t1,t2,t3,v1,n);
    Print()<<"n:"<<n[0]<<"\t"<<n[1]<<"\t"<<n[2]<<"\n";

    STLtools::tri_n(t1,t2,t3,v2,n);
    Print()<<"n:"<<n[0]<<"\t"<<n[1]<<"\t"<<n[2]<<"\n";

}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
        int nghost = 1;
        int max_grid_size=8;
        MultiFab lsphi;
        std::string stl_fname;

        Vector<Real> plo;
        Vector<Real> phi;
        Vector<int> ncells;
        Vector<Real> pointoutside;
        Real dx[3];

        //amrex::Print()<<"Distance test=========\n";
        //test_distance();
        //amrex::Print()<<"======================\n";
        amrex::Print()<<"Intersection test=====\n";
        test_intersection();
        amrex::Print()<<"======================\n";

        ParmParse pp;
        pp.getarr("prob_lo",plo);
        pp.getarr("prob_hi",phi);
        pp.getarr("ncells",ncells);
        pp.get("stl_file",stl_fname);
        pp.getarr("outside_point",pointoutside);

        RealBox real_box({AMREX_D_DECL(plo[0], plo[1], plo[2])},
                {AMREX_D_DECL(phi[0], phi[1], phi[2])});

        Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};

        IntVect domain_lo(AMREX_D_DECL(0,0,0));
        IntVect domain_hi(AMREX_D_DECL(ncells[0],ncells[1],ncells[2]));

        dx[0]=(phi[0]-plo[0])/ncells[0];
        dx[1]=(phi[1]-plo[1])/ncells[1];
        dx[2]=(phi[2]-plo[2])/ncells[2];

        Box domain(domain_lo, domain_hi);
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);

        Geometry geom(domain,real_box,CoordSys::cartesian,is_periodic);
        DistributionMapping dm(ba);

        BoxArray nodal_ba = amrex::convert(ba, IntVect::TheNodeVector());
        lsphi.define(nodal_ba, dm, 1, nghost);

        STLtools::read_stl_file(stl_fname);

        Array<Real,AMREX_SPACEDIM> plo_arr{plo[0],plo[1],plo[2]};
        Array<Real,AMREX_SPACEDIM> len_arr{phi[0]-plo[0],phi[1]-plo[1],phi[2]-plo[2]};
        Array<Real,AMREX_SPACEDIM> po_arr{pointoutside[0],pointoutside[1],pointoutside[2]};
       
        //check orb based search
        //=========================
        Real mindist=1e50;
        int idmin=-1;
        int ndistcalcs; 
        STLtools::searchtriangulation(0,STLtools::num_tri-1,STLtools::sorted_indexarray,
                                      po_arr.data(),mindist,idmin,ndistcalcs);
        amrex::Print()<<"mindist,id,triid:"<<std::sqrt(mindist)<<
        "\t"<<idmin<<"\t"<<STLtools::sorted_indexarray[idmin]<<"\n"; 

        //brute search
        //=========================
        mindist=1e50;
        idmin=-1;
        STLtools::brutesearch(0,STLtools::num_tri-1,STLtools::sorted_indexarray,
                              po_arr.data(),mindist,idmin);
        amrex::Print()<<"brute force mindist,id,triid:"<<mindist<<
        "\t"<<idmin<<"\t"<<STLtools::sorted_indexarray[idmin]<<"\n"; 
        //=========================

        for (MFIter mfi(lsphi); mfi.isValid(); ++mfi) // Loop over grids
        {
            const Box& bx = mfi.validbox();

            auto lsphi_arr=lsphi[mfi].array();

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Real coords[3],centpvec[3];
                Real coords_t[3]={0.0};
                Real t1[3],t2[3],t3[3];
                Real outp[]={po_arr[0],po_arr[1],po_arr[2]};
                int sign01,sign02,sign03,most_exposed_tri;

                coords[0]=plo_arr[0]+i*dx[0];
                coords[1]=plo_arr[1]+j*dx[1];
                coords[2]=plo_arr[2]+k*dx[2];

                for(int dim=0;dim<3;dim++)
                {
                    for(int j=0;j<3;j++)
                    {
                        coords_t[dim] += STLtools::eigdirs[3*dim+j]*coords[j];
                    }
                }
                if(    (coords_t[0]>STLtools::bbox_lo[0]) && (coords_t[0]<STLtools::bbox_hi[0])
                        && (coords_t[1]>STLtools::bbox_lo[1]) && (coords_t[1]<STLtools::bbox_hi[1]) 
                        && (coords_t[2]>STLtools::bbox_lo[2]) && (coords_t[2]<STLtools::bbox_hi[2]) )
                {
                    Real mindist=1e50;
                    Real dist,dist1; 
                    Real maxdotpdt=-1e50;
                    Real mat[9],dotpdt,most_exposed_dotpdt,transf_pts[9];
                    int mintri_id=0;
                    Real ptvec1[3],ptvec2[3],ptvec3[3],norm[3];

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

                        for(int n=0;n<9;n++)
                        {
                            mat[n]=STLtools::tri_mat[tr*9+n];
                            transf_pts[n]=STLtools::rot_tri_pts[tr*9+n];
                        }

                        //dist=STLtools::point_tri_distance(coords,t1,t2,t3,
                        //        STLtools::domndir[tr],mat);
                        dist=STLtools::point_tri_distance(coords,transf_pts,
                                STLtools::domndir[tr],mat);
                        //dist=STLtools::point_tri_distance(coords,t1,t2,t3);

                        norm[0]=STLtools::tri_normals[tr*3+0];
                        norm[1]=STLtools::tri_normals[tr*3+1];
                        norm[2]=STLtools::tri_normals[tr*3+2];

                        centpvec[0]=coords[0]-0.3333*(t1[0]+t2[0]+t3[0]);
                        centpvec[1]=coords[1]-0.3333*(t1[1]+t2[1]+t3[1]);
                        centpvec[2]=coords[2]-0.3333*(t1[2]+t2[2]+t3[2]);

                        dotpdt=norm[0]*centpvec[0] + norm[1]*centpvec[1] + norm[2]*centpvec[2];
                        dotpdt=dotpdt/sqrt(centpvec[0]*centpvec[0] + centpvec[1]*centpvec[1] + centpvec[2]*centpvec[2]);

                        if(dist<mindist)
                        {
                            mindist=dist;
                            mintri_id=tr;
                        }
                        if(fabs(dotpdt)>maxdotpdt)
                        {
                            maxdotpdt=fabs(dotpdt);
                            most_exposed_dotpdt=dotpdt;
                            most_exposed_tri=tr;
                        }
                    }

                    sign01=amrex::Math::copysign(1.0,most_exposed_dotpdt);
                    lsphi_arr(i,j,k)=mindist;

                    int num_intersects=0;
                    int sign1=1.0;
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

                            num_intersects += (1-STLtools::lineseg_tri_intersect(outp,coords,t1,t2,t3));
                        }
                        if(num_intersects%2 == 1)
                        {
                            sign1=-1.0;
                        }
                        else
                        {
                            sign1=1.0;
                        }
                        //lsphi_arr(i,j,k)*=sign1;

                        lsphi_arr(i,j,k)=0.0;
                    }
                    else
                    {
                        lsphi_arr(i,j,k)=1.0;
                        //lsphi_arr(i,j,k)=sqrt(len_arr[0]*len_arr[0]+len_arr[1]*len_arr[1]+len_arr[2]*len_arr[2]);
                    }


                    /*if(sign01!=sign1)
                      {
                      Print()<<"sign1,sign01,sign02,sign03:"<<sign1<<"\t"<<sign01<<"\t"<<sign02<<"\t"<<sign03<<"\t"
                      <<mintri_id<<"\t"<<STLtools::DotProd(ptvec1,norm)<<"\t"<<STLtools::DotProd(ptvec2,norm)
                      <<"\t"<<STLtools::DotProd(ptvec3,norm)<<"\n";
                      }*/
                    }); 
        }

        //write plot file
        const std::string& pltfile = "plt";
        WriteSingleLevelPlotfile(pltfile, lsphi, {"phi"}, geom, 0.0, 0);

    }

    amrex::Finalize();
    return(0);
}
