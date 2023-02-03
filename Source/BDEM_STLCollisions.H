amrex::ParallelFor(np,[=]
       AMREX_GPU_DEVICE (int i) noexcept
{
    Real pos[3],pos_t[3],dist;
    Real mat[9],transf_pts[9];
    int min_tri_id=0;
    bool success,collided;
    Real centpvec[3],norm[3];

    ParticleType& p = pstruct[i];
    for(int pc=0; pc<p.idata(intData::num_comp_sphere); pc++){
        Real pos_body[THREEDIM];
        get_inertial_body_fixed_pos(p, pc, pos, pos_body);

        for(int dim=0;dim<3;dim++)
        {
            pos_t[dim]=0.0;
            for(int j=0;j<3;j++)
            {
                pos_t[dim] += STLtools::eigdirs[3*dim+j]*pos[j];
            }
        }

        //amrex::Print()<<"bbox_lo:"<<STLtools::bbox_lo[0]<<"\t"<<STLtools::bbox_lo[1]<<"\t"<<STLtools::bbox_lo[2]<<"\n";
        //amrex::Print()<<"bbox_hi:"<<STLtools::bbox_hi[0]<<"\t"<<STLtools::bbox_hi[1]<<"\t"<<STLtools::bbox_hi[2]<<"\n";

        Real mindist=BIGVAL;
        if(     (pos_t[0]>STLtools::bbox_lo[0]) && 
                (pos_t[0]<STLtools::bbox_hi[0]) &&
                (pos_t[1]>STLtools::bbox_lo[1]) &&
                (pos_t[1]<STLtools::bbox_hi[1]) &&
                (pos_t[2]>STLtools::bbox_lo[2]) &&
                (pos_t[2]<STLtools::bbox_hi[2]) )
        {
            
            int ndistcalcs;

#ifndef AMREX_USE_GPU
            //cpu version, need to figure out gpu version
            STLtools::searchtriangulation(0,STLtools::num_tri-1,
                            STLtools::sorted_indexarray,
                            pos,mindist,min_tri_id,ndistcalcs);

            mindist=std::sqrt(mindist);
            int min_tr=STLtools::sorted_indexarray[min_tri_id];

#else
            STLtools::brutesearch(0,STLtools::num_tri-1,
                                  STLtools::sorted_indexarray,
                                  pos,mindist,min_tri_id);
            int min_tr=STLtools::sorted_indexarray[min_tri_id];
#endif

            //Print()<<"tr:"<<tr<<"\n";
            centpvec[0]=pos[0]-
                0.3333*(STLtools::tri_pts[min_tr*9+0]+STLtools::tri_pts[min_tr*9+3]+STLtools::tri_pts[min_tr*9+6]);
            centpvec[1]=pos[1]-
                0.3333*(STLtools::tri_pts[min_tr*9+1]+STLtools::tri_pts[min_tr*9+4]+STLtools::tri_pts[min_tr*9+7]);
            centpvec[2]=pos[2]-
                0.3333*(STLtools::tri_pts[min_tr*9+2]+STLtools::tri_pts[min_tr*9+5]+STLtools::tri_pts[min_tr*9+8]);

            norm[0]=STLtools::tri_normals[min_tr*3+0];
            norm[1]=STLtools::tri_normals[min_tr*3+1];
            norm[2]=STLtools::tri_normals[min_tr*3+2];

            Real dotpdt=norm[0]*centpvec[0] + norm[1]*centpvec[1] + norm[2]*centpvec[2];
            mindist*=amrex::Math::copysign(1.0,dotpdt);

            if(mindist<p.rdata(realData::radius) && amrex::Math::abs(mindist)<p.rdata(realData::radius))
            //if(mindist<p.rdata(realData::radius))
            {
                /*if(amrex::Math::abs(mindist)>p.rdata(realData::radius))
                {
                    amrex::Print()<<"mindist,rad:"<<mindist<<"\t"<<p.rdata(realData::radius)
                    <<"\t"<<p.pos(0)<<"\t"<<p.pos(1)<<"\t"<<p.pos(2)<<
                    "\t"<<norm[0]<<"\t"<<norm[1]<<"\t"<<norm[2]<<"\t"<<dotpdt<<"\n";
                }*/
                
                Real fn[3],ft[3],rollfriqtorq_by_d[3],tow[3],contact_vect[3];

                norm[0]*=-1;
                norm[1]*=-1;
                norm[2]*=-1;

                get_wall_forces(mindist,norm,p,dt,fn,ft,rollfriqtorq_by_d,pos,contact_vect);

                // each particle updates its force (no need for atomics)
                p.rdata(realData::fx) += fn[XDIR] + ft[XDIR];
                p.rdata(realData::fy) += fn[YDIR] + ft[YDIR];
                p.rdata(realData::fz) += fn[ZDIR] + ft[ZDIR];

                Real ftot[THREEDIM] = {fn[XDIR]+ft[XDIR], fn[YDIR]+ft[YDIR], fn[ZDIR]+ft[ZDIR]};
                Real torq[THREEDIM];
                crosspdt(contact_vect, ftot, torq);
                p.rdata(realData::taux) += torq[XDIR] + mindist*rollfriqtorq_by_d[XDIR];
                p.rdata(realData::tauy) += torq[YDIR] + mindist*rollfriqtorq_by_d[YDIR];
                p.rdata(realData::tauz) += torq[ZDIR] + mindist*rollfriqtorq_by_d[ZDIR];
            }
        }
    }

});
