#ifndef BDEM_CHECKPAIR_H_
#define BDEM_CHECKPAIR_H_

#include <BDEM_Specs.H>

struct CheckPair
{
    // TODO: May need to modify for glued sphere particles
    template <class P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    bool operator()(const P& p1, const P& p2) const
    {
        bool not_pair = false;
        for(int pc1 = 0; pc1 < p1.idata(intData::num_comp_sphere); pc1++){
            for(int pc2 = 0; pc2 < p2.idata(intData::num_comp_sphere); pc2++){
                Real p1pos_inert[THREEDIM];
                Real p1pos_body[THREEDIM] = {p1.rdata(realData::radius)*((p1.idata(intData::num_comp_sphere) - 1) - 2.0*pc1), 0.0, 0.0};
  
                amrex::Real q0 = p1.rdata(realData::q0);
                amrex::Real q1 = p1.rdata(realData::q1);
                amrex::Real q2 = p1.rdata(realData::q2);
                amrex::Real q3 = p1.rdata(realData::q3);
                 
                p1pos_inert[XDIR] = (q0*q0 - q1*q1 + q2*q2 - q3*q3)*p1pos_body[XDIR] - 2.0*(q1*q2 + q0*q3)*p1pos_body[YDIR] + 2.0*(q2*q3 - q0*q1)*p1pos_body[ZDIR];
                p1pos_inert[YDIR] = -2.0*(q1*q2 - q0*q3)*p1pos_body[XDIR] + (q0*q0 + q1*q1 - q2*q2 - q3*q3)*p1pos_body[YDIR] - 2.0*(q1*q3 + q0*q2)*p1pos_body[ZDIR];
                p1pos_inert[ZDIR] = 2.0*(q2*q3 + q0*q1)*p1pos_body[XDIR] - 2.0*(q1*q3 - q0*q2)*p1pos_body[YDIR] + (q0*q0 - q1*q1 - q2*q2 + q3*q3)*p1pos_body[ZDIR];

                // Calculate the inertial frame position of p1
                Real p1posx = p1.pos(0) + p1pos_inert[XDIR];
                Real p1posy = p1.pos(1) + p1pos_inert[YDIR];
                Real p1posz = p1.pos(2) + p1pos_inert[ZDIR];                

                Real p2pos_inert[THREEDIM];
                Real p2pos_body[THREEDIM] = {p1.rdata(realData::radius)*((p2.idata(intData::num_comp_sphere) - 1) - 2.0*pc2), 0.0, 0.0};

                q0 = p2.rdata(realData::q0);
                q1 = p2.rdata(realData::q1);
                q2 = p2.rdata(realData::q2);
                q3 = p2.rdata(realData::q3);
                 
                p2pos_inert[XDIR] = (q0*q0 - q1*q1 + q2*q2 - q3*q3)*p2pos_body[XDIR] - 2.0*(q1*q2 + q0*q3)*p2pos_body[YDIR] + 2.0*(q2*q3 - q0*q1)*p2pos_body[ZDIR];
                p2pos_inert[YDIR] = -2.0*(q1*q2 - q0*q3)*p2pos_body[XDIR] + (q0*q0 + q1*q1 - q2*q2 - q3*q3)*p2pos_body[YDIR] - 2.0*(q1*q3 + q0*q2)*p2pos_body[ZDIR];
                p2pos_inert[ZDIR] = 2.0*(q2*q3 + q0*q1)*p2pos_body[XDIR] - 2.0*(q1*q3 - q0*q2)*p2pos_body[YDIR] + (q0*q0 - q1*q1 - q2*q2 + q3*q3)*p2pos_body[ZDIR];

                // Calculate the inertial frame position of p2
                Real p2posx = p2.pos(0) + p2pos_inert[XDIR];
                Real p2posy = p2.pos(1) + p2pos_inert[YDIR];
                Real p2posz = p2.pos(2) + p2pos_inert[ZDIR];

                amrex::Real d0 = (p1posx - p2posx);
                amrex::Real d1 = (p1posy - p2posy);
                amrex::Real d2 = (p1posz - p2posz);    
                amrex::Real dsquared = d0*d0 + d1*d1 + d2*d2;   
                if(dsquared <= 9.0*p1.rdata(realData::radius)*p1.rdata(realData::radius)){
                    return true;
                }
            }
        }
        // If we get through all component particles with no pair match, return false
        return not_pair;
    }
};

#endif
