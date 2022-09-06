#ifndef WALLTEMP_H_
#define WALLTEMP_H_

#include <AMReX.H>
#include <Constants.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real get_wall_temperature(amrex::Real x,amrex::Real y,amrex::Real z,int dir,amrex::Real c[3])
{
   amrex::Real temp=NTP_TEMP;
   amrex::Real pos[3]={x,y,z};

   temp=c[0]+c[1]*pos[dir]+c[2]*pos[dir]*pos[dir];
   return(temp);
}

#endif
