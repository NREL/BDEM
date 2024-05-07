#include<MoveUtils.H>

void move_this_point(amrex::Real coord[3],amrex::Real newcoord[3],amrex::Real timestep,
                     int movetype,int movedir,amrex::Real center[3],amrex::Real movevel)
{
  newcoord[0]=coord[0];
  newcoord[1]=coord[1];
  newcoord[2]=coord[2];
   if(movetype==1 || movetype==3)
   {
        newcoord[movedir] += movevel*timestep;
   }
   else
   if(movetype==2)
   {
        int t1=(movedir+1)%3;
        int t2=(movedir+2)%3;
        amrex::Real rad=std::pow((newcoord[t1]-center[t1]),2.0);
        rad+=std::pow((newcoord[t2]-center[t2]),2.0);
        rad=sqrt(rad);
   
        amrex::Real sin_theta=(newcoord[t2]-center[t2])/rad;
        amrex::Real cos_theta=(coord[t1]-center[t1])/rad;
        newcoord[t1] += -movevel*rad*sin_theta*timestep;
        newcoord[t2] +=  movevel*rad*cos_theta*timestep;
        newcoord[movedir] += 0.0;
   }
}
