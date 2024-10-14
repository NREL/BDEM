#include<MoveUtils.H>

void move_this_point(amrex::Real coord[3],amrex::Real newcoord[3],amrex::Real timestep,
                     int movetype,amrex::Real u[3],amrex::Real center[3],amrex::Real movevel)
{
  newcoord[0]=coord[0];
  newcoord[1]=coord[1];
  newcoord[2]=coord[2];
   if(movetype==1 || movetype==3)
   {
        newcoord[0] += movevel*timestep*u[0];
        newcoord[1] += movevel*timestep*u[1];
        newcoord[2] += movevel*timestep*u[2];
  
   }
   else if(movetype==2)
   {

        // Rodrigues' rotation formula
        
        amrex::Real theta = movevel*timestep;
        amrex::Real cosTheta = cos(theta);
        amrex::Real sinTheta = sin(theta);

        amrex::Real px = coord[0] - center[0];
        amrex::Real py = coord[1] - center[1]; 
        amrex::Real pz = coord[2] - center[2]; 
        
        amrex::Real rotatedX = (cosTheta + (1. - cosTheta) * u[0] * u[0]) * px +
                      ((1. - cosTheta) * u[0] * u[1] - u[2] * sinTheta) * py +
                      ((1. - cosTheta) * u[0] * u[2] + u[1] * sinTheta) * pz;

        amrex::Real rotatedY = ((1. - cosTheta) * u[1] * u[0] + u[2] * sinTheta) * px +
                        (cosTheta + (1. - cosTheta) * u[1] * u[1]) * py +
                        ((1. - cosTheta) * u[1] * u[2] - u[0] * sinTheta) * pz;

        amrex::Real rotatedZ = ((1. - cosTheta) * u[2] * u[0] - u[1] * sinTheta) * px +
                        ((1. - cosTheta) * u[2] * u[1] + u[0] * sinTheta) * py +
                        (cosTheta + (1. - cosTheta) * u[2] * u[2]) * pz;
        
        newcoord[0] = rotatedX + center[0];
        newcoord[1] = rotatedY + center[1];
        newcoord[2] = rotatedZ + center[2];

        // int t1=(movedir+1)%3;
        // int t2=(movedir+2)%3;
        // amrex::Real rad=std::pow((newcoord[t1]-center[t1]),2.0);
        // rad+=std::pow((newcoord[t2]-center[t2]),2.0);
        // rad=sqrt(rad);
   
        // amrex::Real sin_theta=(newcoord[t2]-center[t2])/rad;
        // amrex::Real cos_theta=(coord[t1]-center[t1])/rad;
        // newcoord[t1] += -movevel*rad*sin_theta*timestep;
        // newcoord[t2] +=  movevel*rad*cos_theta*timestep;
        // newcoord[movedir] += 0.0;
   }
}
