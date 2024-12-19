#include <BDEM_ParticleContainer.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rotate_vector_to_inertial (amrex::Real quats[4], amrex::Real v_body[THREEDIM], amrex::Real v_inert[THREEDIM])
{
    // Rotates a vector from a particle's body-fixed reference frame to the inertial frame 

    amrex::Real q0 = quats[0];
    amrex::Real q1 = quats[1];
    amrex::Real q2 = quats[2];
    amrex::Real q3 = quats[3];
     
    v_inert[XDIR] = (q0*q0 + q1*q1 - q2*q2 - q3*q3)*v_body[XDIR] + 2.0*(q1*q2  - q0*q3)*v_body[YDIR] + 2.0*(q1*q3 + q0*q2)*v_body[ZDIR];
    v_inert[YDIR] = 2.0*(q1*q2 + q0*q3)*v_body[XDIR] + (q0*q0 - q1*q1 + q2*q2 - q3*q3)*v_body[YDIR] + 2.0*(q2*q3 - q0*q1)*v_body[ZDIR];
    v_inert[ZDIR] = 2.0*(q1*q3 - q0*q2)*v_body[XDIR] + 2.0*(q2*q3 + q0*q1)*v_body[YDIR] + (q0*q0 - q1*q1 - q2*q2 + q3*q3)*v_body[ZDIR];
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_bonded_particle_pos(int type, int comp, amrex::Real radius, 
                             amrex::Real bp_pos[THREEDIM], 
                             amrex::Real quats[4], 
                             amrex::Real pc_pos[THREEDIM]){
    
  amrex::Real pos[THREEDIM];
  int j = comp;

  // Get the nominal origin-centered particle positions
  if(type == 0){
    pos[XDIR] = 0.0; pos[YDIR] = 0.0; pos[ZDIR] = 0.0;
  } 
  if(type == 1){
    pos[XDIR] = 2.0*(0.5 - j); pos[YDIR] = 0.0; pos[ZDIR] = 0.0;
  }
  if(type == 2){
    pos[XDIR] = 2.0*(1.0 - j); pos[YDIR] = 0.0; pos[ZDIR] = 0.0;
  }
  if(type == 3){
    pos[XDIR] = 2.0*(1.5 - j); pos[YDIR] = 0.0; pos[ZDIR] = 0.0;
  }
  if(type == 4){
    pos[XDIR] = 2.0*(2.0 - j); pos[YDIR] = 0.0; pos[ZDIR] = 0.0;
  }
  if(type == 5){
    pos[XDIR] = 2.0*(2.5 - j); pos[YDIR] = 0.0; pos[ZDIR] = 0.0;
  }
  if(type == 6){
    pos[XDIR] = 2.0*(3.0 - j); pos[YDIR] = 0.0; pos[ZDIR] = 0.0;
  }
  if(type == 7){
    pos[XDIR] = 2.0*(3.5 - j); pos[YDIR] = 0.0; pos[ZDIR] = 0.0;
  }
  if(type == 8){
    pos[XDIR] = 2.0*(4.0 - j); pos[YDIR] = 0.0; pos[ZDIR] = 0.0;
  }
  if(type == 9){
    pos[XDIR] = 2.0*(4.5 - j); pos[YDIR] = 0.0; pos[ZDIR] = 0.0;
  }
  if(type == 10){
    amrex::Real xval = (j<3) ? 2.0*(1.0-j):(j<7) ? 2.0*(4.5-j):2.0*(8.0-j);
    amrex::Real yval = (j<3) ?     sqthree:(j<7) ?         0.0:-sqthree;
    pos[XDIR] = xval; pos[YDIR] = yval; pos[ZDIR] = 0.0;
  }
  if(type == 11){
    amrex::Real xval = (j<8) ? 2.0*(3.5-j):2.0*(11.5-j);
    amrex::Real yval = (j<8) ? 1.0:-1.0;
    pos[XDIR] = xval; pos[YDIR] = yval; pos[ZDIR] = 0.0;
  }
  if(type == 12){
    amrex::Real xval = (j<8) ? 2.0*(3.5-j):(j<16) ? 2.0*(11.5-j):(j<24) ? 2.0*(19.5-j):2.0*(27.5-j);
    amrex::Real yval = (j<8 || (j>15 && j<24)) ? 1.0:-1.0;
    amrex::Real zval = (j<16) ? 1.0:-1.0;
    pos[XDIR] = xval; pos[YDIR] = yval; pos[ZDIR] = zval;
  }
  if(type == 13){
    amrex::Real xval = (j<12) ? 2.0*(5.5-j):(j<24) ? 2.0*(17.5-j):2.0*(29.5-j);
    amrex::Real yval = (j<12) ? 2.0:(j<24) ? 0.0:-2.0;
    pos[XDIR] = xval; pos[YDIR] = yval; pos[ZDIR] = 0.0;
  }
  if(type == 14){
    amrex::Real xval = (j<12) ? 2.0*(5.5-j):(j<24) ? 2.0*(17.5-j):(j<36) ? 2.0*(29.5-j):(j<48) ? 2.0*(41.5-j):(j<60) ? 2.0*(53.5-j):2.0*(65.5-j);
    amrex::Real yval = (j<12 || (j>35 && j<48)) ? 2.0:(j<24 || (j>47 && j<60)) ? 0.0:-2.0;
    amrex::Real zval = (j<36) ? 1.0:-1.0;
    pos[XDIR] = xval; pos[YDIR] = yval; pos[ZDIR] = zval;
  }
  if(type == 15){
    amrex::Real xval = (j<11) ? 2.0*(5.0-j):(j<22) ? 2.0*(16.0-j):(j<34) ? 2.0*(27.5-j):(j<46) ? 2.0*(39.5-j):(j<58) ? 2.0*(51.5-j):(j<69) ? 2.0*(63.0-j):2.0*(74.0-j);
    amrex::Real yval = (j<11 || (j>57 && j<69)) ? 1.0:(j<22 || j>68) ? -1.0:(j<34) ? 2.0:(j<46) ? 0.0:-2.0;
    amrex::Real zval = (j<22) ? sqtwo:(j<58) ? 0.0:-sqtwo;
    pos[XDIR] = xval; pos[YDIR] = yval; pos[ZDIR] = zval;
  }
  if(type == 16){
    amrex::Real xval = (j<4) ? 2.0*(1.5-j):2.0*(5.5-j);
    amrex::Real yval = (j<4) ? 1.0:-1.0;
    pos[XDIR] = xval; pos[YDIR] = yval; pos[ZDIR] = 0.0;
  }
  if(type == 17){
    amrex::Real xval = (j<4) ? 2.0*(1.5-j):(j<8) ? 2.0*(5.5-j):2.0*(9.5-j);
    amrex::Real yval = (j<4) ? 2.0:(j<8) ? 0:-2.0;
    pos[XDIR] = xval; pos[YDIR] = yval; pos[ZDIR] = 0.0;
  }
  if(type == 18){
    amrex::Real xval = (j<6) ? 2.0*(2.5-j):2.0*(8.5-j);
    amrex::Real yval = (j<6) ? 1.0:-1.0;
    pos[XDIR] = xval; pos[YDIR] = yval; pos[ZDIR] = 0.0;
  }
  if(type == 19){
    amrex::Real xval = (j<8) ? 2.0*(3.5-j):(j<16) ? 2.0*(11.5-j):(j<24) ? 2.0*(19.5-j):2.0*(27.5-j);
    amrex::Real yval = (j<8) ? 3.0:(j<16) ? 1.0:(j<24) ? -1.0:-3.0;
    pos[XDIR] = xval; pos[YDIR] = yval; pos[ZDIR] = 0.0;
  }
  if(type == 20){
    amrex::Real xval = (j<8) ? 2.0*(3.5-j):(j<16) ? 2.0*(11.5-j):(j<24) ? 2.0*(19.5-j):(j<32) ? 2.0*(27.5-j):2.0*(35.5-j);
    amrex::Real yval = (j<8) ? 4.0:(j<16) ? 2.0:(j<24) ? 0.0:(j<32) ? -2.0:-4.0;
    pos[XDIR] = xval; pos[YDIR] = yval; pos[ZDIR] = 0.0;
  }
  if(type == 21){
    amrex::Real xval = (j<4) ? 2.0*(1.5-j)+0.25:2.0*(5.5-j)-0.25;
    amrex::Real yval = (j<4) ? 0.9682458:-0.9682458;
    pos[XDIR] = xval; pos[YDIR] = yval; pos[ZDIR] = 0.0;
  }
  if(type == 22){
    amrex::Real xval = (j<5) ? 2.0*(2.0-j)+0.25:2.0*(7.0-j)-0.25;
    amrex::Real yval = (j<5) ? 0.9682458:-0.9682458;
    pos[XDIR] = xval; pos[YDIR] = yval; pos[ZDIR] = 0.0;
  }
  if(type == 23){
    amrex::Real xval = (j<6) ? 2.0*(2.5-j)+0.25:2.0*(8.5-j)-0.25;
    amrex::Real yval = (j<6) ? 0.9682458:-0.9682458;
    pos[XDIR] = xval; pos[YDIR] = yval; pos[ZDIR] = 0.0;
  }
  if(type == 24){
    amrex::Real xval = (j<8) ? 2.0*(3.5-j)+0.25:(j<16) ? 2.0*(11.5-j)-0.25:(j<24) ? 2.0*(19.5-j)+0.25:2.0*(27.5-j)-0.25;
    amrex::Real yval = (j<8) ? 2.90473751:(j<16) ? 0.9682458:(j<24) ? -0.9682458:-2.90473751;
    pos[XDIR] = xval; pos[YDIR] = yval; pos[ZDIR] = 0.0;
  }
  if(type == 25){
    amrex::Real xval = (j<10) ? 2.0*(4.5-j)+0.25:(j<20) ? 2.0*(14.5-j)-0.25:(j<30) ? 2.0*(24.5-j)+0.25:2.0*(34.5-j)-0.25;
    amrex::Real yval = (j<10) ? 2.90473751:(j<20) ? 0.9682458:(j<30) ? -0.9682458:-2.90473751;
    pos[XDIR] = xval; pos[YDIR] = yval; pos[ZDIR] = 0.0;
  }
  if(type == 26){
    amrex::Real xval = (j<6) ? 2.0*j:(j==6) ? 0.5:(j==7) ? 9.5:(j==8) ? 1.0:9.0;
    amrex::Real yval = (j<6) ? 0.0:(j<8) ? 1.9365:3.8730;
    pos[XDIR] = xval; pos[YDIR] = yval; pos[ZDIR] = 0.0;
  }
  if(type == 27){
    amrex::Real xval = (j<8) ? 2.0*j:(j==8) ? 0.5:(j==9) ? 13.5:(j==10) ? 1.0:13.0;
    amrex::Real yval = (j<8) ? 0.0:(j<10) ? 1.9365:3.8730;
    pos[XDIR] = xval; pos[YDIR] = yval; pos[ZDIR] = 0.0;
  }
  if(type == 28){
    amrex::Real xval = (j<10) ? 2.0*j:(j==10) ? 0.5:(j==11) ? 17.5:(j==12) ? 1.0:(j==13) ? 17.0:(j==14) ? 1.5:16.5;
    amrex::Real yval = (j<10) ? 0.0:(j<12) ? 1.9365:(j<14) ? 3.8730:5.8095;
    pos[XDIR] = xval; pos[YDIR] = yval; pos[ZDIR] = 0.0;
  }

  // Calculate rotated particle position
  rotate_vector_to_inertial(quats, pos, pc_pos);

  // Scale and translate particle position using radius and bp center
  pc_pos[XDIR] *= radius;
  pc_pos[YDIR] *= radius;
  pc_pos[ZDIR] *= radius;
  pc_pos[XDIR] += bp_pos[XDIR];
  pc_pos[YDIR] += bp_pos[YDIR];
  pc_pos[ZDIR] += bp_pos[ZDIR];
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void crosspdt (const amrex::Real* v1, const amrex::Real* v2, amrex::Real* v3)
{
    v3[XDIR] = v1[YDIR]*v2[ZDIR] - v1[ZDIR]*v2[YDIR];
    v3[YDIR] = v1[ZDIR]*v2[XDIR] - v1[XDIR]*v2[ZDIR];
    v3[ZDIR] = v1[XDIR]*v2[YDIR] - v1[YDIR]*v2[XDIR];
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real dotpdt (const amrex::Real* v1, const amrex::Real* v2)
{
    return(v1[XDIR]*v2[XDIR] + v1[YDIR]*v2[YDIR] + v1[ZDIR]*v2[ZDIR]);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_relative_velocities (const BDEMParticleContainer::ParticleType& p1, 
               const BDEMParticleContainer::ParticleType& p2,
               amrex::Real& vrn, amrex::Real* vslip,
               const amrex::Real* normal, 
               const amrex::Real dist_mag)
{
    // translational relative velocity
    amrex::Real vreltrans[3];
    vreltrans[0] =  p1.rdata(realData::xvel) - p2.rdata(realData::xvel);
    vreltrans[1] =  p1.rdata(realData::yvel) - p2.rdata(realData::yvel);
    vreltrans[2] =  p1.rdata(realData::zvel) - p2.rdata(realData::zvel);

    // distance from the contact point to the particle centers
    amrex::Real dist_c1 = 0.5*(dist_mag*dist_mag 
                    + p1.rdata(realData::radius)*p1.rdata(realData::radius) 
                    - p2.rdata(realData::radius)*p2.rdata(realData::radius)) / dist_mag;

    amrex::Real dist_c2 = dist_mag - dist_c1;

    amrex::Real omega_sum[3];
    omega_sum[0] = p1.rdata(realData::xangvel) * dist_c1 + p2.rdata(realData::xangvel) * dist_c2;
    omega_sum[1] = p1.rdata(realData::yangvel) * dist_c1 + p2.rdata(realData::yangvel) * dist_c2;
    omega_sum[2] = p1.rdata(realData::zangvel) * dist_c1 + p2.rdata(realData::zangvel) * dist_c2;

    // the rotational relative velocity
    amrex::Real v_rot[3];
    crosspdt(omega_sum, normal, v_rot);
    

    // total relative velocity
    vreltrans[0] += v_rot[0];
    vreltrans[1] += v_rot[1];
    vreltrans[2] += v_rot[2];

    // normal component of relative velocity
    vrn = dotpdt(vreltrans, normal);

    vslip[0] = vreltrans[0] - vrn*normal[0];
    vslip[1] = vreltrans[1] - vrn*normal[1];
    vslip[2] = vreltrans[2] - vrn*normal[2];    
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real basisval(int l,int m,int n, int i, int j, int k,
                     amrex::Real xp[AMREX_SPACEDIM], 
                     const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo, 
                     const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx)
{
    amrex::Real lval,mval,nval;
    amrex::Real rx,ry,rz;

    lval=0.0;
    mval=0.0;
    nval=0.0;

    rx 	= (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    ry  = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    rz  = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];

    lval=(l==0)?(one-rx):rx;
    mval=(m==0)?(one-ry):ry;
    nval=(n==0)?(one-rz):rz;

    return(lval*mval*nval);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real basisvalder(int dir, int l,int m,int n, int i, int j, int k, 
                        amrex::Real xp[AMREX_SPACEDIM], 
                        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo, 
                        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx)
{
    amrex::Real lval,mval,nval;
    amrex::Real rx,ry,rz;
    amrex::Real dxinv[AMREX_SPACEDIM];

    dxinv[XDIR]=one/dx[XDIR];
    dxinv[YDIR]=one/dx[YDIR];
    dxinv[ZDIR]=one/dx[ZDIR];

    rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];

    lval=(l==0)?(one-rx):rx;
    mval=(m==0)?(one-ry):ry;
    nval=(n==0)?(one-rz):rz;

    if(dir==0)
    {
        lval=(l==0)?-dxinv[XDIR]:dxinv[XDIR];
    }

    else if(dir==1)
    {
        mval=(m==0)?-dxinv[YDIR]:dxinv[YDIR];
    }
    else
    {
        nval=(n==0)?-dxinv[ZDIR]:dxinv[ZDIR];
    }

    return(lval*mval*nval);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real get_levelset_value (amrex::Real ptxyz[THREEDIM],
                                const int lsref,
                                amrex::Array4<amrex::Real const> const& phi,
                                amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& problo,
                                amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
{
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> ls_dx;
    
    ls_dx[XDIR]=dx[XDIR]/amrex::Real(lsref);
    ls_dx[YDIR]=dx[YDIR]/amrex::Real(lsref);
    ls_dx[ZDIR]=dx[ZDIR]/amrex::Real(lsref);
    
    int i=amrex::Math::floor((ptxyz[XDIR]-problo[XDIR]+TINYVAL)/ls_dx[XDIR]);
    int j=amrex::Math::floor((ptxyz[YDIR]-problo[YDIR]+TINYVAL)/ls_dx[YDIR]);
    int k=amrex::Math::floor((ptxyz[ZDIR]-problo[ZDIR]+TINYVAL)/ls_dx[ZDIR]);

    amrex::Real lsval=0.0;
    for(int n=0;n<2;n++)
    {
        for(int m=0;m<2;m++)
        {
            for(int l=0;l<2;l++)
            {
                lsval += phi(i+l,j+m,k+n)*
                    basisval(l,m,n,i,j,k,ptxyz,problo,ls_dx);
            }
        }
    }

    return(lsval);

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void get_levelset_normal(
        amrex::Real ptxyz[THREEDIM],
        const int lsref,
        amrex::Array4<amrex::Real const> const& phi,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> problo,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
        amrex::Real lsnormal[AMREX_SPACEDIM])
{
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> ls_dx;

    ls_dx[XDIR]=dx[XDIR]/lsref;
    ls_dx[YDIR]=dx[YDIR]/lsref;
    ls_dx[ZDIR]=dx[ZDIR]/lsref;
    
    int i=amrex::Math::floor((ptxyz[XDIR]-problo[XDIR]+TINYVAL)/ls_dx[XDIR]);
    int j=amrex::Math::floor((ptxyz[YDIR]-problo[YDIR]+TINYVAL)/ls_dx[YDIR]);
    int k=amrex::Math::floor((ptxyz[ZDIR]-problo[ZDIR]+TINYVAL)/ls_dx[ZDIR]);

    lsnormal[XDIR]=0.0;
    lsnormal[YDIR]=0.0;
    lsnormal[ZDIR]=0.0;
    for(int n=0;n<2;n++)
    {
        for(int m=0;m<2;m++)
        {
            for(int l=0;l<2;l++)
            {
                lsnormal[XDIR] += phi(i+l,j+m,k+n)*
                    basisvalder(XDIR,l,m,n,i,j,k,ptxyz,problo,ls_dx);
                lsnormal[YDIR] += phi(i+l,j+m,k+n)*
                    basisvalder(YDIR,l,m,n,i,j,k,ptxyz,problo,ls_dx);
                lsnormal[ZDIR] += phi(i+l,j+m,k+n)*
                    basisvalder(ZDIR,l,m,n,i,j,k,ptxyz,problo,ls_dx);
            }
        }
    }
    
    amrex::Real ls_inv_norm = 1.0 / sqrt(dotpdt(lsnormal, lsnormal));
    lsnormal[0] *= ls_inv_norm;
    lsnormal[1] *= ls_inv_norm;
    lsnormal[2] *= ls_inv_norm;
}

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_wall_forces(amrex::Real walldist, amrex::Real normal[THREEDIM],
        BDEMParticleContainer::ParticleType& p,amrex::Real dt,
        amrex::Real fn[THREEDIM],amrex::Real ft[THREEDIM], 
        amrex::Real rollfric_torq_by_d[THREEDIM], amrex::Real pos[THREEDIM],
        int contact_law, int liquid_bridging)
{
    Real rp = p.rdata(realData::radius);
    Real overlap_n = rp - walldist;

    // Location of the point of contact with the wall 
    Real contact_loc[THREEDIM];
    contact_loc[XDIR] = pos[XDIR] + walldist*normal[XDIR];
    contact_loc[YDIR] = pos[YDIR] + walldist*normal[YDIR];
    contact_loc[ZDIR] = pos[ZDIR] + walldist*normal[ZDIR];

    // Vector point from the particle center of mass to the point of contact
    Real contact_vect[THREEDIM];
    contact_vect[XDIR] = contact_loc[XDIR] - p.pos(0);
    contact_vect[YDIR] = contact_loc[YDIR] - p.pos(1);
    contact_vect[ZDIR] = contact_loc[ZDIR] - p.pos(2);

    // Calculate the velocity at the point of contact
    Real vrel[THREEDIM];
    Real angvel[THREEDIM] = {p.rdata(realData::xangvel), p.rdata(realData::yangvel), p.rdata(realData::zangvel)};
    Real cpdt[THREEDIM];
    crosspdt(angvel, contact_vect, cpdt);
    vrel[XDIR] = p.rdata(realData::xvel) + cpdt[XDIR];
    vrel[YDIR] = p.rdata(realData::yvel) + cpdt[YDIR];
    vrel[ZDIR] = p.rdata(realData::zvel) + cpdt[ZDIR];

    Real vreltrans_norm = dotpdt(vrel, normal);

    Real vrel_t[3];
    vrel_t[0] = vrel[0] - vreltrans_norm*normal[0];
    vrel_t[1] = vrel[1] - vreltrans_norm*normal[1];
    vrel_t[2] = vrel[2] - vreltrans_norm*normal[2];

    // Contact law: default - linear, 1 - hertz
    Real meff = p.rdata(realData::mass);
    Real Reff = rp;
    Real kn_dem_w;
    Real etan_dem_w;
    Real kt_dem_w;
    Real etat_dem_w;
    if (contact_law == CONTACT_HERTZ)
    {
        Real E1 = p.rdata(realData::E);
        Real E2 = DEM::E_wall;
        E2 = (E2==100.0) ? E1 : E2; // E_wall=E_p if not specified
        Real nu1 = p.rdata(realData::nu);
        Real nu2 = DEM::nu_wall;
        Real E_eff = E1*E2/((1.0-nu1*nu1)*E2+(1.0-nu2*nu2)*E1);
        Real temp = 2.0*E2*(2.0-nu1)*(1.0+nu1) + 2.0*E1*(2.0-nu2)*(1.0+nu2);
        Real G_eff = E1*E2/temp;
        Real beta = log(DEM::e_n) / sqrt(PI*PI + log(DEM::e_n)*log(DEM::e_n));
        Real temp1 = sqrt(Reff*overlap_n);
        kn_dem_w = 4.0/3.0*E_eff*temp1;
        kt_dem_w = 8.0*G_eff*temp1;
        etan_dem_w = -2.0*sqrt(5.0/3.0)*beta*sqrt(meff*E_eff)*sqrt(temp1);
        etat_dem_w = -4.0*sqrt(5.0/3.0)*beta*sqrt(meff*G_eff)*sqrt(temp1);
    }
    else
    {            
        kn_dem_w = DEM::k_n_wall;
        etan_dem_w = two*sqrt(DEM::k_n_wall*meff)*fabs(log(DEM::e_n_wall));
        etan_dem_w /= sqrt(PI*PI + log(DEM::e_n_wall)*log(DEM::e_n_wall));

        kt_dem_w = DEM::k_t_wall;
        etat_dem_w = two*sqrt(DEM::k_t_wall*meff)*fabs(log(DEM::e_t_wall));
        etat_dem_w /= sqrt(PI*PI + log(DEM::e_t_wall)*log(DEM::e_t_wall));
    }

    Real overlap_t[3];
    Real mag_overlap_t;

    // calculate the normal contact force
    fn[0] = -(kn_dem_w*overlap_n*normal[0]
            + etan_dem_w*vreltrans_norm*normal[0]);
    fn[1] = -(kn_dem_w*overlap_n*normal[1]
            + etan_dem_w*vreltrans_norm*normal[1]);
    fn[2] = -(kn_dem_w*overlap_n*normal[2]
            + etan_dem_w*vreltrans_norm*normal[2]);

    // calculate the tangential displacement
    overlap_t[0] = dt*vrel_t[0];
    overlap_t[1] = dt*vrel_t[1];
    overlap_t[2] = dt*vrel_t[2];

    mag_overlap_t = sqrt(dotpdt(overlap_t, overlap_t));

    if (mag_overlap_t > zero)
    {
        Real mu_wall = DEM::mu_wall;

        // Correct for coarse grained particles
        if ( DEM::particles_in_parcel > 1 )
        {
            Real omega_p  [THREEDIM];
            omega_p[XDIR] = p.rdata(realData::xangvel);
            omega_p[YDIR] = p.rdata(realData::yangvel);
            omega_p[ZDIR] = p.rdata(realData::zangvel); 


            Real omega_p_cp [THREEDIM];
            CrossProd(omega_p, normal, omega_p_cp);

            Real oneOvera = 1./DEM::CG_ratio;
            Real beta = 1./std::sqrt(oneOvera*oneOvera + ( 1. -  oneOvera)*( 1. -  oneOvera));

            Real mag_overlap_orig = overlap_n*oneOvera;
            Real dorig = rp*oneOvera*2.;
            Real denomVec[THREEDIM];

            for (int i = 0; i < 3; i++)
            {
                denomVec[i] = vrel_t[i] +  beta*mag_overlap_orig*omega_p_cp[i];
            }

            Real num = (DEM::CG_ratio - 1.)*beta*dorig*normVec(omega_p_cp);
            Real denom = normVec(denomVec);

            if (denom > TINYVAL) 
            {
                mu_wall += DEM::mu * num/denom;  
            }
        }

        Real fnmd = mu_wall * sqrt(dotpdt(fn, fn));
        Real ftsp[THREEDIM];
        Real ftn[THREEDIM];
        Real mag_ftsp;
        Real mag_ftn;
        Real tangent[THREEDIM];
        tangent[XDIR] = overlap_t[XDIR]/mag_overlap_t;
        tangent[YDIR] = overlap_t[YDIR]/mag_overlap_t;
        tangent[ZDIR] = overlap_t[ZDIR]/mag_overlap_t;
        // calculate the tangential contact force
        Real overlap_t_tan = dotpdt(overlap_t,tangent);
        Real vrel_t_tan = dotpdt(vrel_t,tangent);
        // calculate the tangential contact force
        ftsp[XDIR] = -(kt_dem_w*overlap_t_tan*tangent[XDIR]
                + etat_dem_w*vrel_t_tan*tangent[XDIR]);
        ftsp[YDIR] = -(kt_dem_w*overlap_t_tan*tangent[YDIR]
                + etat_dem_w*vrel_t_tan*tangent[YDIR]);
        ftsp[ZDIR] = -(kt_dem_w*overlap_t_tan*tangent[ZDIR]
                + etat_dem_w*vrel_t_tan*tangent[ZDIR]);

        ftn[XDIR] = -fnmd * tangent[XDIR];
        ftn[YDIR] = -fnmd * tangent[YDIR];
        ftn[ZDIR] = -fnmd * tangent[ZDIR];
        mag_ftsp = sqrt(dotpdt(ftsp,ftsp));
        mag_ftn = sqrt(dotpdt(ftn,ftn));
        if (mag_ftsp < mag_ftn)
        {
            ft[XDIR] = ftsp[XDIR];
            ft[YDIR] = ftsp[YDIR];
            ft[ZDIR] = ftsp[ZDIR];
        }
        else
        {
            ft[XDIR] = ftn[XDIR];
            ft[YDIR] = ftn[YDIR];
            ft[ZDIR] = ftn[ZDIR];
        }
    }  
    else 
    {
        ft[XDIR] = zero;
        ft[YDIR] = zero;
        ft[ZDIR] = zero;
    }

    rollfric_torq_by_d[XDIR]=zero;
    rollfric_torq_by_d[YDIR]=zero;
    rollfric_torq_by_d[ZDIR]=zero;
    Real omega_p1[THREEDIM];
    omega_p1[XDIR] = p.rdata(realData::xangvel);
    omega_p1[YDIR] = p.rdata(realData::yangvel);
    omega_p1[ZDIR] = p.rdata(realData::zangvel);
    Real mag_omegap1 = sqrt(dotpdt(omega_p1,omega_p1));
    if (mag_omegap1 > zero)
    {
        Real mag_fn = sqrt(dotpdt(fn,fn));
        Real fric_force   = -DEM::muR * mag_fn;
        rollfric_torq_by_d[XDIR]=fric_force*omega_p1[XDIR]/mag_omegap1;
        rollfric_torq_by_d[YDIR]=fric_force*omega_p1[YDIR]/mag_omegap1;
        rollfric_torq_by_d[ZDIR]=fric_force*omega_p1[ZDIR]/mag_omegap1;
    }

    p.rdata(realData::fx) += fn[XDIR] + ft[XDIR];
    p.rdata(realData::fy) += fn[YDIR] + ft[YDIR];
    p.rdata(realData::fz) += fn[ZDIR] + ft[ZDIR];

    Real ftot[THREEDIM] = {fn[XDIR]+ft[XDIR], fn[YDIR]+ft[YDIR], fn[ZDIR]+ft[ZDIR]};
    Real torq[THREEDIM];
    crosspdt(contact_vect, ftot, torq);
    p.rdata(realData::taux) += torq[XDIR] + walldist*rollfric_torq_by_d[XDIR];
    p.rdata(realData::tauy) += torq[YDIR] + walldist*rollfric_torq_by_d[YDIR];
    p.rdata(realData::tauz) += torq[ZDIR] + walldist*rollfric_torq_by_d[ZDIR];

    /* Calculating liquid bridge force for particles in contact with wall
          - Assumes bridges are only present when particle is in direct contact with wall
          - Assumes there is enough free liquid for bridge to form
          - FIXME: Currently, bridge is calculated assuming the wall surface acts as a sphere with the same radius
    */
    if(liquid_bridging){ 
      // Bridge volume 
      Real Vtot = (p.rdata(realData::liquid_volume))*pow(0.5,0.5);

      // Calculate the fitting coefficients from Soulie et al. (2006)
      Real acoef = -1.1*pow( Vtot/pow(rp,3.0) ,-0.53);
      Real bcoef = (-0.148*log(Vtot/pow(rp,3.0)) - 0.96)*pow(DEM::contact_angle,two) - 0.0082*log(Vtot/pow(rp,3.0)) + 0.48;
      Real ccoef = 0.0018*log(Vtot/pow(rp,3.0)) + 0.078;

      // Calculate the normal capillary force
      Real Fcap = PI*DEM::wall_gamma*rp*(ccoef + exp(bcoef));

      p.rdata(realData::fx) += Fcap*normal[XDIR];
      p.rdata(realData::fy) += Fcap*normal[YDIR];
      p.rdata(realData::fz) += Fcap*normal[ZDIR];
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void solve_fibrils_equation
(
    Real vp[3],
    Real vn[3],
    Real rad1,
    Real rad2,
    Real &fractionOfFibrils
)
{
    

}
