#include <BDEM_ParticleContainer.H>


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rotate_vector_to_inertial (const BDEMParticleContainer::ParticleType& p, 
               amrex::Real v_body[THREEDIM], amrex::Real v_inert[THREEDIM])
{
    // Rotates a vector from a particle's body-fixed reference frame to the inertial frame 

    amrex::Real q0 = p.rdata(realData::q0);
    amrex::Real q1 = p.rdata(realData::q1);
    amrex::Real q2 = p.rdata(realData::q2);
    amrex::Real q3 = p.rdata(realData::q3);
     
    v_inert[XDIR] = (q0*q0 + q1*q1 - q2*q2 - q3*q3)*v_body[XDIR] + 2.0*(q1*q2  - q0*q3)*v_body[YDIR] + 2.0*(q1*q3 + q0*q2)*v_body[ZDIR];
    v_inert[YDIR] = 2.0*(q1*q2 + q0*q3)*v_body[XDIR] + (q0*q0 - q1*q1 + q2*q2 - q3*q3)*v_body[YDIR] + 2.0*(q2*q3 - q0*q1)*v_body[ZDIR];
    v_inert[ZDIR] = 2.0*(q1*q3 - q0*q2)*v_body[XDIR] + 2.0*(q2*q3 + q0*q1)*v_body[YDIR] + (q0*q0 - q1*q1 - q2*q2 + q3*q3)*v_body[ZDIR];
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rotate_vector_to_body (const BDEMParticleContainer::ParticleType& p, 
               amrex::Real v_inert[THREEDIM], amrex::Real v_body[THREEDIM])
{
    // Rotates a vector from the inertial reference frame to the particle's body-fixed reference frame
 
    amrex::Real q0 = p.rdata(realData::q0);
    amrex::Real q1 = p.rdata(realData::q1);
    amrex::Real q2 = p.rdata(realData::q2);
    amrex::Real q3 = p.rdata(realData::q3);
     
    v_body[XDIR] = (q0*q0 + q1*q1 - q2*q2 - q3*q3)*v_inert[XDIR] + 2.0*(q1*q2 + q0*q3)*v_inert[YDIR] + 2.0*(q1*q3 - q0*q2)*v_inert[ZDIR];
    v_body[YDIR] = 2.0*(q1*q2  - q0*q3)*v_inert[XDIR] + (q0*q0 - q1*q1 + q2*q2 - q3*q3)*v_inert[YDIR] + 2.0*(q2*q3 + q0*q1)*v_inert[ZDIR];
    v_body[ZDIR] = 2.0*(q1*q3 + q0*q2)*v_inert[XDIR] + 2.0*(q2*q3 - q0*q1)*v_inert[YDIR] + (q0*q0 - q1*q1 - q2*q2 + q3*q3)*v_inert[ZDIR];
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_inertial_body_fixed_pos(const BDEMParticleContainer::ParticleType& p, 
               const int pc, amrex::Real p_inert[THREEDIM], amrex::Real p_body[THREEDIM])
{
    // Return the position components of component sphere pc of particle p
    // for both body-fixed and inertial reference frame

    p_body[XDIR] = p.rdata(realData::radius)*((p.idata(intData::num_comp_sphere) - 1) - 2.0*pc);
    p_body[YDIR] = 0.0;
    p_body[ZDIR] = 0.0;
    rotate_vector_to_inertial(p, p_body, p_inert);
    p_inert[XDIR] += p.pos(0);
    p_inert[YDIR] += p.pos(1);
    p_inert[ZDIR] += p.pos(2);
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
        amrex::Real rollfric_torq_by_d[THREEDIM], amrex::Real pos[THREEDIM])
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

    int phase = p.idata(intData::phase);

    Real kn_dem_w = DEM::k_n_wall;
    Real meff = p.rdata(realData::mass);
    Real etan_dem_w = two*sqrt(DEM::k_n_wall*meff)*fabs(log(DEM::e_n_wall));
    etan_dem_w /= sqrt(PI*PI + log(DEM::e_n_wall)*log(DEM::e_n_wall));

    Real kt_dem_w = DEM::k_t_wall;
    Real Reff = rp;
    Real rad1 = rp;
    Real etat_dem_w = two*sqrt(DEM::k_t_wall*meff)*fabs(log(DEM::e_t_wall));
    etat_dem_w /= sqrt(PI*PI + log(DEM::e_t_wall)*log(DEM::e_t_wall));

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
        Real fnmd = DEM::mu_wall * sqrt(dotpdt(fn, fn));
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

    // Only calculate nonzero rolling friction when particle has a single component sphere
    rollfric_torq_by_d[XDIR]=zero;
    rollfric_torq_by_d[YDIR]=zero;
    rollfric_torq_by_d[ZDIR]=zero;
    if(p.idata(intData::num_comp_sphere) == 1){
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
}
