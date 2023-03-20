#include <BDEM_ParticleContainer.H>

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
amrex::Real get_levelset_value (BDEMParticleContainer::ParticleType const& p,
                                const int lsref,
                                amrex::Array4<amrex::Real const> const& phi,
                                amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& problo,
                                amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
{
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> ls_dx;
    
    amrex::Real ptxyz[3] = {p.pos(0),p.pos(1),p.pos(2)};

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
        BDEMParticleContainer::ParticleType const& p,
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
    
    amrex::Real ptxyz[3] = {p.pos(0),p.pos(1),p.pos(2)};

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
        amrex::Real rollfric_torq_by_d[THREEDIM],
        int contact_law)
{
    Real v_rot[3];
    v_rot[0] = walldist * p.rdata(realData::xangvel);
    v_rot[1] = walldist * p.rdata(realData::yangvel);
    v_rot[2] = walldist * p.rdata(realData::zangvel);

    Real vreltrans[3];
    Real cprod[3];
    Real rp = p.rdata(realData::radius);
    Real overlap_n = rp - walldist;

    crosspdt(v_rot, normal, cprod);
    vreltrans[0] = p.rdata(realData::xvel) + cprod[0];
    vreltrans[1] = p.rdata(realData::yvel) + cprod[1];
    vreltrans[2] = p.rdata(realData::zvel) + cprod[2];

    Real vreltrans_norm = dotpdt(vreltrans, normal);

    Real vrel_t[3];
    vrel_t[0] = vreltrans[0] - vreltrans_norm*normal[0];
    vrel_t[1] = vreltrans[1] - vreltrans_norm*normal[1];
    vrel_t[2] = vreltrans[2] - vreltrans_norm*normal[2];

    int phase = p.idata(intData::phase);

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
        etan_dem_w = 2.0*sqrt(5.0/3.0)*beta*sqrt(meff*E_eff)*sqrt(temp1);
        etat_dem_w = 4.0*sqrt(5.0/3.0)*beta*sqrt(meff*G_eff)*sqrt(temp1);
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

    Real omega_p1[THREEDIM];
    omega_p1[XDIR] = p.rdata(realData::xangvel);
    omega_p1[YDIR] = p.rdata(realData::yangvel);
    omega_p1[ZDIR] = p.rdata(realData::zangvel);
    Real mag_omegap1 = sqrt(dotpdt(omega_p1,omega_p1));
    rollfric_torq_by_d[XDIR]=zero;
    rollfric_torq_by_d[YDIR]=zero;
    rollfric_torq_by_d[ZDIR]=zero;
    if (mag_omegap1 > zero)
    {
        Real mag_fn = sqrt(dotpdt(fn,fn));
        Real fric_force   = -DEM::muR * mag_fn;
        rollfric_torq_by_d[XDIR]=fric_force*omega_p1[XDIR]/mag_omegap1;
        rollfric_torq_by_d[YDIR]=fric_force*omega_p1[YDIR]/mag_omegap1;
        rollfric_torq_by_d[ZDIR]=fric_force*omega_p1[ZDIR]/mag_omegap1;
    }
}
