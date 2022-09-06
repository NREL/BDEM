#include <BDEM_ParticleContainer.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void cross_product (const amrex::Real* a, const amrex::Real* b, amrex::Real* c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real dot_product (const amrex::Real* a, const amrex::Real* b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void cfrelvel (const BDEMParticleContainer::ParticleType& p1, 
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
    amrex::Real dist_c1 = (dist_mag*dist_mag 
                    + p1.rdata(realData::radius)*p1.rdata(realData::radius) 
                    - p2.rdata(realData::radius)*p2.rdata(realData::radius)) / (2.0*dist_mag);

    amrex::Real dist_c2 = dist_mag - dist_c1;

    amrex::Real omega_sum[3];
    omega_sum[0] = p1.rdata(realData::xangvel) * dist_c1 + p2.rdata(realData::xangvel) * dist_c2;
    omega_sum[1] = p1.rdata(realData::yangvel) * dist_c1 + p2.rdata(realData::yangvel) * dist_c2;
    omega_sum[2] = p1.rdata(realData::zangvel) * dist_c1 + p2.rdata(realData::zangvel) * dist_c2;

    // the rotational relative velocity
    amrex::Real v_rot[3];
    cross_product(omega_sum, normal, v_rot);
    

    // total relative velocity
    vreltrans[0] += v_rot[0];
    vreltrans[1] += v_rot[1];
    vreltrans[2] += v_rot[2];

    // normal component of relative velocity
    vrn = dot_product(vreltrans, normal);

    vslip[0] = vreltrans[0] - vrn*normal[0];
    vslip[1] = vreltrans[1] - vrn*normal[1];
    vslip[2] = vreltrans[2] - vrn*normal[2];    
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real interp_level_set (BDEMParticleContainer::ParticleType const& p,
                              const int n_refine,
                              amrex::Array4<amrex::Real const> const& phi,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
{
    amrex::Real scaled_dxi[3];
    scaled_dxi[0] = n_refine * dxi[0];
    scaled_dxi[1] = n_refine * dxi[1];
    scaled_dxi[2] = n_refine * dxi[2];
    
    amrex::Real x = (p.pos(0) - plo[0]) * scaled_dxi[0];
    amrex::Real y = (p.pos(1) - plo[1]) * scaled_dxi[1];
    amrex::Real z = (p.pos(2) - plo[2]) * scaled_dxi[2];

    int i = static_cast<int>(std::floor(x));
    int j = static_cast<int>(std::floor(y));
    int k = static_cast<int>(std::floor(z));
    
    amrex::Real wx_hi = x - i;
    amrex::Real wy_hi = y - j;
    amrex::Real wz_hi = z - k;
    
    amrex::Real wx_lo = 1.0 - wx_hi;
    amrex::Real wy_lo = 1.0 - wy_hi;
    amrex::Real wz_lo = 1.0 - wz_hi;

    amrex::Real ls_value = 0;
    ls_value += phi(i,   j  , k  ) * wx_lo * wy_lo * wz_lo;
    ls_value += phi(i+1, j  , k  ) * wx_hi * wy_lo * wz_lo;
    ls_value += phi(i,   j+1, k  ) * wx_lo * wy_hi * wz_lo;
    ls_value += phi(i+1, j+1, k  ) * wx_hi * wy_hi * wz_lo;
    ls_value += phi(i,   j  , k+1) * wx_lo * wy_lo * wz_hi;
    ls_value += phi(i+1, j  , k+1) * wx_hi * wy_lo * wz_hi;
    ls_value += phi(i  , j+1, k+1) * wx_lo * wy_hi * wz_hi;
    ls_value += phi(i+1, j+1, k+1) * wx_hi * wy_hi * wz_hi;
    
    return ls_value;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void level_set_normal (BDEMParticleContainer::ParticleType const& p,
                       const int n_refine,
                       amrex::Real* normal,
                       amrex::Array4<amrex::Real const> const& phi,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
{
    amrex::Real scaled_dxi[3];
    scaled_dxi[0] = n_refine * dxi[0];
    scaled_dxi[1] = n_refine * dxi[1];
    scaled_dxi[2] = n_refine * dxi[2];
    
    amrex::Real x = (p.pos(0) - plo[0]) * scaled_dxi[0];
    amrex::Real y = (p.pos(1) - plo[1]) * scaled_dxi[1];
    amrex::Real z = (p.pos(2) - plo[2]) * scaled_dxi[2];

    int i = static_cast<int>(std::floor(x));
    int j = static_cast<int>(std::floor(y));
    int k = static_cast<int>(std::floor(z));
    
    amrex::Real wx_hi = x - i;
    amrex::Real wy_hi = y - j;
    amrex::Real wz_hi = z - k;
    
    amrex::Real wx_lo = 1.0 - wx_hi;
    amrex::Real wy_lo = 1.0 - wy_hi;
    amrex::Real wz_lo = 1.0 - wz_hi;

    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;

    normal[0] -= phi(i,   j  , k  ) * scaled_dxi[0] * wy_lo * wz_lo;
    normal[0] += phi(i+1, j  , k  ) * scaled_dxi[0] * wy_lo * wz_lo;
    normal[0] -= phi(i,   j+1, k  ) * scaled_dxi[0] * wy_hi * wz_lo;
    normal[0] += phi(i+1, j+1, k  ) * scaled_dxi[0] * wy_hi * wz_lo;
    normal[0] -= phi(i,   j  , k+1) * scaled_dxi[0] * wy_lo * wz_hi;
    normal[0] += phi(i+1, j  , k+1) * scaled_dxi[0] * wy_lo * wz_hi;
    normal[0] -= phi(i  , j+1, k+1) * scaled_dxi[0] * wy_hi * wz_hi;
    normal[0] += phi(i+1, j+1, k+1) * scaled_dxi[0] * wy_hi * wz_hi;

    normal[1] -= phi(i,   j  , k  ) * scaled_dxi[1] * wx_lo * wz_lo;
    normal[1] += phi(i  , j+1, k  ) * scaled_dxi[1] * wx_lo * wz_lo;
    normal[1] -= phi(i+1, j  , k  ) * scaled_dxi[1] * wx_hi * wz_lo;
    normal[1] += phi(i+1, j+1, k  ) * scaled_dxi[1] * wx_hi * wz_lo;
    normal[1] -= phi(i,   j  , k+1) * scaled_dxi[1] * wx_lo * wz_hi;
    normal[1] += phi(i  , j+1, k+1) * scaled_dxi[1] * wx_lo * wz_hi;
    normal[1] -= phi(i+1, j  , k+1) * scaled_dxi[1] * wx_hi * wz_hi;
    normal[1] += phi(i+1, j+1, k+1) * scaled_dxi[1] * wx_hi * wz_hi;

    normal[2] -= phi(i  , j  , k  ) * scaled_dxi[2] * wx_lo * wy_lo;
    normal[2] += phi(i  , j  , k+1) * scaled_dxi[2] * wx_lo * wy_lo;
    normal[2] -= phi(i+1, j  , k  ) * scaled_dxi[2] * wx_hi * wy_lo;
    normal[2] += phi(i+1, j  , k+1) * scaled_dxi[2] * wx_hi * wy_lo;
    normal[2] -= phi(i,   j+1, k  ) * scaled_dxi[2] * wx_lo * wy_hi;
    normal[2] += phi(i  , j+1, k+1) * scaled_dxi[2] * wx_lo * wy_hi;
    normal[2] -= phi(i+1, j+1, k  ) * scaled_dxi[2] * wx_hi * wy_hi;
    normal[2] += phi(i+1, j+1, k+1) * scaled_dxi[2] * wx_hi * wy_hi;
    
    amrex::Real inv_norm = 1.0 / sqrt(dot_product(normal, normal));
    normal[0] *= inv_norm;
    normal[1] *= inv_norm;
    normal[2] *= inv_norm;
}

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_wall_forces(amrex::Real walldist, amrex::Real normal[3],
        BDEMParticleContainer::ParticleType& p,amrex::Real dt,
        amrex::Real fn[3],amrex::Real ft[3], amrex::Real nFr[3])
{
    Real v_rot[3];
    v_rot[0] = walldist * p.rdata(realData::xangvel);
    v_rot[1] = walldist * p.rdata(realData::yangvel);
    v_rot[2] = walldist * p.rdata(realData::zangvel);

    Real vreltrans[3];
    Real cprod[3];
    Real rp = p.rdata(realData::radius);
    Real overlap_n = rp - walldist;

    cross_product(v_rot, normal, cprod);
    vreltrans[0] = p.rdata(realData::xvel) + cprod[0];
    vreltrans[1] = p.rdata(realData::yvel) + cprod[1];
    vreltrans[2] = p.rdata(realData::zvel) + cprod[2];

    Real vreltrans_norm = dot_product(vreltrans, normal);

    Real vrel_t[3];
    vrel_t[0] = vreltrans[0] - vreltrans_norm*normal[0];
    vrel_t[1] = vreltrans[1] - vreltrans_norm*normal[1];
    vrel_t[2] = vreltrans[2] - vreltrans_norm*normal[2];

    int phase = p.idata(intData::phase);

    Real kn_des_w = DEM::k_n_wall;
    Real meff = p.rdata(realData::mass);
    Real etan_des_w = two*sqrt(DEM::k_n_wall*meff)*fabs(log(DEM::e_n_wall));
    etan_des_w /= sqrt(PI*PI + log(DEM::e_n_wall)*log(DEM::e_n_wall));

    Real kt_des_w = DEM::k_t_wall;
    Real Reff = rp;
    Real rad1 = rp;
    Real etat_des_w = two*sqrt(DEM::k_t_wall*meff)*fabs(log(DEM::e_t_wall));
    etat_des_w /= sqrt(PI*PI + log(DEM::e_t_wall)*log(DEM::e_t_wall));

    Real overlap_t[3];
    Real mag_overlap_t;

    // calculate the normal contact force
    fn[0] = -(kn_des_w*overlap_n*normal[0]
            + etan_des_w*vreltrans_norm*normal[0]);
    fn[1] = -(kn_des_w*overlap_n*normal[1]
            + etan_des_w*vreltrans_norm*normal[1]);
    fn[2] = -(kn_des_w*overlap_n*normal[2]
            + etan_des_w*vreltrans_norm*normal[2]);

    // calculate the tangential displacement
    overlap_t[0] = dt*vrel_t[0];
    overlap_t[1] = dt*vrel_t[1];
    overlap_t[2] = dt*vrel_t[2];

    mag_overlap_t = sqrt(dot_product(overlap_t, overlap_t));

    if (mag_overlap_t > zero)
    {
        Real fnmd = DEM::mu_wall * sqrt(dot_product(fn, fn));
        Real ftsp[THREEDIM];
        Real ftn[THREEDIM];
        Real mag_ftsp;
        Real mag_ftn;
        Real tangent[THREEDIM];
        tangent[XDIR] = overlap_t[XDIR]/mag_overlap_t;
        tangent[YDIR] = overlap_t[YDIR]/mag_overlap_t;
        tangent[ZDIR] = overlap_t[ZDIR]/mag_overlap_t;
        // calculate the tangential contact force
        Real overlap_t_tan = dot_product(overlap_t,tangent);
        Real vrel_t_tan = dot_product(vrel_t,tangent);
        // calculate the tangential contact force
        ftsp[XDIR] = -(kt_des_w*overlap_t_tan*tangent[XDIR]
                + etat_des_w*vrel_t_tan*tangent[XDIR]);
        ftsp[YDIR] = -(kt_des_w*overlap_t_tan*tangent[YDIR]
                + etat_des_w*vrel_t_tan*tangent[YDIR]);
        ftsp[ZDIR] = -(kt_des_w*overlap_t_tan*tangent[ZDIR]
                + etat_des_w*vrel_t_tan*tangent[ZDIR]);

        ftn[XDIR] = -fnmd * tangent[XDIR];
        ftn[YDIR] = -fnmd * tangent[YDIR];
        ftn[ZDIR] = -fnmd * tangent[ZDIR];
        mag_ftsp = sqrt(dot_product(ftsp,ftsp));
        mag_ftn = sqrt(dot_product(ftn,ftn));
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

    Real omega_r[THREEDIM];
    omega_r[XDIR] = p.rdata(realData::xangvel);
    omega_r[YDIR] = p.rdata(realData::yangvel);
    omega_r[ZDIR] = p.rdata(realData::zangvel);
    Real mag_omegar = sqrt(dot_product(omega_r,omega_r));
    Real mag_fn = sqrt(dot_product(fn,fn));
    if (mag_omegar > zero)
    {
        Real nomega[THREEDIM];
        Real vRatio[THREEDIM];
        Real Froll[THREEDIM];
        vRatio[XDIR] = omega_r[XDIR]/mag_omegar;
        vRatio[YDIR] = omega_r[YDIR]/mag_omegar;
        vRatio[ZDIR] = omega_r[ZDIR]/mag_omegar;
        Froll[XDIR] = -DEM::muR_wall * mag_fn * vRatio[XDIR];
        Froll[YDIR] = -DEM::muR_wall * mag_fn * vRatio[YDIR];
        Froll[ZDIR] = -DEM::muR_wall * mag_fn * vRatio[ZDIR];
        cross_product(normal, Froll, nFr);
    }
    else
    {
        nFr[XDIR] = zero;
        nFr[YDIR] = zero;
        nFr[ZDIR] = zero;
    }
}


