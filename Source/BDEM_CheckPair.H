#ifndef BDEM_CHECKPAIR_H_
#define BDEM_CHECKPAIR_H_

#include <BDEM_Specs.H>

struct CheckPair
{
    // FIXME: Make sure check pair accounts for liquid bridges
    // FIXME: Neighbor detection seems to fail when there is a large particle radius ratio, or small particle radii used
    // TODO: Find more efficient check for neighbors, updated code is more robust but also substantially slower
    template <class P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    bool operator()(const P& p1, const P& p2) const
    {
        bool not_pair = false;

        Real max_len; 
        max_len = p1.rdata(realData::radius) + p2.rdata(realData::radius);
        amrex::Real d0 = (p1.pos(0) - p2.pos(0));
        amrex::Real d1 = (p1.pos(1) - p2.pos(1));
        amrex::Real d2 = (p1.pos(2) - p2.pos(2));
        amrex::Real d = sqrt(d0*d0 + d1*d1 + d2*d2);
        return (d <= 1.25*max_len);
        // return (d <= 1.1*max_len);

        // If we get through all component particles with no pair match, return false
        return not_pair;
    }
};

#endif
