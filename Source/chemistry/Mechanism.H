#ifndef MECHANISM_H_
#define MECHANISM_H_

#include <AMReX.H>
#include <Constants.H>

AMREX_GPU_DEVICE AMREX_INLINE void getProductionRate(int nspecies,int nsolidspecs,
            int nreac,amrex::Real* spec,amrex::Real* molwts,
            amrex::Real rhos, amrex::Real rad,amrex::Real radinit,
            amrex::Real Temp,
            int* solidspec_ids, 
            amrex::Real* reactmatrix, amrex::Real *arrh_A,
            amrex::Real *arrh_Ea, amrex::Real* wdot)
{

    //note: spec[sp] is concentration in kg/m3

    //zero production rate
    for(int sp=0;sp<(nspecies+1);sp++)
    {
        wdot[sp]=0.0;
    }

    for(int r=0;r<nreac;r++)
    {
        amrex::Real concmult=1.0;
        for(int sp=0;sp<nspecies;sp++)
        {
            if(reactmatrix[r*nspecies+sp] < 0.0)
            {
                amrex::Real stoichcoeff=amrex::Math::abs(reactmatrix[r*nspecies+sp]);
                concmult*=std::pow(spec[sp]/molwts[sp],stoichcoeff);
            }
        }
        amrex::Real rateconst=arrh_A[r]*std::exp(-arrh_Ea[r]/RUNIV/Temp);
        amrex::Real reactrate=rateconst*concmult;

        //update production/removal
        for(int sp=0;sp<nspecies;sp++)
        {
            if(amrex::Math::abs(reactmatrix[r*nspecies+sp]) > 0.0)
            {
                wdot[sp] += reactmatrix[r*nspecies+sp]*reactrate;
            }
        }
    }

    //convert wdot from mol/m3/s to kg/m3/s
    for(int sp=0;sp<nspecies;sp++)
    {
        wdot[sp]*=molwts[sp];
    }

    //radius change
    amrex::Real radsource=0.0;
    for(int i=0;i<nsolidspecs;i++)
    {
        radsource+=wdot[solidspec_ids[i]];
    }

    wdot[nspecies]=radsource/rhos*(radinit*radinit*radinit/(3.0*rad*rad));
}

class mechanism
{
    public:

        int nspecies;
        int nreac;
        amrex::Vector<amrex::Real> reactmatrix;
        amrex::Vector<amrex::Real> arrh_A;
        amrex::Vector<amrex::Real> arrh_Ea;
        amrex::Vector<std::string> specnames;
        amrex::Vector<amrex::Real> molwts;
        amrex::Vector<int> solidspec_ids;

        void define()
        {
            setparams();
        
            reactmatrix.resize(nspecies*nreac,0.0);
            arrh_A.resize(nreac,0.0);
            arrh_Ea.resize(nreac,0.0);
            specnames.resize(nspecies);
            molwts.resize(nspecies);

            solidspec_ids.resize(0);

            setreactmatrix();
            setspecnames();
            setsolidspecs();
            setmolwts();

        }

        virtual void setparams() = 0;
        virtual void setreactmatrix() = 0;
        virtual void setspecnames() = 0;
        virtual void setsolidspecs() = 0;
        virtual void setmolwts() = 0;

};

#endif
