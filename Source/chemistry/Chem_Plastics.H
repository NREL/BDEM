#ifndef CHEMPLASTICS_H_
#define CHEMPLASTICS_H_

#include <Mechanism.H>

class plasticschem:
    public mechanism
{
    private:
        const int PL=0;
        const int TR=1;
        const int GS=2;

    public:

    void setparams()
    {
        nspecies=3;
        nreac=1;

    }

    void setmolwts()
    {
        //in kg/mol
        //this is made up!
        molwts[PL]=100.0;
        molwts[TR]=100.0;
        molwts[GS]=100.0;
    }

    void setreactmatrix()
    {
        int r;
        
        //reaction 1
        r=0;
        reactmatrix[r*nspecies+PL] = -1.0;
        reactmatrix[r*nspecies+TR] =  0.6;
        reactmatrix[r*nspecies+GS] =  0.4;
        arrh_A[r] = 2.3e18;
        arrh_Ea[r] = 285.7e3; //J/mol 

    }

    void setspecnames()
    {
        specnames[PL]="Plastics";
        specnames[TR]="Tar";
        specnames[GS]="Gas";
    }

    void setsolidspecs()
    {
        solidspec_ids.push_back(PL);
    }
}; 

#endif
