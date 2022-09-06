#ifndef CHEMPLASTICSMIX_H_
#define CHEMPLASTICSMIX_H_

#include <Mechanism.H>

class plasticsmixchem:
    public mechanism
{
    private:
        const int PE=0;
        const int PET=1;
        const int PS=2;
        const int TR=3;
        const int GS=4;

    public:

    void setparams()
    {
        nspecies=5;
        nreac=3;

    }

    void setmolwts()
    {
        //in kg/mol
        //this is made up!
        molwts[PE]=100.0;
        molwts[PET]=100.0;
        molwts[PS]=100.0;
        molwts[TR]=100.0;
        molwts[GS]=100.0;
    }

    void setreactmatrix()
    {
        int r;
        
        //reaction 1
        r=0;
        reactmatrix[r*nspecies+PE] = -1.0;
        reactmatrix[r*nspecies+TR] =  0.6;
        reactmatrix[r*nspecies+GS] =  0.4;
        arrh_A[r] = 2.3e18;
        arrh_Ea[r] = 285.7e3; //J/mol 
        
        //reaction 2
        r=1;
        reactmatrix[r*nspecies+PET] = -1.0;
        reactmatrix[r*nspecies+TR]  =  0.667;
        reactmatrix[r*nspecies+GS]  =  0.333;
        arrh_A[r] = 3.85e9;
        arrh_Ea[r] = 161.2e3; //J/mol 
        
        //reaction 3
        r=2;
        reactmatrix[r*nspecies+PS] = -1.0;
        reactmatrix[r*nspecies+TR] =  0.938;
        reactmatrix[r*nspecies+GS] =  0.062;
        arrh_A[r] = 1.61e8;
        arrh_Ea[r] = 136.6e3; //J/mol 

    }

    void setspecnames()
    {
        specnames[PE]="PE";
        specnames[PET]="PET";
        specnames[PS]="PS";
        specnames[TR]="Tar";
        specnames[GS]="Gas";
    }

    void setsolidspecs()
    {
        solidspec_ids.push_back(PE);
        solidspec_ids.push_back(PET);
        solidspec_ids.push_back(PS);
    }
}; 

#endif
