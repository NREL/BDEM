#ifndef CHEMNULL_H_
#define CHEMNULL_H_

#include <Mechanism.H>

class nullchem:
    public mechanism
{
    private:
        const int SLD=0;

    public:

    void setparams()
    {
        nspecies=1;
        nreac=0;
    }

    void setreactmatrix()
    {

    }

    void setspecnames()
    {
        specnames[SLD]="Solid";
    }

    void setmolwts()
    {
        //in kg/mol
        molwts[SLD]=1.0;
    }

    void setsolidspecs()
    {
        solidspec_ids.push_back(SLD);
    }
}; 

#endif
