#ifndef MOVUTILS_H_
#define MOVUTILS_H_
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include<AMReX_REAL.H>
#include<AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>

void move_this_point(amrex::Real coord[3],amrex::Real newcoord[3],amrex::Real timestep,
                     int movetype,int movedir,amrex::Real center[3],amrex::Real movevel);

#endif
