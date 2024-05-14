# BDEM
<div align="center">
<img src="https://github.com/NREL/BDEM/blob/main/images/screw_feeder_discharge_types_noletter_vertical.png" alt="BDEM Logo">
</div>

## A discrete-element method tool for modeling granular flows

BDEM is a massively-parallel discrete-element method (DEM) solver capable of simulating a wide array of granular flow problems. It has been extensively used for biomass material, hence the B in the name. BDEM is implemented using the AMReX framework and can support parallel particle data structires on hybrid CPU/GPU HPC architectures. Our solver can represent highly-variable granular material using models that capture complex particle shapes, moisture content, and particle force interactions. Level-set and STL boundary formulations allow for the inclusing of complex moving geometries. 

<div align="center">
<img src="https://github.com/NREL/BDEM/blob/main/images/BDEM_functionality.png" alt="BDEM Functionality">
</div>

## Models and Features

- Linear spring-dashpot (LSD) and Hertz-Mindlin models for particle contact forces
- Bonded sphere model for representing complex particle shapes
- Liquid bridge representation for moisture-laden particles
- Level-set and STL descriptions of irregular moving geometry
- Reactive particle chemistry and temporally evolving particle sizes
- Parallelization via OpenMPI/MPICH and GPU Acceleration with CUDA (NVidia) and HIP (AMD)
- Parallel I/O
- Plotfile format supported by Amrvis, VisIt, ParaView and yt

# Build instructions
* gcc and an MPI library (openMPI/MPICH) for CPU builds. cuda-11.0 is also required for GPU builds
* This tool depends on the AMReX library (https://github.com/AMReX-Codes/amrex)
* Go to the build folder and make sure the AMREX_HOME variable is set to your AMReX repository
* Build executable using the GNUMakefile (set USE_MPI and USE_CUDA=TRUE/FALSE depending on architecture and desired parallel execution) and run "make"
* Several test cases can be found in the test directory for getting started using the code

# Visualization instructions

* The outputs for a case are in the form of AMReX plotfiles
* These plot files can be open usine AMReX grid reader in ParaView (see https://amrex-codes.github.io/amrex/docs_html/Visualization.html#paraview)
* Alternatively visit can be used after converting to vtp format. see https://amrex-codes.github.io/amrex/docs_html/Visualization_Chapter.html
* We also provide python scripts to access particle data as numpy arrays. see https://github.com/NREL/BDEM/blob/main/tests/read_and_avg_particle_data.py
