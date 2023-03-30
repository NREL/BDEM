#include <BDEM_ParticleContainer.H>

struct bondedParticleData {
  // Bonded particle positions for each particle type assumes center is 
  // positioned at origin, with a constant nominal radius
  amrex::Real Type0_positions[BP_NP0][THREEDIM] =
    {{0.0, 0.0, 0.0}};

  amrex::Real Type1_positions[BP_NP1][THREEDIM] =
    {{2.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {-2.0, 0.0, 0.0}};

  amrex::Real Type2_positions[BP_NP2][THREEDIM] =
    {{4.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {-2.0, 0.0, 0.0}, {-4.0, 0.0, 0.0}};

  amrex::Real Type3_positions[BP_NP3][THREEDIM] =
    {{9.0, 0.0, 0.0}, {7.0, 0.0, 0.0}, {5.0, 0.0, 0.0}, {3.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
     {-1.0, 0.0, 0.0}, {-3.0, 0.0, 0.0}, {-5.0, 0.0, 0.0}, {-7.0, 0.0, 0.0}, {-9.0, 0.0, 0.0}};

  amrex::Real Type4_positions[BP_NP4][THREEDIM] =
    {{2.0, 1.73205080756, 0.0}, {0.0, 1.73205080756, 0.0}, {-2.0, 1.73205080756, 0.0},        // Top row
     {3.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}, {-3.0, 0.0, 0.0},                    // Middle row
     {2.0, -1.73205080756, 0.0}, {0.0, -1.73205080756, 0.0}, {-2.0, -1.73205080756, 0.0}};    // Bottom row

  // Initialize arrays that contain bond data
  int Type0_bonds[BP_NP0][MAXBONDS];
  int Type1_bonds[BP_NP1][MAXBONDS];
  int Type2_bonds[BP_NP2][MAXBONDS];
  int Type3_bonds[BP_NP3][MAXBONDS];
  int Type4_bonds[BP_NP4][MAXBONDS];

  // Constructor to fill in bond info
  bondedParticleData(){
      // Manually adding particle bonds 
      // TODO: Probably a better way to do this...

      // Value of -1 indicates no bond
      for(int i=0; i<BP_NP0; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type0_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP1; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type1_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP2; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type2_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP3; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type3_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP4; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type4_bonds[i][j] = -1;
          }
      }

      // Particle type 1
      Type1_bonds[0][0] = 1;
      Type1_bonds[1][0] = 0; Type1_bonds[1][1] = 2;
      Type1_bonds[2][0] = 1;

      // Particle type 2
      Type2_bonds[0][0] = 1;
      Type2_bonds[1][0] = 0; Type2_bonds[1][1] = 2;
      Type2_bonds[2][0] = 1; Type2_bonds[2][1] = 3;
      Type2_bonds[3][0] = 2; Type2_bonds[3][1] = 4;
      Type2_bonds[4][0] = 3;

      // Particle type 3
      Type3_bonds[0][0] = 1;
      Type3_bonds[1][0] = 0; Type3_bonds[1][1] = 2;
      Type3_bonds[2][0] = 1; Type3_bonds[2][1] = 3;
      Type3_bonds[3][0] = 2; Type3_bonds[3][1] = 4;
      Type3_bonds[4][0] = 3; Type3_bonds[4][1] = 5;
      Type3_bonds[5][0] = 4; Type3_bonds[5][1] = 6;
      Type3_bonds[6][0] = 5; Type3_bonds[6][1] = 7;
      Type3_bonds[7][0] = 6; Type3_bonds[7][1] = 8;
      Type3_bonds[8][0] = 7; Type3_bonds[8][1] = 9;
      Type3_bonds[9][0] = 8;

      // Particle type 4
      Type4_bonds[0][0] = 1; Type4_bonds[0][1] = 3; Type4_bonds[0][2] = 4;
      Type4_bonds[1][0] = 0; Type4_bonds[1][1] = 2; Type4_bonds[1][2] = 4; Type4_bonds[1][3] = 5;
      Type4_bonds[2][0] = 1; Type4_bonds[2][1] = 5; Type4_bonds[2][2] = 6;
      Type4_bonds[3][0] = 0; Type4_bonds[3][1] = 4; Type4_bonds[3][2] = 7;
      Type4_bonds[4][0] = 0; Type4_bonds[4][1] = 1; Type4_bonds[4][2] = 3; Type4_bonds[4][3] = 5; Type4_bonds[4][4] = 7; Type4_bonds[4][5] = 8;
      Type4_bonds[5][0] = 1; Type4_bonds[5][1] = 2; Type4_bonds[5][2] = 4; Type4_bonds[5][3] = 6; Type4_bonds[5][4] = 8; Type4_bonds[5][5] = 9;
      Type4_bonds[6][0] = 2; Type4_bonds[6][1] = 5; Type4_bonds[6][2] = 9;
      Type4_bonds[7][0] = 3; Type4_bonds[7][1] = 4; Type4_bonds[7][2] = 8;
      Type4_bonds[8][0] = 4; Type4_bonds[8][1] = 5; Type4_bonds[8][2] = 7; Type4_bonds[8][3] = 9;
      Type4_bonds[9][0] = 5; Type4_bonds[9][1] = 6; Type4_bonds[9][2] = 8;
  }
};

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_bonded_particle_pos(const bondedParticleData bp_data,
                             int type, int comp, amrex::Real radius, 
                             amrex::Real bp_pos[THREEDIM], 
                             amrex::Real bp_euler_angles[THREEDIM], 
                             amrex::Real pc_pos[THREEDIM]){
    
  amrex::Real posnorm[THREEDIM];

  // Get the nominal origin-centered particle positions
  if(type == 0){
    posnorm[XDIR] = bp_data.Type0_positions[comp][0]; 
    posnorm[YDIR] = bp_data.Type0_positions[comp][1]; 
    posnorm[ZDIR] = bp_data.Type0_positions[comp][2]; 
  } 
  if(type == 1){
    posnorm[XDIR] = bp_data.Type1_positions[comp][0]; 
    posnorm[YDIR] = bp_data.Type1_positions[comp][1]; 
    posnorm[ZDIR] = bp_data.Type1_positions[comp][2]; 
  }
  if(type == 2){
    posnorm[XDIR] = bp_data.Type2_positions[comp][0]; 
    posnorm[YDIR] = bp_data.Type2_positions[comp][1]; 
    posnorm[ZDIR] = bp_data.Type2_positions[comp][2]; 
  }
  if(type == 3){
    posnorm[XDIR] = bp_data.Type3_positions[comp][0]; 
    posnorm[YDIR] = bp_data.Type3_positions[comp][1]; 
    posnorm[ZDIR] = bp_data.Type3_positions[comp][2]; 
  }
  if(type == 4){
    posnorm[XDIR] = bp_data.Type4_positions[comp][0]; 
    posnorm[YDIR] = bp_data.Type4_positions[comp][1]; 
    posnorm[ZDIR] = bp_data.Type4_positions[comp][2]; 
  }

  // Scale and translate particle position using radius and bp center
  // TODO: Use euler angle data to rotate bonded particles
  pc_pos[XDIR] = posnorm[XDIR]*radius + bp_pos[XDIR];
  pc_pos[YDIR] = posnorm[YDIR]*radius + bp_pos[YDIR];
  pc_pos[ZDIR] = posnorm[ZDIR]*radius + bp_pos[ZDIR];
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void bp_init(BDEMParticleContainer::ParticleType& p,
             const bondedParticleData bp_data, 
             int phase, amrex::Real ppos[THREEDIM],
             amrex::Real radius, amrex::Real density, 
             amrex::Real vel[THREEDIM], amrex::Real temperature,
             int comp, int type, int* ids){

  // Set unique ID and assign CPU
  p.cpu() = ParallelDescriptor::MyProc();
  
  // Set particle quantities that are read in from particle file
  p.idata(intData::phase) = phase;
  p.pos(0) = ppos[XDIR];
  p.pos(1) = ppos[YDIR];
  p.pos(2) = ppos[ZDIR];
  p.rdata(realData::radius) = radius;
  p.rdata(realData::density) = density;
  p.rdata(realData::xvel) = vel[XDIR];
  p.rdata(realData::yvel) = vel[YDIR];
  p.rdata(realData::zvel) = vel[ZDIR];
  p.rdata(realData::temperature) = temperature; 

  // Set rotation-related quantities to zero
  p.idata(intData::num_comp_sphere) = 1;
  p.rdata(realData::euler_angle_x) = zero;
  p.rdata(realData::euler_angle_y) = zero;
  p.rdata(realData::euler_angle_z) = zero;
  p.rdata(realData::q0) = 1.0;
  p.rdata(realData::q1) = zero;
  p.rdata(realData::q2) = zero;
  p.rdata(realData::q3) = zero;
  p.rdata(realData::pax) = zero;
  p.rdata(realData::pay) = zero;
  p.rdata(realData::paz) = zero;

  //set initial radius
  p.rdata(realData::radinit)=p.rdata(realData::radius);
  
  // Set other particle properties
  p.rdata(realData::volume)      = fourbythree*PI*pow(p.rdata(realData::radius),three)*p.idata(intData::num_comp_sphere);
  p.rdata(realData::mass)        = p.rdata(realData::density)*p.rdata(realData::volume);
  p.rdata(realData::Iinv)        = 2.5/(p.rdata(realData::mass)*pow(p.rdata(realData::radius),two));
  
  // Assuming no initial component sphere rotation
  p.rdata(realData::xangvel)     = zero;
  p.rdata(realData::yangvel)     = zero;
  p.rdata(realData::zangvel)     = zero;
  
  p.rdata(realData::Ixinv) = 2.5/(p.rdata(realData::mass)*pow(p.rdata(realData::radius),two));
  p.rdata(realData::Iyinv) = 2.5/(p.rdata(realData::mass)*pow(p.rdata(realData::radius),two));
  p.rdata(realData::Izinv) = 2.5/(p.rdata(realData::mass)*pow(p.rdata(realData::radius),two));
  
  p.rdata(realData::fx) = zero;
  p.rdata(realData::fy) = zero;
  p.rdata(realData::fz) = zero;
  p.rdata(realData::taux) = zero;
  p.rdata(realData::tauy) = zero;
  p.rdata(realData::tauz) = zero;  
  p.rdata(realData::fx_bond) = zero;
  p.rdata(realData::fy_bond) = zero;
  p.rdata(realData::fz_bond) = zero;
  p.rdata(realData::taux_bond) = zero;
  p.rdata(realData::tauy_bond) = zero;
  p.rdata(realData::tauz_bond) = zero;  
  p.rdata(realData::theta_x) = zero;  

  // Set bridge indices to -1 to indicate no existing bridges
  for(int br=0; br<MAXBRIDGES; br++){
      p.idata(intData::first_bridge+3*br) = -1;
      p.idata(intData::first_bridge+3*br+1) = -1;
      p.idata(intData::first_bridge+3*br+2) = -1;
  }
  p.rdata(realData::liquid_volume) = zero;
  p.rdata(realData::total_bridge_volume) = zero;

  // Fill in the particle each particle ID that p shares a bond with
  if(type == 0){
      for(int b=0; b<MAXBONDS; b++) p.idata(intData::first_bond + b) = (bp_data.Type0_bonds[comp][b] != -1) ? ids[bp_data.Type0_bonds[comp][b]]:0;
  }
  if(type == 1){
      for(int b=0; b<MAXBONDS; b++) p.idata(intData::first_bond + b) = (bp_data.Type1_bonds[comp][b] != -1) ? ids[bp_data.Type1_bonds[comp][b]]:0;
  }
  if(type == 2){
      for(int b=0; b<MAXBONDS; b++) p.idata(intData::first_bond + b) = (bp_data.Type2_bonds[comp][b] != -1) ? ids[bp_data.Type2_bonds[comp][b]]:0;
  }
  if(type == 3){
      for(int b=0; b<MAXBONDS; b++) p.idata(intData::first_bond + b) = (bp_data.Type3_bonds[comp][b] != -1) ? ids[bp_data.Type3_bonds[comp][b]]:0;
  }
  if(type == 4){
      for(int b=0; b<MAXBONDS; b++) p.idata(intData::first_bond + b) = (bp_data.Type4_bonds[comp][b] != -1) ? ids[bp_data.Type4_bonds[comp][b]]:0;
  }
  
  //FIXME: get chemistry data from inputs file
  for(int sp=0;sp<MAXSPECIES;sp++)
  {
      p.rdata(realData::firstspec+sp)=0.0;
  }
}
