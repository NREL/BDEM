#include <BDEM_ParticleContainer.H>

void bp_init(BDEMParticleContainer::ParticleType& p,
             const ParticleTypeData bp_data, 
             int phase, amrex::Real ppos[THREEDIM],
             amrex::Real radius, amrex::Real density, 
             amrex::Real vel[THREEDIM], amrex::Real temperature,
             int comp, int type, int* ids){

  // Set unique ID and assign CPU
  p.cpu() = ParallelDescriptor::MyProc();
  
  // Set particle quantities that are read in from particle file
  p.idata(intData::phase) = phase;
  p.idata(intData::type_id) = type;
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
  if(type == 5){
      for(int b=0; b<MAXBONDS; b++) p.idata(intData::first_bond + b) = (bp_data.Type5_bonds[comp][b] != -1) ? ids[bp_data.Type5_bonds[comp][b]]:0;
  }
  if(type == 6){
      for(int b=0; b<MAXBONDS; b++) p.idata(intData::first_bond + b) = (bp_data.Type6_bonds[comp][b] != -1) ? ids[bp_data.Type6_bonds[comp][b]]:0;
  }
  if(type == 7){
      for(int b=0; b<MAXBONDS; b++) p.idata(intData::first_bond + b) = (bp_data.Type7_bonds[comp][b] != -1) ? ids[bp_data.Type7_bonds[comp][b]]:0;
  }
  if(type == 8){
      for(int b=0; b<MAXBONDS; b++) p.idata(intData::first_bond + b) = (bp_data.Type8_bonds[comp][b] != -1) ? ids[bp_data.Type8_bonds[comp][b]]:0;
  }
  if(type == 9){
      for(int b=0; b<MAXBONDS; b++) p.idata(intData::first_bond + b) = (bp_data.Type9_bonds[comp][b] != -1) ? ids[bp_data.Type9_bonds[comp][b]]:0;
  }
  if(type == 10){
      for(int b=0; b<MAXBONDS; b++) p.idata(intData::first_bond + b) = (bp_data.Type10_bonds[comp][b] != -1) ? ids[bp_data.Type10_bonds[comp][b]]:0;
  }
  
  //FIXME: get chemistry data from inputs file
  for(int sp=0;sp<MAXSPECIES;sp++)
  {
      p.rdata(realData::firstspec+sp)=0.0;
  }
}

