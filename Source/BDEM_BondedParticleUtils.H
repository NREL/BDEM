#include <BDEM_ParticleContainer.H>

struct bondedParticleData {
  // Bonded particle positions for each particle type assumes center is 
  // positioned at origin, with a constant nominal radius

  // TODO: Find better/automated way of setting particle position and bond index info

  /* --------------------------------------------------
   *      ROD-LIKE PARTICLE TYPES
   * -------------------------------------------------- */
  amrex::Real Type0_positions[BP_NP0][THREEDIM] =
    {{0.0, 0.0, 0.0}};

  amrex::Real Type1_positions[BP_NP1][THREEDIM] =
    {{2.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {-2.0, 0.0, 0.0}};

  amrex::Real Type2_positions[BP_NP2][THREEDIM] =
    {{3.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}, {-3.0, 0.0, 0.0}};

  amrex::Real Type3_positions[BP_NP3][THREEDIM] =
    {{4.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {-2.0, 0.0, 0.0}, {-4.0, 0.0, 0.0}};

  amrex::Real Type4_positions[BP_NP4][THREEDIM] =
    {{9.0, 0.0, 0.0}, {7.0, 0.0, 0.0}, {5.0, 0.0, 0.0}, {3.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
     {-1.0, 0.0, 0.0}, {-3.0, 0.0, 0.0}, {-5.0, 0.0, 0.0}, {-7.0, 0.0, 0.0}, {-9.0, 0.0, 0.0}};


  /* --------------------------------------------------
   *      RECTANGULAR PARTICLE TYPES
   * -------------------------------------------------- */

  amrex::Real Type5_positions[BP_NP5][THREEDIM] =
    {{2.0, sqthree, 0.0}, {0.0, sqthree, 0.0}, {-2.0, sqthree, 0.0},        // Top row
     {3.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}, {-3.0, 0.0, 0.0},                    // Middle row
     {2.0, -sqthree, 0.0}, {0.0, -sqthree, 0.0}, {-2.0, -sqthree, 0.0}};    // Bottom row

  // 8x2x1 rectangle
  amrex::Real Type6_positions[BP_NP6][THREEDIM] =
    {{7.0,1.0,0.0},{5.0,1.0,0.0},{3.0,1.0,0.0},{1.0,1.0,0.0},{-1.0,1.0,0.0},{-3.0,1.0,0.0},{-5.0,1.0,0.0},{-7.0,1.0,0.0},
     {7.0,-1.0,0.0},{5.0,-1.0,0.0},{3.0,-1.0,0.0},{1.0,-1.0,0.0},{-1.0,-1.0,0.0},{-3.0,-1.0,0.0},{-5.0,-1.0,0.0},{-7.0,-1.0,0.0}};

  // 8x2x2 rectangle
  amrex::Real Type7_positions[BP_NP7][THREEDIM] =
    {{7.0,1.0,1.0},{5.0,1.0,1.0},{3.0,1.0,1.0},{1.0,1.0,1.0},{-1.0,1.0,1.0},{-3.0,1.0,1.0},{-5.0,1.0,1.0},{-7.0,1.0,1.0},
     {7.0,-1.0,1.0},{5.0,-1.0,1.0},{3.0,-1.0,1.0},{1.0,-1.0,1.0},{-1.0,-1.0,1.0},{-3.0,-1.0,1.0},{-5.0,-1.0,1.0},{-7.0,-1.0,1.0},
     {7.0,1.0,-1.0},{5.0,1.0,-1.0},{3.0,1.0,-1.0},{1.0,1.0,-1.0},{-1.0,1.0,-1.0},{-3.0,1.0,-1.0},{-5.0,1.0,-1.0},{-7.0,1.0,-1.0},
     {7.0,-1.0,-1.0},{5.0,-1.0,-1.0},{3.0,-1.0,-1.0},{1.0,-1.0,-1.0},{-1.0,-1.0,-1.0},{-3.0,-1.0,-1.0},{-5.0,-1.0,-1.0},{-7.0,-1.0,-1.0}
    };

  // 12x3x1 rectangle
  amrex::Real Type8_positions[BP_NP8][THREEDIM] =
    {{11.0,2.0,0.0},{9.0,2.0,0.0},{7.0,2.0,0.0},{5.0,2.0,0.0},{3.0,2.0,0.0},{1.0,2.0,0.0},
     {-1.0,2.0,0.0},{-3.0,2.0,0.0},{-5.0,2.0,0.0},{-7.0,2.0,0.0},{-9.0,2.0,0.0},{-11.0,2.0,0.0},
     {11.0,0.0,0.0},{9.0,0.0,0.0},{7.0,0.0,0.0},{5.0,0.0,0.0},{3.0,0.0,0.0},{1.0,0.0,0.0},
     {-1.0,0.0,0.0},{-3.0,0.0,0.0},{-5.0,0.0,0.0},{-7.0,0.0,0.0},{-9.0,0.0,0.0},{-11.0,0.0,0.0},
     {11.0,-2.0,0.0},{9.0,-2.0,0.0},{7.0,-2.0,0.0},{5.0,-2.0,0.0},{3.0,-2.0,0.0},{1.0,-2.0,0.0},
     {-1.0,-2.0,0.0},{-3.0,-2.0,0.0},{-5.0,-2.0,0.0},{-7.0,-2.0,0.0},{-9.0,-2.0,0.0},{-11.0,-2.0,0.0}
    };

  // 12x3x2 rectangle
  amrex::Real Type9_positions[BP_NP9][THREEDIM] =
    {
     {11.0,2.0,1.0},{9.0,2.0,1.0},{7.0,2.0,1.0},{5.0,2.0,1.0},{3.0,2.0,1.0},{1.0,2.0,1.0},
     {-1.0,2.0,1.0},{-3.0,2.0,1.0},{-5.0,2.0,1.0},{-7.0,2.0,1.0},{-9.0,2.0,1.0},{-11.0,2.0,1.0},
     {11.0,0.0,1.0},{9.0,0.0,1.0},{7.0,0.0,1.0},{5.0,0.0,1.0},{3.0,0.0,1.0},{1.0,0.0,1.0},
     {-1.0,0.0,1.0},{-3.0,0.0,1.0},{-5.0,0.0,1.0},{-7.0,0.0,1.0},{-9.0,0.0,1.0},{-11.0,0.0,1.0},
     {11.0,-2.0,1.0},{9.0,-2.0,1.0},{7.0,-2.0,1.0},{5.0,-2.0,1.0},{3.0,-2.0,1.0},{1.0,-2.0,1.0},
     {-1.0,-2.0,1.0},{-3.0,-2.0,1.0},{-5.0,-2.0,1.0},{-7.0,-2.0,1.0},{-9.0,-2.0,1.0},{-11.0,-2.0,1.0},
     {11.0,2.0,-1.0},{9.0,2.0,-1.0},{7.0,2.0,-1.0},{5.0,2.0,-1.0},{3.0,2.0,-1.0},{1.0,2.0,-1.0},
     {-1.0,2.0,-1.0},{-3.0,2.0,-1.0},{-5.0,2.0,-1.0},{-7.0,2.0,-1.0},{-9.0,2.0,-1.0},{-11.0,2.0,-1.0},
     {11.0,0.0,-1.0},{9.0,0.0,-1.0},{7.0,0.0,-1.0},{5.0,0.0,-1.0},{3.0,0.0,-1.0},{1.0,0.0,-1.0},
     {-1.0,0.0,-1.0},{-3.0,0.0,-1.0},{-5.0,0.0,-1.0},{-7.0,0.0,-1.0},{-9.0,0.0,-1.0},{-11.0,0.0,-1.0},
     {11.0,-2.0,-1.0},{9.0,-2.0,-1.0},{7.0,-2.0,-1.0},{5.0,-2.0,-1.0},{3.0,-2.0,-1.0},{1.0,-2.0,-1.0},
     {-1.0,-2.0,-1.0},{-3.0,-2.0,-1.0},{-5.0,-2.0,-1.0},{-7.0,-2.0,-1.0},{-9.0,-2.0,-1.0},{-11.0,-2.0,-1.0}
    };

  // 12x3x3 offset rectangle
  amrex::Real Type10_positions[BP_NP10][THREEDIM] =
    {
     {10.0,1.0,sqtwo},{8.0,1.0,sqtwo},{6.0,1.0,sqtwo},{4.0,1.0,sqtwo},{2.0,1.0,sqtwo},{0.0,1.0,sqtwo},
     {-2.0,1.0,sqtwo},{-4.0,1.0,sqtwo},{-6.0,1.0,sqtwo},{-8.0,1.0,sqtwo},{-10.0,1.0,sqtwo},
     {10.0,-1.0,sqtwo},{8.0,-1.0,sqtwo},{6.0,-1.0,sqtwo},{4.0,-1.0,sqtwo},{2.0,-1.0,sqtwo},{0.0,-1.0,sqtwo},
     {-2.0,-1.0,sqtwo},{-4.0,-1.0,sqtwo},{-6.0,-1.0,sqtwo},{-8.0,-1.0,sqtwo},{-10.0,-1.0,sqtwo},
      
     {11.0,2.0,0.0},{9.0,2.0,0.0},{7.0,2.0,0.0},{5.0,2.0,0.0},{3.0,2.0,0.0},{1.0,2.0,0.0},
     {-1.0,2.0,0.0},{-3.0,2.0,0.0},{-5.0,2.0,0.0},{-7.0,2.0,0.0},{-9.0,2.0,0.0},{-11.0,2.0,0.0},
     {11.0,0.0,0.0},{9.0,0.0,0.0},{7.0,0.0,0.0},{5.0,0.0,0.0},{3.0,0.0,0.0},{1.0,0.0,0.0},
     {-1.0,0.0,0.0},{-3.0,0.0,0.0},{-5.0,0.0,0.0},{-7.0,0.0,0.0},{-9.0,0.0,0.0},{-11.0,0.0,0.0},
     {11.0,-2.0,0.0},{9.0,-2.0,0.0},{7.0,-2.0,0.0},{5.0,-2.0,0.0},{3.0,-2.0,0.0},{1.0,-2.0,0.0},
     {-1.0,-2.0,0.0},{-3.0,-2.0,0.0},{-5.0,-2.0,0.0},{-7.0,-2.0,0.0},{-9.0,-2.0,0.0},{-11.0,-2.0,0.0},

     {10.0,1.0,-sqtwo},{8.0,1.0,-sqtwo},{6.0,1.0,-sqtwo},{4.0,1.0,-sqtwo},{2.0,1.0,-sqtwo},{0.0,1.0,-sqtwo},
     {-2.0,1.0,-sqtwo},{-4.0,1.0,-sqtwo},{-6.0,1.0,-sqtwo},{-8.0,1.0,-sqtwo},{-10.0,1.0,-sqtwo},
     {10.0,-1.0,-sqtwo},{8.0,-1.0,-sqtwo},{6.0,-1.0,-sqtwo},{4.0,-1.0,-sqtwo},{2.0,-1.0,-sqtwo},{0.0,-1.0,-sqtwo},
     {-2.0,-1.0,-sqtwo},{-4.0,-1.0,-sqtwo},{-6.0,-1.0,-sqtwo},{-8.0,-1.0,-sqtwo},{-10.0,-1.0,-sqtwo}
    };

  // Initialize arrays that contain bond data
  int Type0_bonds[BP_NP0][MAXBONDS];
  int Type1_bonds[BP_NP1][MAXBONDS];
  int Type2_bonds[BP_NP2][MAXBONDS];
  int Type3_bonds[BP_NP3][MAXBONDS];
  int Type4_bonds[BP_NP4][MAXBONDS];
  int Type5_bonds[BP_NP5][MAXBONDS];
  int Type6_bonds[BP_NP6][MAXBONDS];
  int Type7_bonds[BP_NP7][MAXBONDS];
  int Type8_bonds[BP_NP8][MAXBONDS];
  int Type9_bonds[BP_NP9][MAXBONDS];
  int Type10_bonds[BP_NP10][MAXBONDS];

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
      for(int i=0; i<BP_NP5; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type5_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP6; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type6_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP7; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type7_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP8; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type8_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP9; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type9_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP10; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type10_bonds[i][j] = -1;
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
      Type2_bonds[3][0] = 2;

      // Particle type 3
      Type3_bonds[0][0] = 1;
      Type3_bonds[1][0] = 0; Type3_bonds[1][1] = 2;
      Type3_bonds[2][0] = 1; Type3_bonds[2][1] = 3;
      Type3_bonds[3][0] = 2; Type3_bonds[3][1] = 4;
      Type3_bonds[4][0] = 3;

      // Particle type 4
      Type4_bonds[0][0] = 1;
      Type4_bonds[1][0] = 0; Type4_bonds[1][1] = 2;
      Type4_bonds[2][0] = 1; Type4_bonds[2][1] = 3;
      Type4_bonds[3][0] = 2; Type4_bonds[3][1] = 4;
      Type4_bonds[4][0] = 3; Type4_bonds[4][1] = 5;
      Type4_bonds[5][0] = 4; Type4_bonds[5][1] = 6;
      Type4_bonds[6][0] = 5; Type4_bonds[6][1] = 7;
      Type4_bonds[7][0] = 6; Type4_bonds[7][1] = 8;
      Type4_bonds[8][0] = 7; Type4_bonds[8][1] = 9;
      Type4_bonds[9][0] = 8;

      // Particle type 5
      Type5_bonds[0][0] = 1; Type5_bonds[0][1] = 3; Type5_bonds[0][2] = 4;
      Type5_bonds[1][0] = 0; Type5_bonds[1][1] = 2; Type5_bonds[1][2] = 4; Type5_bonds[1][3] = 5;
      Type5_bonds[2][0] = 1; Type5_bonds[2][1] = 5; Type5_bonds[2][2] = 6;
      Type5_bonds[3][0] = 0; Type5_bonds[3][1] = 4; Type5_bonds[3][2] = 7;
      Type5_bonds[4][0] = 0; Type5_bonds[4][1] = 1; Type5_bonds[4][2] = 3; Type5_bonds[4][3] = 5; Type5_bonds[4][4] = 7; Type5_bonds[4][5] = 8;
      Type5_bonds[5][0] = 1; Type5_bonds[5][1] = 2; Type5_bonds[5][2] = 4; Type5_bonds[5][3] = 6; Type5_bonds[5][4] = 8; Type5_bonds[5][5] = 9;
      Type5_bonds[6][0] = 2; Type5_bonds[6][1] = 5; Type5_bonds[6][2] = 9;
      Type5_bonds[7][0] = 3; Type5_bonds[7][1] = 4; Type5_bonds[7][2] = 8;
      Type5_bonds[8][0] = 4; Type5_bonds[8][1] = 5; Type5_bonds[8][2] = 7; Type5_bonds[8][3] = 9;
      Type5_bonds[9][0] = 5; Type5_bonds[9][1] = 6; Type5_bonds[9][2] = 8;

      // Particle type 6
      Type6_bonds[0][0] = 1; Type6_bonds[0][1] = 8; 
      Type6_bonds[1][0] = 0; Type6_bonds[1][1] = 2; Type6_bonds[1][2] = 9;  
      Type6_bonds[2][0] = 1; Type6_bonds[2][1] = 3; Type6_bonds[2][2] = 10;  
      Type6_bonds[3][0] = 2; Type6_bonds[3][1] = 4; Type6_bonds[3][2] = 11;
      Type6_bonds[4][0] = 3; Type6_bonds[4][1] = 5; Type6_bonds[4][2] = 12;  
      Type6_bonds[5][0] = 4; Type6_bonds[5][1] = 6; Type6_bonds[5][2] = 13;  
      Type6_bonds[6][0] = 5; Type6_bonds[6][1] = 7; Type6_bonds[6][2] = 14;  
      Type6_bonds[7][0] = 6; Type6_bonds[7][1] = 15; 
      Type6_bonds[8][0] = 0; Type6_bonds[8][1] = 9; 
      Type6_bonds[9][0] = 1; Type6_bonds[9][1] = 8; Type6_bonds[9][2] = 10;  
      Type6_bonds[10][0] = 2; Type6_bonds[10][1] = 9; Type6_bonds[10][2] = 11;  
      Type6_bonds[11][0] = 3; Type6_bonds[11][1] = 10; Type6_bonds[11][2] = 12;  
      Type6_bonds[12][0] = 4; Type6_bonds[12][1] = 11; Type6_bonds[12][2] = 13;  
      Type6_bonds[13][0] = 5; Type6_bonds[13][1] = 12; Type6_bonds[13][2] = 14;  
      Type6_bonds[14][0] = 6; Type6_bonds[14][1] = 13; Type6_bonds[14][2] = 15;  
      Type6_bonds[15][0] = 7; Type6_bonds[15][1] = 14;

      // Particle type 7
      Type7_bonds[0][0] = 1;  Type7_bonds[0][1] = 8;  Type7_bonds[0][2] = 16;
      Type7_bonds[7][0] = 6;  Type7_bonds[7][1] = 15;  Type7_bonds[7][2] = 23;
      Type7_bonds[8][0] = 0;  Type7_bonds[8][1] = 9;  Type7_bonds[8][2] = 24;
      Type7_bonds[15][0] = 7; Type7_bonds[15][1] = 14; Type7_bonds[15][2] = 31; 
      Type7_bonds[16][0] = 0; Type7_bonds[16][1] = 17; Type7_bonds[16][2] = 24;    
      Type7_bonds[23][0] = 7; Type7_bonds[23][1] = 22; Type7_bonds[23][2] = 31;
      Type7_bonds[24][0] = 8; Type7_bonds[24][1] = 16; Type7_bonds[24][2] = 25;   
      Type7_bonds[31][0] = 15; Type7_bonds[31][1] = 23; Type7_bonds[31][2] = 30;

      for(int i=0; i<BP_NP8; i++){
        if((i>0 && i<7) || (i>16 && i<23)){ // Top layer (y-dir)
          Type7_bonds[i][0] = i-1; Type7_bonds[i][1] = i+1; Type7_bonds[i][2] = i+8; Type7_bonds[i][3] = (i+16) % BP_NP7; 
        }      
        if((i>8 && i<15) || (i>24 && i<31)){ // Bottom layer (y-dir)
          Type7_bonds[i][0] = i-8; Type7_bonds[i][1] = i-1; Type7_bonds[i][2] = i+1; Type7_bonds[i][3] = (i+16) % BP_NP7; 
        }      
      } 

      // Particle type 8
      // Manually adding all corner point bonds and middle layer y-edge point bonds
      Type8_bonds[0][0] = 1;  Type8_bonds[0][1] = 12; 
      Type8_bonds[11][0] = 10; Type8_bonds[11][1] = 23;   
      Type8_bonds[12][0] = 0; Type8_bonds[12][1] = 13; Type8_bonds[12][2] = 24;
      Type8_bonds[23][0] = 11; Type8_bonds[23][1] = 22; Type8_bonds[23][2] = 35; 
      Type8_bonds[24][0] = 12; Type8_bonds[24][1] = 25; 
      Type8_bonds[35][0] = 23; Type8_bonds[35][1] = 34;

      for(int i=0; i<BP_NP8; i++){
        if(i>0 && i<11){ // Top layer (y-dir) 
          Type8_bonds[i][0] = i-1; Type8_bonds[i][1] = i+1; Type8_bonds[i][2] = i+12;  
        }
        if(i>12 && i<23){ // Middle layer (y-dir)
          Type8_bonds[i][0] = i-12; Type8_bonds[i][1] = i-1; Type8_bonds[i][2] = i+1; Type8_bonds[i][3] = i+12; 
        }
        if(i>24 && i<35){ // Bottom layer (y-dir) 
          Type8_bonds[i][0] = i-12; Type8_bonds[i][1] = i-1; Type8_bonds[i][2] = i+1; 
        }
      }

      // Particle type 9
      // Manually adding all corner point bonds and middle layer y-edge point bonds
      Type9_bonds[0][0] = 1;  Type9_bonds[0][1] = 12;  Type9_bonds[0][2] = 36;
      Type9_bonds[11][0] = 10; Type9_bonds[11][1] = 23; Type9_bonds[11][2] = 47;
      Type9_bonds[12][0] = 0; Type9_bonds[12][1] = 13; Type9_bonds[12][2] = 24; Type9_bonds[12][3] = 48;
      Type9_bonds[23][0] = 11; Type9_bonds[23][1] = 22; Type9_bonds[23][2] = 35; Type9_bonds[23][3] = 59;
      Type9_bonds[24][0] = 12; Type9_bonds[24][1] = 25; Type9_bonds[24][2] = 60;
      Type9_bonds[35][0] = 23; Type9_bonds[35][1] = 34; Type9_bonds[35][2] = 71;
      Type9_bonds[36][0] = 0; Type9_bonds[36][1] = 37; Type9_bonds[36][2] = 48;
      Type9_bonds[47][0] = 11; Type9_bonds[47][1] = 46; Type9_bonds[47][2] = 59;
      Type9_bonds[48][0] = 12; Type9_bonds[48][1] = 36; Type9_bonds[48][2] = 49; Type9_bonds[48][3] = 60;
      Type9_bonds[59][0] = 23; Type9_bonds[59][1] = 47; Type9_bonds[59][2] = 58; Type9_bonds[59][3] = 71;
      Type9_bonds[60][0] = 24; Type9_bonds[60][1] = 48; Type9_bonds[60][2] = 61;
      Type9_bonds[71][0] = 35; Type9_bonds[71][1] = 59; Type9_bonds[71][2] = 70;

      for(int i=0; i<BP_NP9; i++){
        if((i>0 && i<11) || (i>36 && i<47)){ // Top layer (y-dir) 
          Type9_bonds[i][0] = i-1; Type9_bonds[i][1] = i+1; Type9_bonds[i][2] = i+12; Type9_bonds[i][3] = (i+36) % BP_NP9; 
        }
        if((i>12 && i<23) || (i>48 && i<59)){ // Middle layer (y-dir)
          Type9_bonds[i][0] = i-12; Type9_bonds[i][1] = i-1; Type9_bonds[i][2] = i+1; Type9_bonds[i][3] = i+12; Type9_bonds[i][4] = (i+36) % BP_NP9; 
        }
        if((i>24 && i<35) || (i>60 && i<71)){ // Bottom layer (y-dir) 
          Type9_bonds[i][0] = i-12; Type9_bonds[i][1] = i-1; Type9_bonds[i][2] = i+1; Type9_bonds[i][3] = (i+36) % BP_NP9; 
        }
      }

      // Particle type 10
      // Manually adding all corner point bonds and middle layer y-edge point bonds
      Type10_bonds[0][0] = 1;  Type10_bonds[0][1] = 11;  Type10_bonds[0][2] = 22;  Type10_bonds[0][3] = 23;  Type10_bonds[0][4] = 34;  Type10_bonds[0][5] = 35;  
      Type10_bonds[10][0] = 9; Type10_bonds[10][1] = 21; Type10_bonds[10][2] = 32; Type10_bonds[10][3] = 33; Type10_bonds[10][4] = 44; Type10_bonds[10][5] = 45; 
      Type10_bonds[11][0] = 0; Type10_bonds[11][1] = 12; Type10_bonds[11][2] = 34; Type10_bonds[11][3] = 35; Type10_bonds[11][4] = 46; Type10_bonds[11][5] = 47; 
      Type10_bonds[21][0] = 10; Type10_bonds[21][1] = 20; Type10_bonds[21][2] = 44; Type10_bonds[21][3] = 45; Type10_bonds[21][4] = 56; Type10_bonds[21][5] = 57; 
      Type10_bonds[22][0] = 0; Type10_bonds[22][1] = 23; Type10_bonds[22][2] = 34; Type10_bonds[22][3] = 58; 
      Type10_bonds[33][0] = 10; Type10_bonds[33][1] = 32; Type10_bonds[33][2] = 45; Type10_bonds[33][3] = 68; 
      Type10_bonds[46][0] = 11; Type10_bonds[46][1] = 34; Type10_bonds[46][2] = 47; Type10_bonds[46][3] = 69; 
      Type10_bonds[57][0] = 21; Type10_bonds[57][1] = 45; Type10_bonds[57][2] = 56; Type10_bonds[57][3] = 79; 
      Type10_bonds[58][0] = 22; Type10_bonds[58][1] = 23; Type10_bonds[58][2] = 34; Type10_bonds[58][3] = 35; Type10_bonds[58][4] = 59; Type10_bonds[58][5] = 69; 
      Type10_bonds[68][0] = 32; Type10_bonds[68][1] = 33; Type10_bonds[68][2] = 44; Type10_bonds[68][3] = 45; Type10_bonds[68][4] = 67; Type10_bonds[68][5] = 79; 
      Type10_bonds[69][0] = 34; Type10_bonds[69][1] = 35; Type10_bonds[69][2] = 46; Type10_bonds[69][3] = 47; Type10_bonds[69][4] = 58; Type10_bonds[69][5] = 70; 
      Type10_bonds[79][0] = 44; Type10_bonds[79][1] = 45; Type10_bonds[79][2] = 56; Type10_bonds[79][3] = 57; Type10_bonds[79][4] = 68; Type10_bonds[79][5] = 78; 
      Type10_bonds[34][0] = 0; Type10_bonds[34][1] = 11; Type10_bonds[34][2] = 22; Type10_bonds[34][3] = 35; Type10_bonds[34][4] = 46; Type10_bonds[34][5] = 58; 
      Type10_bonds[34][6] = 69; 
      Type10_bonds[45][0] = 10; Type10_bonds[45][1] = 21; Type10_bonds[45][2] = 33; Type10_bonds[45][3] = 44; Type10_bonds[45][4] = 57; Type10_bonds[45][5] = 68; 
      Type10_bonds[45][6] = 79; 
      for(int i=0; i<BP_NP10; i++){
        // Top layer (z-dir) interior bonds
        if((i>0 && i<10)){ 
          Type10_bonds[i][0] = i-1; Type10_bonds[i][1] = i+1; Type10_bonds[i][2] = i+11; 
          Type10_bonds[i][3] = i+22; Type10_bonds[i][4] = i+23; Type10_bonds[i][5] = i+34; Type10_bonds[i][6] = i+35;
        } 
        if (i>11 || i<21){
          Type10_bonds[i][0] = i-11; Type10_bonds[i][1] = i-1; Type10_bonds[i][2] = i+1; 
          Type10_bonds[i][3] = i+23; Type10_bonds[i][4] = i+24; Type10_bonds[i][5] = i+45; Type10_bonds[i][6] = i+46;
        }
        // Middle layer (z-dir) interior bonds
        if(i>22 && i<33){
          Type10_bonds[i][0] = i-23; Type10_bonds[i][1] = i-22; Type10_bonds[i][2] = i-1; 
          Type10_bonds[i][3] = i+1; Type10_bonds[i][4] = i+12; Type10_bonds[i][5] = i+35; Type10_bonds[i][6] = i+36;
        }
        if(i>34 && i<45){
          Type10_bonds[i][0] = i-35; Type10_bonds[i][1] = i-34; Type10_bonds[i][2] = i-24; Type10_bonds[i][3] = i-23; 
          Type10_bonds[i][4] = i-12; Type10_bonds[i][5] = i-1; Type10_bonds[i][6] = i+1; Type10_bonds[i][7] = i+12; 
          Type10_bonds[i][8] = i+23; Type10_bonds[i][9] = i+24; Type10_bonds[i][10] = i+34; Type10_bonds[i][11] = i+35; 
        }
        if(i>46 && i<57){
          Type10_bonds[i][0] = i-36; Type10_bonds[i][1] = i-35; Type10_bonds[i][2] = i-12; Type10_bonds[i][3] = i-1; 
          Type10_bonds[i][4] = i+1; Type10_bonds[i][5] = i+22; Type10_bonds[i][6] = i+23;
        }
        // Bottom layer (z-dir) interior bonds
        if((i>58 && i<68)){ 
          Type10_bonds[i][0] = i-36; Type10_bonds[i][1] = i-35; Type10_bonds[i][2] = i-24; Type10_bonds[i][3] = i-23; 
          Type10_bonds[i][4] = i-1; Type10_bonds[i][5] = i+1; Type10_bonds[i][6] = i+11;
        } 
        if (i>69 && i<79){
          Type10_bonds[i][0] = i-35; Type10_bonds[i][1] = i-34; Type10_bonds[i][2] = i-23; Type10_bonds[i][3] = i-22; 
          Type10_bonds[i][4] = i-11; Type10_bonds[i][5] = i-1; Type10_bonds[i][6] = i+1;
        }
      }
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
  if(type == 5){
    posnorm[XDIR] = bp_data.Type5_positions[comp][0]; 
    posnorm[YDIR] = bp_data.Type5_positions[comp][1]; 
    posnorm[ZDIR] = bp_data.Type5_positions[comp][2]; 
  }
  if(type == 6){
    posnorm[XDIR] = bp_data.Type6_positions[comp][0]; 
    posnorm[YDIR] = bp_data.Type6_positions[comp][1]; 
    posnorm[ZDIR] = bp_data.Type6_positions[comp][2]; 
  }
  if(type == 7){
    posnorm[XDIR] = bp_data.Type7_positions[comp][0]; 
    posnorm[YDIR] = bp_data.Type7_positions[comp][1]; 
    posnorm[ZDIR] = bp_data.Type7_positions[comp][2]; 
  }
  if(type == 8){
    posnorm[XDIR] = bp_data.Type8_positions[comp][0]; 
    posnorm[YDIR] = bp_data.Type8_positions[comp][1]; 
    posnorm[ZDIR] = bp_data.Type8_positions[comp][2]; 
  }
  if(type == 9){
    posnorm[XDIR] = bp_data.Type9_positions[comp][0]; 
    posnorm[YDIR] = bp_data.Type9_positions[comp][1]; 
    posnorm[ZDIR] = bp_data.Type9_positions[comp][2]; 
  }
  if(type == 10){
    posnorm[XDIR] = bp_data.Type10_positions[comp][0]; 
    posnorm[YDIR] = bp_data.Type10_positions[comp][1]; 
    posnorm[ZDIR] = bp_data.Type10_positions[comp][2]; 
  }

  // Use Euler angles to calculate quaternion components
  amrex::Real eax = bp_euler_angles[XDIR];
  amrex::Real eay = bp_euler_angles[YDIR];
  amrex::Real eaz = bp_euler_angles[ZDIR];
  amrex::Real quats[4] ={cos(eaz/2.0)*cos(eay/2.0)*cos(eax/2.0) + sin(eaz/2.0)*sin(eay/2.0)*sin(eax/2.0),
                        -sin(eaz/2.0)*sin(eay/2.0)*cos(eax/2.0) + cos(eaz/2.0)*cos(eay/2.0)*sin(eax/2.0),
                         cos(eaz/2.0)*sin(eay/2.0)*cos(eax/2.0) + sin(eaz/2.0)*cos(eay/2.0)*sin(eax/2.0),
                         sin(eaz/2.0)*cos(eay/2.0)*cos(eax/2.0) - cos(eaz/2.0)*sin(eay/2.0)*sin(eax/2.0)};
  
  // Calculate rotated particle position
  rotate_vector_to_inertial(quats, posnorm, pc_pos);

  // Scale and translate particle position using radius and bp center
  pc_pos[XDIR] *= radius;
  pc_pos[YDIR] *= radius;
  pc_pos[ZDIR] *= radius;
  pc_pos[XDIR] += bp_pos[XDIR];
  pc_pos[YDIR] += bp_pos[YDIR];
  pc_pos[ZDIR] += bp_pos[ZDIR];
}

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
