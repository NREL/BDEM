struct ParticleBondData {
  // Initialize arrays that contain bond data
  int Type0_bonds[BP_NP0][MAXBONDS];    // Single particle

  // Rod-like particle
  int Type1_bonds[BP_NP1][MAXBONDS];    // 2x1
  int Type2_bonds[BP_NP2][MAXBONDS];    // 3x1
  int Type3_bonds[BP_NP3][MAXBONDS];    // 4x1
  int Type4_bonds[BP_NP4][MAXBONDS];    // 5x1
  int Type5_bonds[BP_NP5][MAXBONDS];    // 6x1
  int Type6_bonds[BP_NP6][MAXBONDS];    // 7x1
  int Type7_bonds[BP_NP7][MAXBONDS];    // 8x1
  int Type8_bonds[BP_NP8][MAXBONDS];    // 9x1
  int Type9_bonds[BP_NP9][MAXBONDS];    // 10x1

  int Type10_bonds[BP_NP10][MAXBONDS];    // 3,4,3 hex, single layer

  // AR = 4 rectangular particles
  int Type11_bonds[BP_NP11][MAXBONDS];    // 8x2, single layer
  int Type12_bonds[BP_NP12][MAXBONDS];    // 8x2, double layer
  int Type13_bonds[BP_NP13][MAXBONDS];    // 12x3, single layer
  int Type14_bonds[BP_NP14][MAXBONDS];    // 12x3, double layer
  int Type15_bonds[BP_NP15][MAXBONDS];  // 12x3, triple layer offset

  // Particle types for corn stover, approximating particle length = 8 mm
  int Type16_bonds[BP_NP16][MAXBONDS];    // 4x2, single layer
  int Type17_bonds[BP_NP17][MAXBONDS];    // 4x3, single layer
  int Type18_bonds[BP_NP18][MAXBONDS];    // 6x2, single layer
  int Type19_bonds[BP_NP19][MAXBONDS];    // 8x4, single layer
  int Type20_bonds[BP_NP20][MAXBONDS];    // 8x5, single layer
  int Type21_bonds[BP_NP21][MAXBONDS];    // 4x2, single layer, offset
  int Type22_bonds[BP_NP22][MAXBONDS];    // 5x2, single layer, offset

  // Constructor to fill in bond info
  ParticleBondData(){
      // Manually adding particle bonds 
      // TODO: Probably a better way to do this...

      // Value of -1 indicates no bond
      // TODO: find less ugly/manual method for initializing these...
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
      for(int i=0; i<BP_NP11; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type11_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP12; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type12_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP13; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type13_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP14; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type14_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP15; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type15_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP16; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type16_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP17; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type17_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP18; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type18_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP19; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type19_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP20; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type20_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP21; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type21_bonds[i][j] = -1;
          }
      }
      for(int i=0; i<BP_NP22; i++){
          for(int j=0; j<MAXBONDS; j++){
              Type22_bonds[i][j] = -1;
          }
      }

      // Particle type 1
      Type1_bonds[0][0] = 1;
      Type1_bonds[1][0] = 0;

      // Particle type 2
      Type2_bonds[0][0] = 1;
      Type2_bonds[1][0] = 0; Type2_bonds[1][1] = 2;
      Type2_bonds[2][0] = 1;

      // Particle type 3
      Type3_bonds[0][0] = 1;
      Type3_bonds[1][0] = 0; Type3_bonds[1][1] = 2;
      Type3_bonds[2][0] = 1; Type3_bonds[2][1] = 3;
      Type3_bonds[3][0] = 2;

      // Particle type 4
      Type4_bonds[0][0] = 1;
      Type4_bonds[1][0] = 0; Type4_bonds[1][1] = 2;
      Type4_bonds[2][0] = 1; Type4_bonds[2][1] = 3;
      Type4_bonds[3][0] = 2; Type4_bonds[3][1] = 4;
      Type4_bonds[4][0] = 3;

      // Particle type 9
      Type5_bonds[0][0] = 1;
      Type5_bonds[1][0] = 0; Type5_bonds[1][1] = 2;
      Type5_bonds[2][0] = 1; Type5_bonds[2][1] = 3;
      Type5_bonds[3][0] = 2; Type5_bonds[3][1] = 4;
      Type5_bonds[4][0] = 3; Type5_bonds[4][1] = 5;
      Type5_bonds[5][0] = 4;
      
      // Particle type 9
      Type6_bonds[0][0] = 1;
      Type6_bonds[1][0] = 0; Type6_bonds[1][1] = 2;
      Type6_bonds[2][0] = 1; Type6_bonds[2][1] = 3;
      Type6_bonds[3][0] = 2; Type6_bonds[3][1] = 4;
      Type6_bonds[4][0] = 3; Type6_bonds[4][1] = 5;
      Type6_bonds[5][0] = 4; Type6_bonds[5][1] = 6;
      Type6_bonds[6][0] = 5;
      
      // Particle type 9
      Type7_bonds[0][0] = 1;
      Type7_bonds[1][0] = 0; Type7_bonds[1][1] = 2;
      Type7_bonds[2][0] = 1; Type7_bonds[2][1] = 3;
      Type7_bonds[3][0] = 2; Type7_bonds[3][1] = 4;
      Type7_bonds[4][0] = 3; Type7_bonds[4][1] = 5;
      Type7_bonds[5][0] = 4; Type7_bonds[5][1] = 6;
      Type7_bonds[6][0] = 5; Type7_bonds[6][1] = 7;
      Type7_bonds[7][0] = 6;
      
      // Particle type 9
      Type8_bonds[0][0] = 1;
      Type8_bonds[1][0] = 0; Type8_bonds[1][1] = 2;
      Type8_bonds[2][0] = 1; Type8_bonds[2][1] = 3;
      Type8_bonds[3][0] = 2; Type8_bonds[3][1] = 4;
      Type8_bonds[4][0] = 3; Type8_bonds[4][1] = 5;
      Type8_bonds[5][0] = 4; Type8_bonds[5][1] = 6;
      Type8_bonds[6][0] = 5; Type8_bonds[6][1] = 7;
      Type8_bonds[7][0] = 6; Type8_bonds[7][1] = 8;
      Type8_bonds[8][0] = 7;

      // Particle type 9
      Type9_bonds[0][0] = 1;
      Type9_bonds[1][0] = 0; Type9_bonds[1][1] = 2;
      Type9_bonds[2][0] = 1; Type9_bonds[2][1] = 3;
      Type9_bonds[3][0] = 2; Type9_bonds[3][1] = 4;
      Type9_bonds[4][0] = 3; Type9_bonds[4][1] = 5;
      Type9_bonds[5][0] = 4; Type9_bonds[5][1] = 6;
      Type9_bonds[6][0] = 5; Type9_bonds[6][1] = 7;
      Type9_bonds[7][0] = 6; Type9_bonds[7][1] = 8;
      Type9_bonds[8][0] = 7; Type9_bonds[8][1] = 9;
      Type9_bonds[9][0] = 8;

      // Particle type 10
      Type10_bonds[0][0] = 1; Type10_bonds[0][1] = 3; Type10_bonds[0][2] = 4;
      Type10_bonds[1][0] = 0; Type10_bonds[1][1] = 2; Type10_bonds[1][2] = 4; Type10_bonds[1][3] = 5;
      Type10_bonds[2][0] = 1; Type10_bonds[2][1] = 5; Type10_bonds[2][2] = 6;
      Type10_bonds[3][0] = 0; Type10_bonds[3][1] = 4; Type10_bonds[3][2] = 7;
      Type10_bonds[4][0] = 0; Type10_bonds[4][1] = 1; Type10_bonds[4][2] = 3; Type10_bonds[4][3] = 5; Type10_bonds[4][4] = 7; Type10_bonds[4][5] = 8;
      Type10_bonds[5][0] = 1; Type10_bonds[5][1] = 2; Type10_bonds[5][2] = 4; Type10_bonds[5][3] = 6; Type10_bonds[5][4] = 8; Type10_bonds[5][5] = 9;
      Type10_bonds[6][0] = 2; Type10_bonds[6][1] = 5; Type10_bonds[6][2] = 9;
      Type10_bonds[7][0] = 3; Type10_bonds[7][1] = 4; Type10_bonds[7][2] = 8;
      Type10_bonds[8][0] = 4; Type10_bonds[8][1] = 5; Type10_bonds[8][2] = 7; Type10_bonds[8][3] = 9;
      Type10_bonds[9][0] = 5; Type10_bonds[9][1] = 6; Type10_bonds[9][2] = 8;

      // Particle type 11
      Type11_bonds[0][0] = 1; Type11_bonds[0][1] = 8; 
      Type11_bonds[1][0] = 0; Type11_bonds[1][1] = 2; Type11_bonds[1][2] = 9;  
      Type11_bonds[2][0] = 1; Type11_bonds[2][1] = 3; Type11_bonds[2][2] = 10;  
      Type11_bonds[3][0] = 2; Type11_bonds[3][1] = 4; Type11_bonds[3][2] = 11;
      Type11_bonds[4][0] = 3; Type11_bonds[4][1] = 5; Type11_bonds[4][2] = 12;  
      Type11_bonds[5][0] = 4; Type11_bonds[5][1] = 6; Type11_bonds[5][2] = 13;  
      Type11_bonds[6][0] = 5; Type11_bonds[6][1] = 7; Type11_bonds[6][2] = 14;  
      Type11_bonds[7][0] = 6; Type11_bonds[7][1] = 15; 
      Type11_bonds[8][0] = 0; Type11_bonds[8][1] = 9; 
      Type11_bonds[9][0] = 1; Type11_bonds[9][1] = 8; Type11_bonds[9][2] = 10;  
      Type11_bonds[10][0] = 2; Type11_bonds[10][1] = 9; Type11_bonds[10][2] = 11;  
      Type11_bonds[11][0] = 3; Type11_bonds[11][1] = 10; Type11_bonds[11][2] = 12;  
      Type11_bonds[12][0] = 4; Type11_bonds[12][1] = 11; Type11_bonds[12][2] = 13;  
      Type11_bonds[13][0] = 5; Type11_bonds[13][1] = 12; Type11_bonds[13][2] = 14;  
      Type11_bonds[14][0] = 6; Type11_bonds[14][1] = 13; Type11_bonds[14][2] = 15;  
      Type11_bonds[15][0] = 7; Type11_bonds[15][1] = 14;

      // Particle type 12
      Type12_bonds[0][0] = 1;  Type12_bonds[0][1] = 8;  Type12_bonds[0][2] = 16;
      Type12_bonds[7][0] = 6;  Type12_bonds[7][1] = 15;  Type12_bonds[7][2] = 23;
      Type12_bonds[8][0] = 0;  Type12_bonds[8][1] = 9;  Type12_bonds[8][2] = 24;
      Type12_bonds[15][0] = 7; Type12_bonds[15][1] = 14; Type12_bonds[15][2] = 31; 
      Type12_bonds[16][0] = 0; Type12_bonds[16][1] = 17; Type12_bonds[16][2] = 24;    
      Type12_bonds[23][0] = 7; Type12_bonds[23][1] = 22; Type12_bonds[23][2] = 31;
      Type12_bonds[24][0] = 8; Type12_bonds[24][1] = 16; Type12_bonds[24][2] = 25;   
      Type12_bonds[31][0] = 15; Type12_bonds[31][1] = 23; Type12_bonds[31][2] = 30;

      for(int i=0; i<BP_NP12; i++){
        if((i>0 && i<7) || (i>16 && i<23)){ // Top layer (y-dir)
          Type12_bonds[i][0] = i-1; Type12_bonds[i][1] = i+1; Type12_bonds[i][2] = i+8; Type12_bonds[i][3] = (i+16) % BP_NP12; 
        }      
        if((i>8 && i<15) || (i>24 && i<31)){ // Bottom layer (y-dir)
          Type12_bonds[i][0] = i-8; Type12_bonds[i][1] = i-1; Type12_bonds[i][2] = i+1; Type12_bonds[i][3] = (i+16) % BP_NP12; 
        }      
      } 

      // Particle type 13
      // Manually adding all corner point bonds and middle layer y-edge point bonds
      Type13_bonds[0][0] = 1;  Type13_bonds[0][1] = 12; 
      Type13_bonds[11][0] = 10; Type13_bonds[11][1] = 23;   
      Type13_bonds[12][0] = 0; Type13_bonds[12][1] = 13; Type13_bonds[12][2] = 24;
      Type13_bonds[23][0] = 11; Type13_bonds[23][1] = 22; Type13_bonds[23][2] = 35; 
      Type13_bonds[24][0] = 12; Type13_bonds[24][1] = 25; 
      Type13_bonds[35][0] = 23; Type13_bonds[35][1] = 34;

      for(int i=0; i<BP_NP13; i++){
        if(i>0 && i<11){ // Top layer (y-dir) 
          Type13_bonds[i][0] = i-1; Type13_bonds[i][1] = i+1; Type13_bonds[i][2] = i+12;  
        }
        if(i>12 && i<23){ // Middle layer (y-dir)
          Type13_bonds[i][0] = i-12; Type13_bonds[i][1] = i-1; Type13_bonds[i][2] = i+1; Type13_bonds[i][3] = i+12; 
        }
        if(i>24 && i<35){ // Bottom layer (y-dir) 
          Type13_bonds[i][0] = i-12; Type13_bonds[i][1] = i-1; Type13_bonds[i][2] = i+1; 
        }
      }

      // Particle type 14
      // Manually adding all corner point bonds and middle layer y-edge point bonds
      Type14_bonds[0][0] = 1;  Type14_bonds[0][1] = 12;  Type14_bonds[0][2] = 36;
      Type14_bonds[11][0] = 10; Type14_bonds[11][1] = 23; Type14_bonds[11][2] = 47;
      Type14_bonds[12][0] = 0; Type14_bonds[12][1] = 13; Type14_bonds[12][2] = 24; Type14_bonds[12][3] = 48;
      Type14_bonds[23][0] = 11; Type14_bonds[23][1] = 22; Type14_bonds[23][2] = 35; Type14_bonds[23][3] = 59;
      Type14_bonds[24][0] = 12; Type14_bonds[24][1] = 25; Type14_bonds[24][2] = 60;
      Type14_bonds[35][0] = 23; Type14_bonds[35][1] = 34; Type14_bonds[35][2] = 71;
      Type14_bonds[36][0] = 0; Type14_bonds[36][1] = 37; Type14_bonds[36][2] = 48;
      Type14_bonds[47][0] = 11; Type14_bonds[47][1] = 46; Type14_bonds[47][2] = 59;
      Type14_bonds[48][0] = 12; Type14_bonds[48][1] = 36; Type14_bonds[48][2] = 49; Type14_bonds[48][3] = 60;
      Type14_bonds[59][0] = 23; Type14_bonds[59][1] = 47; Type14_bonds[59][2] = 58; Type14_bonds[59][3] = 71;
      Type14_bonds[60][0] = 24; Type14_bonds[60][1] = 48; Type14_bonds[60][2] = 61;
      Type14_bonds[71][0] = 35; Type14_bonds[71][1] = 59; Type14_bonds[71][2] = 70;

      for(int i=0; i<BP_NP14; i++){
        if((i>0 && i<11) || (i>36 && i<47)){ // Top layer (y-dir) 
          Type14_bonds[i][0] = i-1; Type14_bonds[i][1] = i+1; Type14_bonds[i][2] = i+12; Type14_bonds[i][3] = (i+36) % BP_NP14; 
        }
        if((i>12 && i<23) || (i>48 && i<59)){ // Middle layer (y-dir)
          Type14_bonds[i][0] = i-12; Type14_bonds[i][1] = i-1; Type14_bonds[i][2] = i+1; Type14_bonds[i][3] = i+12; Type14_bonds[i][4] = (i+36) % BP_NP14; 
        }
        if((i>24 && i<35) || (i>60 && i<71)){ // Bottom layer (y-dir) 
          Type14_bonds[i][0] = i-12; Type14_bonds[i][1] = i-1; Type14_bonds[i][2] = i+1; Type14_bonds[i][3] = (i+36) % BP_NP14; 
        }
      }

      // Particle type 15
      // Manually adding all corner point bonds and middle layer y-edge point bonds
      Type15_bonds[0][0] = 1;   Type15_bonds[0][1] = 11;  Type15_bonds[0][2] = 22;  Type15_bonds[0][3] = 23;  Type15_bonds[0][4] = 34;  Type15_bonds[0][5] = 35;  
      Type15_bonds[10][0] = 9;  Type15_bonds[10][1] = 21; Type15_bonds[10][2] = 32; Type15_bonds[10][3] = 33; Type15_bonds[10][4] = 44; Type15_bonds[10][5] = 45; 
      Type15_bonds[11][0] = 0;  Type15_bonds[11][1] = 12; Type15_bonds[11][2] = 34; Type15_bonds[11][3] = 35; Type15_bonds[11][4] = 46; Type15_bonds[11][5] = 47; 
      Type15_bonds[21][0] = 10; Type15_bonds[21][1] = 20; Type15_bonds[21][2] = 44; Type15_bonds[21][3] = 45; Type15_bonds[21][4] = 56; Type15_bonds[21][5] = 57; 
      Type15_bonds[22][0] = 0;  Type15_bonds[22][1] = 23; Type15_bonds[22][2] = 34; Type15_bonds[22][3] = 58; 
      Type15_bonds[33][0] = 10; Type15_bonds[33][1] = 32; Type15_bonds[33][2] = 45; Type15_bonds[33][3] = 68; 
      Type15_bonds[46][0] = 11; Type15_bonds[46][1] = 34; Type15_bonds[46][2] = 47; Type15_bonds[46][3] = 69; 
      Type15_bonds[57][0] = 21; Type15_bonds[57][1] = 45; Type15_bonds[57][2] = 56; Type15_bonds[57][3] = 79; 
      Type15_bonds[58][0] = 22; Type15_bonds[58][1] = 23; Type15_bonds[58][2] = 34; Type15_bonds[58][3] = 35; Type15_bonds[58][4] = 59; Type15_bonds[58][5] = 69; 
      Type15_bonds[68][0] = 32; Type15_bonds[68][1] = 33; Type15_bonds[68][2] = 44; Type15_bonds[68][3] = 45; Type15_bonds[68][4] = 67; Type15_bonds[68][5] = 79; 
      Type15_bonds[69][0] = 34; Type15_bonds[69][1] = 35; Type15_bonds[69][2] = 46; Type15_bonds[69][3] = 47; Type15_bonds[69][4] = 58; Type15_bonds[69][5] = 70; 
      Type15_bonds[79][0] = 44; Type15_bonds[79][1] = 45; Type15_bonds[79][2] = 56; Type15_bonds[79][3] = 57; Type15_bonds[79][4] = 68; Type15_bonds[79][5] = 78; 
      Type15_bonds[34][0] = 0;  Type15_bonds[34][1] = 11; Type15_bonds[34][2] = 22; Type15_bonds[34][3] = 35; Type15_bonds[34][4] = 46; Type15_bonds[34][5] = 58; 
      Type15_bonds[34][6] = 69; 
      Type15_bonds[45][0] = 10; Type15_bonds[45][1] = 21; Type15_bonds[45][2] = 33; Type15_bonds[45][3] = 44; Type15_bonds[45][4] = 57; Type15_bonds[45][5] = 68; 
      Type15_bonds[45][6] = 79; 
      for(int i=0; i<BP_NP15; i++){
        // Top layer (z-dir) interior bonds
        if((i>0 && i<10)){ 
          Type15_bonds[i][0] = i-1; Type15_bonds[i][1] = i+1; Type15_bonds[i][2] = i+11; 
          Type15_bonds[i][3] = i+22; Type15_bonds[i][4] = i+23; Type15_bonds[i][5] = i+34; Type15_bonds[i][6] = i+35;
        } 
        if (i>11 && i<21){
          Type15_bonds[i][0] = i-11; Type15_bonds[i][1] = i-1; Type15_bonds[i][2] = i+1; 
          Type15_bonds[i][3] = i+23; Type15_bonds[i][4] = i+24; Type15_bonds[i][5] = i+45; Type15_bonds[i][6] = i+46;
        }
        // Middle layer (z-dir) interior bonds
        if(i>22 && i<33){
          Type15_bonds[i][0] = i-23; Type15_bonds[i][1] = i-22; Type15_bonds[i][2] = i-1; 
          Type15_bonds[i][3] = i+1; Type15_bonds[i][4] = i+12; Type15_bonds[i][5] = i+35; Type15_bonds[i][6] = i+36;
        }
        if(i>34 && i<45){
          Type15_bonds[i][0] = i-35; Type15_bonds[i][1] = i-34; Type15_bonds[i][2] = i-24; Type15_bonds[i][3] = i-23; 
          Type15_bonds[i][4] = i-12; Type15_bonds[i][5] = i-1; Type15_bonds[i][6] = i+1; Type15_bonds[i][7] = i+12; 
          Type15_bonds[i][8] = i+23; Type15_bonds[i][9] = i+24; Type15_bonds[i][10] = i+34; Type15_bonds[i][11] = i+35; 
        }
        if(i>46 && i<57){
          Type15_bonds[i][0] = i-36; Type15_bonds[i][1] = i-35; Type15_bonds[i][2] = i-12; Type15_bonds[i][3] = i-1; 
          Type15_bonds[i][4] = i+1; Type15_bonds[i][5] = i+22; Type15_bonds[i][6] = i+23;
        }
        // Bottom layer (z-dir) interior bonds
        if((i>58 && i<68)){ 
          Type15_bonds[i][0] = i-36; Type15_bonds[i][1] = i-35; Type15_bonds[i][2] = i-24; Type15_bonds[i][3] = i-23; 
          Type15_bonds[i][4] = i-1; Type15_bonds[i][5] = i+1; Type15_bonds[i][6] = i+11;
        } 
        if (i>69 && i<79){
          Type15_bonds[i][0] = i-35; Type15_bonds[i][1] = i-34; Type15_bonds[i][2] = i-23; Type15_bonds[i][3] = i-22; 
          Type15_bonds[i][4] = i-11; Type15_bonds[i][5] = i-1; Type15_bonds[i][6] = i+1;
        }
      }

      // Particle type 16
      Type16_bonds[0][0] = 1;  Type16_bonds[0][1] = 4;
      Type16_bonds[1][0] = 0;  Type16_bonds[1][1] = 2;  Type16_bonds[1][2] = 5;
      Type16_bonds[2][0] = 1;  Type16_bonds[2][1] = 3;  Type16_bonds[2][2] = 6;
      Type16_bonds[3][0] = 2;  Type16_bonds[3][1] = 7;
      Type16_bonds[4][0] = 0;  Type16_bonds[4][1] = 5;
      Type16_bonds[5][0] = 1;  Type16_bonds[5][1] = 4;  Type16_bonds[5][2] = 6;
      Type16_bonds[6][0] = 2;  Type16_bonds[6][1] = 5;  Type16_bonds[6][2] = 7;
      Type16_bonds[7][0] = 3;  Type16_bonds[7][1] = 6;

      // Particle type 17
      Type17_bonds[0][0] = 1;  Type17_bonds[0][1] = 4;
      Type17_bonds[1][0] = 0;  Type17_bonds[1][1] = 2;  Type17_bonds[1][2] = 5;
      Type17_bonds[2][0] = 1;  Type17_bonds[2][1] = 3;  Type17_bonds[2][2] = 6;
      Type17_bonds[3][0] = 2;  Type17_bonds[3][1] = 7;
      Type17_bonds[4][0] = 0;  Type17_bonds[4][1] = 5;  Type17_bonds[4][2] = 8;
      Type17_bonds[5][0] = 1;  Type17_bonds[5][1] = 4;  Type17_bonds[5][2] = 6; Type17_bonds[5][3] = 9;
      Type17_bonds[6][0] = 2;  Type17_bonds[6][1] = 5;  Type17_bonds[6][2] = 7; Type17_bonds[6][3] = 10;
      Type17_bonds[7][0] = 3;  Type17_bonds[7][1] = 6;  Type17_bonds[7][2] = 11;
      Type17_bonds[8][0] = 4;  Type17_bonds[8][1] = 9;
      Type17_bonds[9][0] = 5;  Type17_bonds[9][1] = 8;  Type17_bonds[9][2] = 10;
      Type17_bonds[10][0] = 6;  Type17_bonds[10][1] = 9;  Type17_bonds[10][2] = 11;
      Type17_bonds[11][0] = 7;  Type17_bonds[11][1] = 10;
      
      // Particle type 18
      Type18_bonds[0][0] = 1;   Type18_bonds[0][1] = 6;
      Type18_bonds[1][0] = 0;   Type18_bonds[1][1] = 2;   Type18_bonds[1][2] = 7;
      Type18_bonds[2][0] = 1;   Type18_bonds[2][1] = 3;   Type18_bonds[2][2] = 8;
      Type18_bonds[3][0] = 2;   Type18_bonds[3][1] = 4;   Type18_bonds[3][2] = 9;
      Type18_bonds[4][0] = 3;   Type18_bonds[4][1] = 5;   Type18_bonds[4][2] = 10;
      Type18_bonds[5][0] = 4;   Type18_bonds[5][1] = 11;
      Type18_bonds[6][0] = 0;   Type18_bonds[6][1] = 7;
      Type18_bonds[7][0] = 1;   Type18_bonds[7][1] = 6;   Type18_bonds[7][2] = 8;
      Type18_bonds[8][0] = 2;   Type18_bonds[8][1] = 7;   Type18_bonds[8][2] = 9;
      Type18_bonds[9][0] = 3;   Type18_bonds[9][1] = 8;   Type18_bonds[9][2] = 10;
      Type18_bonds[10][0] = 4;  Type18_bonds[10][1] = 9;  Type18_bonds[10][2] = 11;
      Type18_bonds[11][0] = 5;  Type18_bonds[11][1] = 10;

      // Particle type 19
      // Manually adding all corner point bonds and middle layer y-edge point bonds
      Type19_bonds[0][0] = 1;   Type19_bonds[0][1] = 8;
      Type19_bonds[7][0] = 6;  Type19_bonds[7][1] = 15;
      Type19_bonds[8][0] = 0;   Type19_bonds[8][1] = 9;  Type19_bonds[8][2] = 16;
      Type19_bonds[15][0] = 7; Type19_bonds[15][1] = 14; Type19_bonds[15][2] = 23;
      Type19_bonds[16][0] = 8; Type19_bonds[16][1] = 17; Type19_bonds[16][2] = 24;
      Type19_bonds[23][0] = 15; Type19_bonds[23][1] = 22; Type19_bonds[23][2] = 31;
      Type19_bonds[24][0] = 16;  Type19_bonds[24][1] = 25;
      Type19_bonds[31][0] = 23; Type19_bonds[31][1] = 30;

      for(int i=0; i<BP_NP19; i++){
        if(i>0 && i<7){ Type19_bonds[i][0] = i-1; Type19_bonds[i][1] = i+1; Type19_bonds[i][2] = i+8; }
        if((i>8 && i<15) || (i>16 && i<23)){ Type19_bonds[i][0] = i-8; Type19_bonds[i][1] = i-1; Type19_bonds[i][2] = i+1; Type19_bonds[i][3] = i+8;}
        if(i>24 && i<31){ Type19_bonds[i][0] = i-8; Type19_bonds[i][1] = i-1; Type19_bonds[i][2] = i+1; }
      }

      // Particle type 20
      // Manually adding all corner point bonds and middle layer y-edge point bonds
      Type20_bonds[0][0] = 1;   Type20_bonds[0][1] = 8;
      Type20_bonds[7][0] = 6;   Type20_bonds[7][1] = 15;
      Type20_bonds[8][0] = 0;   Type20_bonds[8][1] = 9;   Type20_bonds[8][2] = 16;
      Type20_bonds[15][0] = 7;  Type20_bonds[15][1] = 14; Type20_bonds[15][2] = 23;
      Type20_bonds[16][0] = 8;  Type20_bonds[16][1] = 17; Type20_bonds[16][2] = 24;
      Type20_bonds[23][0] = 15; Type20_bonds[23][1] = 22; Type20_bonds[23][2] = 31;
      Type20_bonds[24][0] = 16; Type20_bonds[24][1] = 25; Type20_bonds[24][2] = 32;
      Type20_bonds[31][0] = 23; Type20_bonds[31][1] = 30; Type20_bonds[31][2] = 39;
      Type20_bonds[32][0] = 24;   Type20_bonds[32][1] = 33;
      Type20_bonds[39][0] = 31;   Type20_bonds[39][1] = 38;

      for(int i=0; i<BP_NP20; i++){
        if(i>0 && i<7){ Type20_bonds[i][0] = i-1; Type20_bonds[i][1] = i+1; Type20_bonds[i][2] = i+8; }
        if((i>8 && i<15) || (i>16 && i<23) || (i>24 && i<31)){ Type20_bonds[i][0] = i-8; Type20_bonds[i][1] = i-1; Type20_bonds[i][2] = i+1; Type20_bonds[i][3] = i+8;}
        if(i>32 && i<39){ Type20_bonds[i][0] = i-8; Type20_bonds[i][1] = i-1; Type20_bonds[i][2] = i+1; }
      }

      // Particle type 21
      Type21_bonds[0][0] = 1;  Type21_bonds[0][1] = 4;
      Type21_bonds[1][0] = 0;  Type21_bonds[1][1] = 2;  Type21_bonds[1][2] = 5;
      Type21_bonds[2][0] = 1;  Type21_bonds[2][1] = 3;  Type21_bonds[2][2] = 6;
      Type21_bonds[3][0] = 2;  Type21_bonds[3][1] = 7;
      Type21_bonds[4][0] = 0;  Type21_bonds[4][1] = 5;
      Type21_bonds[5][0] = 1;  Type21_bonds[5][1] = 4;  Type21_bonds[5][2] = 6;
      Type21_bonds[6][0] = 2;  Type21_bonds[6][1] = 5;  Type21_bonds[6][2] = 7;
      Type21_bonds[7][0] = 3;  Type21_bonds[7][1] = 6;
      
      // Particle type 22
      Type22_bonds[0][0] = 1;  Type22_bonds[0][1] = 5;
      Type22_bonds[1][0] = 0;  Type22_bonds[1][1] = 2;  Type22_bonds[1][2] = 6;
      Type22_bonds[2][0] = 1;  Type22_bonds[2][1] = 3;  Type22_bonds[2][2] = 7;
      Type22_bonds[3][0] = 2;  Type22_bonds[3][1] = 4;  Type22_bonds[3][2] = 8;
      Type22_bonds[4][0] = 3;  Type22_bonds[4][1] = 9;
      Type22_bonds[5][0] = 0;  Type22_bonds[5][1] = 6;
      Type22_bonds[6][0] = 1;  Type22_bonds[6][1] = 5;  Type22_bonds[6][2] = 7;
      Type22_bonds[7][0] = 2;  Type22_bonds[7][1] = 6;  Type22_bonds[7][2] = 8;
      Type22_bonds[8][0] = 3;  Type22_bonds[8][1] = 7;  Type22_bonds[8][2] = 9;
      Type22_bonds[9][0] = 4;  Type22_bonds[9][1] = 8;
  }
};
