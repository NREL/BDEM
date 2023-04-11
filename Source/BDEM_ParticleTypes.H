struct ParticleBondData {
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
  ParticleBondData(){
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
