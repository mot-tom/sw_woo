/******************************************************************************
                                watanabe.c

                        --- WATANABE potential ---

         Copyright (c) Takanobu Watanabe 2002  All Rights Reserved

         This software is freely distributed without licensing fees
         and is provided without warrantee of any kind.
*******************************************************************************

 Breaf history

  2002 May      New 3 body term depending on coordination number was invented.
  2002 Aug. 22  Formula of force was found.
  2002 Aug. 24  Hamiltonian was preserved under temperature control.
  2002 Aug. 27  Internal stress tensor calculation code was completed.
                Hamiltonian was preserved under pressure control.
  2002 Aug. 28  Program was refined.

  2002 Sep.  2  Three-body function was modified to the double gamma type.

*******************************************************************************/
#include "mddef.h"

namespace WASEDA_MD_LABO {

/*
#define DEBUG
*/
static void
  count_coordinate(void),
  two_body(void),
  three_body(void),
  force_Si_Si(int i, int j, double dxij, double dyij, double dzij),
  force_Si_O (int i, int j, double dxij, double dyij, double dzij),
  force_O_Si (int i, int j, double dxij, double dyij, double dzij),
  force_O_O  (int i, int j, double dxij, double dyij, double dzij),
  force_simple_simple(int i, int j, double dxij, double dyij, double dzij),
  two_body2(void),
  two_body3(void),
  three_body(void),
  force_Si_O2 (int i, int j, double dxij, double dyij, double dzij),
  force_O_Si2 (int i, int j, double dxij, double dyij, double dzij),
  force_Si_other(int i),
  force_O_other(int i),
  force_Si_Si3 (int i, int j, double dxij, double dyij, double dzij);

static double
  fc(double r, WatanabeCutoff cutoff_param),
  d_fc(double r, WatanabeCutoff cutoff_param),
  gSi(double co),
  d_gSi(double co),
  gO(double co),
  d_gO(double co),
  lambda(double co, Watanabe3Body * parameter),
  d_lambda(double co, Watanabe3Body * parameter),
  Fc(double r, double R),
  d_Fc_r(double r, double R),
  d_Fc_R(double r, double R),
  Rc(double co),
  d_Rc(double co),
  bsSi(double co),
  d_bsSi(double co);

void set_WATANABE1_parameter(void)
{
  int i,j,k;

  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      wb2[i][j].A = 0.0;
      wb2[i][j].B = 0.0;
      wb2[i][j].p = 0.0;
      wb2[i][j].q = 0.0;
      wb2[i][j].a = 1.0;
    }
  }
  wb2[Si][Si].A = 7.04955;
  wb2[Si][Si].B = 0.60222;
  wb2[Si][Si].p = 4.0;
  wb2[Si][Si].q = 0.0;
  wb2[Si][Si].a = 1.80;

  wb2[Si][O].A = 21.0;
  wb2[Si][O].B = 0.038;
  wb2[Si][O].p = 5.3;
  wb2[Si][O].q = -1.1;
  wb2[Si][O].a = 1.30;

  wb2[O][O].A = -12.29;
  wb2[O][O].B = 0.00;
  wb2[O][O].p = 0.00;
  wb2[O][O].q = 2.24;
  wb2[O][O].a = 1.25;

  wb2[O][Si].A = wb2[Si][O].A;
  wb2[O][Si].B = wb2[Si][O].B;
  wb2[O][Si].p = wb2[Si][O].p;
  wb2[O][Si].q = wb2[Si][O].q;
  wb2[O][Si].a = wb2[Si][O].a;

  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      for(k=0;k<3;k++){
        wb3[i][j][k].lam    = 0.0;
        wb3[i][j][k].gam[0] = 0.0;
        wb3[i][j][k].gam[1] = 0.0;
        wb3[i][j][k].a[0]   = 0.1;
        wb3[i][j][k].a[1]   = 0.1;
        wb3[i][j][k].mu     = 0.0;
        wb3[i][j][k].nu[0]  = 0.0;
        wb3[i][j][k].nu[1]  = 0.0;
        wb3[i][j][k].b[0]   = 0.1;
        wb3[i][j][k].b[1]   = 0.1;
        wb3[i][j][k].cos0   = 0.0;
        wb3[i][j][k].alpha  = 0.0;
        wb3[i][j][k].beta   = 0.0;
        wb3[i][j][k].g[0]   = 0.0;
        wb3[i][j][k].g[1]   = 0.0;
        wb3[i][j][k].g[2]   = 0.0;
      }
    }
  }

  wb3[Si][Si][Si].lam    = 5.878;
  wb3[Si][Si][Si].gam[0] = 1.01;
  wb3[Si][Si][Si].gam[1] = 1.01;
  wb3[Si][Si][Si].a[0]   = 1.80;
  wb3[Si][Si][Si].a[1]   = 1.80;
  wb3[Si][Si][Si].mu     = 4.0;
  wb3[Si][Si][Si].nu[0]  = 0.088;
  wb3[Si][Si][Si].nu[1]  = 0.088;
  wb3[Si][Si][Si].b[0]   = 1.2;
  wb3[Si][Si][Si].b[1]   = 1.2;
  wb3[Si][Si][Si].cos0   = -0.3333333;
  wb3[Si][Si][Si].alpha  = 0.3;
  wb3[Si][Si][Si].beta   = -1.35;
  wb3[Si][Si][Si].g[0]   = 1.6;
  wb3[Si][Si][Si].g[1]   = 6.0;
  wb3[Si][Si][Si].g[2]   = 0.0;
/*  local01 version 
  wb3[Si][Si][Si].lam    = 5.878;
  wb3[Si][Si][Si].gam[0] = 1.01;
  wb3[Si][Si][Si].gam[1] = 1.01;
  wb3[Si][Si][Si].a[0]   = 1.80;
  wb3[Si][Si][Si].a[1]   = 1.80;
  wb3[Si][Si][Si].mu     = 4.0;
  wb3[Si][Si][Si].nu[0]  = 0.2;
  wb3[Si][Si][Si].nu[1]  = 0.2;
  wb3[Si][Si][Si].b[0]   = 1.3;
  wb3[Si][Si][Si].b[1]   = 1.3;
  wb3[Si][Si][Si].cos0   = -0.3333333;
  wb3[Si][Si][Si].alpha  = 0.3;
  wb3[Si][Si][Si].beta   = -1.35;
  wb3[Si][Si][Si].g[0]   = 1.5;
  wb3[Si][Si][Si].g[1]   = 6.0;
  wb3[Si][Si][Si].g[2]   = 0.0;
*/

  wb3[Si][Si][O].lam    = 3.0;
  wb3[Si][Si][O].gam[0] = 0.032;
  wb3[Si][Si][O].gam[1] = 0.124;
  wb3[Si][Si][O].a[0]   = 1.15;
  wb3[Si][Si][O].a[1]   = 1.10;
  wb3[Si][Si][O].mu     = 5.878;
  wb3[Si][Si][O].nu[0]  = 1.01;
  wb3[Si][Si][O].nu[1]  = 1.082079;
  wb3[Si][Si][O].b[0]   = 1.8;
  wb3[Si][Si][O].b[1]   = 1.5;
  wb3[Si][Si][O].cos0   = -0.3333333;
  wb3[Si][Si][O].alpha  = 0.3;
  wb3[Si][Si][O].beta   = 0.3;
  wb3[Si][Si][O].g[0]   = 3.6;
  wb3[Si][Si][O].g[1]   = 2.0;
  wb3[Si][Si][O].g[2]   = 2.6;
/* local01 version 
  wb3[Si][Si][O].lam    = 3.0;
  wb3[Si][Si][O].gam[0] = 0.2;
  wb3[Si][Si][O].gam[1] = 0.2;
  wb3[Si][Si][O].a[0]   = 1.30;
  wb3[Si][Si][O].a[1]   = 1.30;
  wb3[Si][Si][O].mu     = 5.878;
  wb3[Si][Si][O].nu[0]  = 1.01;
  wb3[Si][Si][O].nu[1]  = 0.7843;
  wb3[Si][Si][O].b[0]   = 1.8;
  wb3[Si][Si][O].b[1]   = 1.3;
  wb3[Si][Si][O].cos0   = -0.3333333;
  wb3[Si][Si][O].alpha  = 0.3;
  wb3[Si][Si][O].beta   = 0.3;
  wb3[Si][Si][O].g[0]   = 2.33;
  wb3[Si][Si][O].g[1]   = 1.9;
  wb3[Si][Si][O].g[2]   = 3.0;
*/

  wb3[Si][O][Si].lam    = 2.50;
  wb3[Si][O][Si].gam[0] = 0.15;
  wb3[Si][O][Si].gam[1] = 0.15;
  wb3[Si][O][Si].a[0]   = 1.2;
  wb3[Si][O][Si].a[1]   = 1.2;
  wb3[Si][O][Si].mu     = 0.0;
  wb3[Si][O][Si].nu[0]  = 0.0;
  wb3[Si][O][Si].nu[1]  = 0.0;
  wb3[Si][O][Si].b[0]   = 0.1;
  wb3[Si][O][Si].b[1]   = 0.1;
  wb3[Si][O][Si].cos0   = -0.812;
  wb3[Si][O][Si].alpha  = 1.6;
  wb3[Si][O][Si].beta   = 1.6;
  wb3[Si][O][Si].g[0]   = 0.0;
  wb3[Si][O][Si].g[1]   = 0.0;
  wb3[Si][O][Si].g[2]   = 0.0;
/* local01 version 
  wb3[Si][O][Si].lam    = 1.40;
  wb3[Si][O][Si].gam[0] = 0.34;
  wb3[Si][O][Si].gam[1] = 0.34;
  wb3[Si][O][Si].a[0]   = 1.20;
  wb3[Si][O][Si].a[1]   = 1.20;
  wb3[Si][O][Si].mu     = 0.0;
  wb3[Si][O][Si].nu[0]  = 0.0;
  wb3[Si][O][Si].nu[1]  = 0.0;
  wb3[Si][O][Si].b[0]   = 0.1;
  wb3[Si][O][Si].b[1]   = 0.1;
  wb3[Si][O][Si].cos0   = -0.812;
  wb3[Si][O][Si].alpha  = 2.2;
  wb3[Si][O][Si].beta   = 2.2;
  wb3[Si][O][Si].g[0]   = 0.0;
  wb3[Si][O][Si].g[1]   = 0.0;
  wb3[Si][O][Si].g[2]   = 0.0;
*/

  wb3[O][Si][Si].lam    = wb3[Si][Si][O].lam;
  wb3[O][Si][Si].gam[0] = wb3[Si][Si][O].gam[1];
  wb3[O][Si][Si].gam[1] = wb3[Si][Si][O].gam[0];
  wb3[O][Si][Si].a[0]   = wb3[Si][Si][O].a[1];
  wb3[O][Si][Si].a[1]   = wb3[Si][Si][O].a[0];
  wb3[O][Si][Si].mu     = wb3[Si][Si][O].mu;
  wb3[O][Si][Si].nu[0]  = wb3[Si][Si][O].nu[1];
  wb3[O][Si][Si].nu[1]  = wb3[Si][Si][O].nu[0];
  wb3[O][Si][Si].b[0]   = wb3[Si][Si][O].b[1];
  wb3[O][Si][Si].b[1]   = wb3[Si][Si][O].b[0];
  wb3[O][Si][Si].cos0   = wb3[Si][Si][O].cos0;
  wb3[O][Si][Si].alpha  = wb3[Si][Si][O].alpha;
  wb3[O][Si][Si].beta   = wb3[Si][Si][O].beta;
  wb3[O][Si][Si].g[0]   = wb3[Si][Si][O].g[0];
  wb3[O][Si][Si].g[1]   = wb3[Si][Si][O].g[1];
  wb3[O][Si][Si].g[2]   = wb3[Si][Si][O].g[2];

  wb3[O][Si][O].lam    = 10.5; 
  wb3[O][Si][O].gam[0] = 0.310;
  wb3[O][Si][O].gam[1] = 0.310;
  wb3[O][Si][O].a[0]   = 1.10;
  wb3[O][Si][O].a[1]   = 1.10;
  wb3[O][Si][O].mu     = 0.0;
  wb3[O][Si][O].nu[0]  = 0.0;
  wb3[O][Si][O].nu[1]  = 0.0;
  wb3[O][Si][O].b[0]   = 0.5;
  wb3[O][Si][O].b[1]   = 0.5;
  wb3[O][Si][O].cos0   = -0.333333333;
  wb3[O][Si][O].alpha  = 0.0;
  wb3[O][Si][O].beta   = 0.0;
  wb3[O][Si][O].g[0]   = 0.0;
  wb3[O][Si][O].g[1]   = 0.0;
  wb3[O][Si][O].g[2]   = 0.0;

/* local01 version 
  wb3[O][Si][O].lam    = 10.5;
  wb3[O][Si][O].gam[0] = 0.50;
  wb3[O][Si][O].gam[1] = 0.50;
  wb3[O][Si][O].a[0]   = 1.30;
  wb3[O][Si][O].a[1]   = 1.30;
  wb3[O][Si][O].mu     = 0.0;
  wb3[O][Si][O].nu[0]  = 0.0;
  wb3[O][Si][O].nu[1]  = 0.0;
  wb3[O][Si][O].b[0]   = 0.5;
  wb3[O][Si][O].b[1]   = 0.5;
  wb3[O][Si][O].cos0   = -0.333333333;
  wb3[O][Si][O].alpha  = 0.0;
  wb3[O][Si][O].beta   = 0.0;
  wb3[O][Si][O].g[0]   = 0.0;
  wb3[O][Si][O].g[1]   = 0.0;
  wb3[O][Si][O].g[2]   = 0.0;
*/

  if(charac[Omol].number > 0){
    // Ar is assumed as O in O2 molecule here
    // added 2006.3.28 for nanoscale oxidation process simulator
    wb2[Ar][Ar].A = 11.57;
    wb2[Ar][Ar].B = 0.4598;
    wb2[Ar][Ar].p = 3.6;
    wb2[Ar][Ar].q = 2.63;
    wb2[Ar][Ar].a = 1.25;

    wb2[Ar][Si].A = -46.0 * eV2eps * pow(1.0/sig, 3.6);
    wb2[Ar][Si].B = 0.0;
    wb2[Ar][Si].p = 0.0;
    wb2[Ar][Si].q = 3.6;
    wb2[Ar][Si].a = 3.8/sig;

    wb2[Ar][O].A = -38.0 * eV2eps * pow(1.0/sig, 3.0);
    wb2[Ar][O].B = 0.0;
    wb2[Ar][O].p = 0.0;
    wb2[Ar][O].q = 3.0;
    wb2[Ar][O].a = 3.0/sig;
  }

#ifndef SILENT
  printf("###### Constants in the potential function ######\n");
  printf("<< Tow-Body Term >>\n");
  printf(" i  j        Aij            Bij           pij        qij        aij\n");
  printf("----------------------------------------------------------------------\n");
  for(i = 0; i < 2; i++){
    for(j = i; j < 2; j++){
      printf("%s %s | %14.10f %14.10f %10.6f %10.6f %10.6f\n",
             charac[i].name, charac[j].name,
             wb2[i][j].A,wb2[i][j].B,wb2[i][j].p,
             wb2[i][j].q,wb2[i][j].a);
    }
  }
  printf("----------------------------------------------------------------------\n\n");
  printf("<< Three-Body Term >>\n");
  printf(" j  i  k      lam        gam^ij     gam^ik     a^ij       a^ik     cos(theta^0)\n");
  printf("-------------------------------------------------------------------------------\n");
  for(i = 0; i < 2; i++){
    for(j = 0; j < 2; j++){
      printf("%s %s %s | %10.6f %10.6f %10.6f %10.6f %10.6f %13.7f\n",
             charac[i].name, charac[i].name,
             charac[j].name,
             wb3[i][i][j].lam,wb3[i][i][j].gam[0],wb3[i][i][j].gam[1],
             wb3[i][i][j].a[0],wb3[i][i][j].a[1],wb3[i][i][j].cos0);
      if(i != j){
        printf("%s %s %s | %10.6f %10.6f %10.6f %10.6f %10.6f %13.7f\n",
               charac[i].name, charac[j].name,
               charac[i].name,
               wb3[i][j][i].lam,wb3[i][j][i].gam[0],wb3[i][j][i].gam[1],
               wb3[i][j][i].a[0],wb3[i][j][i].a[1],wb3[i][j][i].cos0);
      }
    }
  }
  printf("-------------------------------------------------------------------------------\n");
  printf(" j  i  k      mu         nu^ij      nu^ik      b^ij       b^ik  \n");
  printf("----------------------------------------------------------------\n");
  for(i = 0; i < 2; i++){
    for(j = 0; j < 2; j++){
      printf("%s %s %s | %10.6f %10.6f %10.6f %10.6f %10.6f\n",
             charac[i].name, charac[i].name,
             charac[j].name,
             wb3[i][i][j].mu,wb3[i][i][j].nu[0],wb3[i][i][j].nu[1],
             wb3[i][i][j].b[0],wb3[i][i][j].b[1]);
      if(i != j){
        printf("%s %s %s | %10.6f %10.6f %10.6f %10.6f %10.6f\n",
               charac[i].name, charac[j].name,
               charac[i].name,
               wb3[i][j][i].mu,wb3[i][j][i].nu[0],wb3[i][j][i].nu[1],
               wb3[i][j][i].b[0],wb3[i][j][i].b[1]);
      }
    }
  }

  printf("----------------------------------------------------------------\n");
  printf(" j  i  k      alpha      g0         g1         g2\n");
  printf("------------------------------------------------------\n");
  for(i = 0; i < 2; i++){
    for(j = 0; j < 2; j++){
      printf("%s %s %s | %10.6f %10.6f %10.6f %10.6f\n",
             charac[i].name, charac[i].name,
             charac[j].name,
             wb3[i][i][j].alpha,wb3[i][i][j].g[0],wb3[i][i][j].g[1],
             wb3[i][i][j].g[2]);
      if(i != j){
        printf("%s %s %s | %10.6f %10.6f %10.6f %10.6f\n",
               charac[i].name, charac[j].name,
               charac[i].name,
               wb3[i][j][i].alpha,wb3[i][j][i].g[0],wb3[i][j][i].g[1],
               wb3[i][j][i].g[2]);
      }
    }
  }
  printf("------------------------------------------------------\n\n");
#endif

  wbcutsio.R   = 1.20;     
  wbcutsio.D   = 0.05;    
  wbcutsio.RpD = wbcutsio.R + wbcutsio.D;
/*
  wbcutsio.R   = 0.92;
  wbcutsio.D   = 0.05;
  wbcutsio.RpD = wbcutsio.R + wbcutsio.D;

  wbcutsio.R   = 1.04;
  wbcutsio.D   = 0.26;
  wbcutsio.RpD = wbcutsio.R + wbcutsio.D;
*/
  wbcutsisi.R   = 1.60;
  wbcutsisi.D   = 0.05;
  wbcutsisi.RpD = wbcutsisi.R + wbcutsisi.D;

  wbsat.na =  0.712;
  wbsat.nb =  0.522;
  wbsat.nc =  5.274;
  wbsat.nd = -0.0372;
  wbsat.ne =  4.52;

  wbsat.sia = 0.303672895420570121;
  wbsat.sib = 3.93232680731067283;
  wbsat.sic = 0.253450103965147155;
  wbsat.sid = 0.856344826609515408;
  wbsat.ha  = 0.0;

#ifndef SILENT
  printf("<< Si-O Bond-Softening Function >>\n");
  printf("    m1         m2         m3         m4         m5\n");
  printf("-------------------------------------------------------\n");
  printf("%10.6f %10.6f %10.6f %10.6f %10.6f\n",
         wbsat.nc, wbsat.na, wbsat.nb, wbsat.nd, wbsat.ne);
  printf("-------------------------------------------------------\n\n");
  printf("<< Si-O Cutoff Function >>\n");
  printf("     R          D\n");
  printf("------------------------\n");
  printf("%10.6f %10.6f\n", wbcutsio.R, wbcutsio.D);
  printf("------------------------\n\n");
#endif

  double cut, a;

  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      cut = wb2[i][j].a;
      for(k = 0; k < 3; k++) {
        a = wb3[i][j][k].a[0];
        if(cut < a) cut = a;
        a = wb3[j][i][k].a[0];
        if(cut < a) cut = a;
        a = wb3[i][j][k].b[0];
        if(cut < a) cut = a;
        a = wb3[j][i][k].b[0];
        if(cut < a) cut = a;
      }
      cutoff[i][j] = cut;
    }
  }
}

void set_WATANABE2_parameter(void)
{
  int i,j,k;

  for(i=0;i<2;i++){
    for(j=0;j<2;j++){
      wb2[i][j].A = 0.0;
      wb2[i][j].B = 0.0;
      wb2[i][j].p = 0.0;
      wb2[i][j].q = 0.0;
      wb2[i][j].a = 1.0;
    }
  }

/* 222222222222222222222  */

  wb2[Si][Si].A = 7.04955;
  wb2[Si][Si].B = 0.60222;
  wb2[Si][Si].p = 4.0;
  wb2[Si][Si].q = 0.0;
  wb2[Si][Si].a = 1.80;

  wb2[Si][O].A = 21.0;
  wb2[Si][O].B = 0.038;
  wb2[Si][O].p = 5.3;
  wb2[Si][O].q = -1.1;
  wb2[Si][O].a = 1.30;

  wb2[O][O].A = -12.29;
  wb2[O][O].B = 0.00;
  wb2[O][O].p = 0.00;
  wb2[O][O].q = 2.24;
  wb2[O][O].a = 1.25;

  wb2[O][Si].A = wb2[Si][O].A;
  wb2[O][Si].B = wb2[Si][O].B;
  wb2[O][Si].p = wb2[Si][O].p;
  wb2[O][Si].q = wb2[Si][O].q;
  wb2[O][Si].a = wb2[Si][O].a;

  for(i=0;i<2;i++){
    for(j=0;j<2;j++){
      for(k=0;k<2;k++){
        wb3[i][j][k].lam    = 0.0;
        wb3[i][j][k].gam[0] = 0.0;
        wb3[i][j][k].gam[1] = 0.0;
        wb3[i][j][k].a[0]   = 0.1;
        wb3[i][j][k].a[1]   = 0.1;
        wb3[i][j][k].mu     = 0.0;
        wb3[i][j][k].nu[0]  = 0.0;
        wb3[i][j][k].nu[1]  = 0.0;
        wb3[i][j][k].b[0]   = 0.1;
        wb3[i][j][k].b[1]   = 0.1;
        wb3[i][j][k].cos0   = 0.0;
        wb3[i][j][k].alpha  = 0.0;
        wb3[i][j][k].beta   = 0.0;
        wb3[i][j][k].g[0]   = 0.0;
        wb3[i][j][k].g[1]   = 0.0;
        wb3[i][j][k].g[2]   = 0.0;
      }
    }
  }

  wb3[Si][Si][Si].lam    = 5.878;
  wb3[Si][Si][Si].gam[0] = 1.01;
  wb3[Si][Si][Si].gam[1] = 1.01;
  wb3[Si][Si][Si].a[0]   = 1.80;
  wb3[Si][Si][Si].a[1]   = 1.80;
  wb3[Si][Si][Si].mu     = 4.0;
  wb3[Si][Si][Si].nu[0]  = 0.2;
  wb3[Si][Si][Si].nu[1]  = 0.2;
  wb3[Si][Si][Si].b[0]   = 1.3;
  wb3[Si][Si][Si].b[1]   = 1.3;
  wb3[Si][Si][Si].cos0   = -0.3333333;
  wb3[Si][Si][Si].alpha  = 0.3;
  wb3[Si][Si][Si].beta   = -1.35;
  wb3[Si][Si][Si].g[0]   = 1.5;
  wb3[Si][Si][Si].g[1]   = 6.0;
  wb3[Si][Si][Si].g[2]   = 0.0;

  wb3[Si][Si][O].lam    = 3.0;
  wb3[Si][Si][O].gam[0] = 0.2;
  wb3[Si][Si][O].gam[1] = 0.2;
  wb3[Si][Si][O].a[0]   = 1.3;
  wb3[Si][Si][O].a[1]   = 1.3;
  wb3[Si][Si][O].mu     = 5.878;
  wb3[Si][Si][O].nu[0]  = 1.01;
  wb3[Si][Si][O].nu[1]  = 0.7843;
  wb3[Si][Si][O].b[0]   = 1.8;
  wb3[Si][Si][O].b[1]   = 1.3;
  wb3[Si][Si][O].cos0   = -0.3333333;
  wb3[Si][Si][O].alpha  = 0.3;
  wb3[Si][Si][O].beta   = 0.3;
  wb3[Si][Si][O].g[0]   = 2.33;
  wb3[Si][Si][O].g[1]   = 1.9;
  wb3[Si][Si][O].g[2]   = 3.0;

  wb3[Si][O][Si].lam    = 2.5;
  wb3[Si][O][Si].gam[0] = 0.15;
  wb3[Si][O][Si].gam[1] = 0.15;
  wb3[Si][O][Si].a[0]   = 1.20;
  wb3[Si][O][Si].a[1]   = 1.20;
  wb3[Si][O][Si].mu     = 0.0;
  wb3[Si][O][Si].nu[0]  = 0.0;
  wb3[Si][O][Si].nu[1]  = 0.0;
  wb3[Si][O][Si].b[0]   = 0.1;
  wb3[Si][O][Si].b[1]   = 0.1;
  wb3[Si][O][Si].cos0   = -0.812;
  wb3[Si][O][Si].alpha  = 1.6;
  wb3[Si][O][Si].beta   = 1.6;
  wb3[Si][O][Si].g[0]   = 0.0;
  wb3[Si][O][Si].g[1]   = 0.0;
  wb3[Si][O][Si].g[2]   = 0.0;

  wb3[O][Si][Si].lam    = wb3[Si][Si][O].lam;
  wb3[O][Si][Si].gam[0] = wb3[Si][Si][O].gam[1];
  wb3[O][Si][Si].gam[1] = wb3[Si][Si][O].gam[0];
  wb3[O][Si][Si].a[0]   = wb3[Si][Si][O].a[1];
  wb3[O][Si][Si].a[1]   = wb3[Si][Si][O].a[0];
  wb3[O][Si][Si].mu     = wb3[Si][Si][O].mu;
  wb3[O][Si][Si].nu[0]  = wb3[Si][Si][O].nu[1];
  wb3[O][Si][Si].nu[1]  = wb3[Si][Si][O].nu[0];
  wb3[O][Si][Si].b[0]   = wb3[Si][Si][O].b[1];
  wb3[O][Si][Si].b[1]   = wb3[Si][Si][O].b[0];
  wb3[O][Si][Si].cos0   = wb3[Si][Si][O].cos0;
  wb3[O][Si][Si].alpha  = wb3[Si][Si][O].alpha;
  wb3[O][Si][Si].beta   = wb3[Si][Si][O].beta;
  wb3[O][Si][Si].g[0]   = wb3[Si][Si][O].g[0];
  wb3[O][Si][Si].g[1]   = wb3[Si][Si][O].g[1];
  wb3[O][Si][Si].g[2]   = wb3[Si][Si][O].g[2];

  wb3[O][Si][O].lam    = 10.5;
  wb3[O][Si][O].gam[0] = 0.5;
  wb3[O][Si][O].gam[1] = 0.5;
  wb3[O][Si][O].a[0]   = 1.30;
  wb3[O][Si][O].a[1]   = 1.30;
  wb3[O][Si][O].mu     = 0.0;
  wb3[O][Si][O].nu[0]  = 0.0;
  wb3[O][Si][O].nu[1]  = 0.0;
  wb3[O][Si][O].b[0]   = 0.5;
  wb3[O][Si][O].b[1]   = 0.5;
  wb3[O][Si][O].cos0   = -0.333333333;
  wb3[O][Si][O].alpha  = 0.0;
  wb3[O][Si][O].beta   = 0.0;
  wb3[O][Si][O].g[0]   = 0.0;
  wb3[O][Si][O].g[1]   = 0.0;
  wb3[O][Si][O].g[2]   = 0.0;

#ifndef SILENT
  printf("###### Constants in the potential function ######\n");
  printf("<< Tow-Body Term >>\n");
  printf(" i  j        Aij            Bij           pij        qij        aij\n");
  printf("----------------------------------------------------------------------\n");
  for(i = 0; i < 2; i++){
    for(j = i; j < 2; j++){
      printf("%s %s | %14.10f %14.10f %10.6f %10.6f %10.6f\n",
             charac[i].name, charac[j].name,
             wb2[i][j].A,wb2[i][j].B,wb2[i][j].p,
             wb2[i][j].q,wb2[i][j].a);
    }
  }
  printf("----------------------------------------------------------------------\n\n");
  printf("<< Three-Body Term >>\n");
  printf(" j  i  k      lam        gam^ij     gam^ik     a^ij       a^ik     cos(theta^0)\n");
  printf("-------------------------------------------------------------------------------\n");
  for(i = 0; i < 2; i++){
    for(j = 0; j < 2; j++){
      printf("%s %s %s | %10.6f %10.6f %10.6f %10.6f %10.6f %13.7f\n",
             charac[i].name, charac[i].name,
             charac[j].name,
             wb3[i][i][j].lam,wb3[i][i][j].gam[0],wb3[i][i][j].gam[1],
             wb3[i][i][j].a[0],wb3[i][i][j].a[1],wb3[i][i][j].cos0);
      if(i != j){
        printf("%s %s %s | %10.6f %10.6f %10.6f %10.6f %10.6f %13.7f\n",
               charac[i].name, charac[j].name,
               charac[i].name,
               wb3[i][j][i].lam,wb3[i][j][i].gam[0],wb3[i][j][i].gam[1],
               wb3[i][j][i].a[0],wb3[i][j][i].a[1],wb3[i][j][i].cos0);
      }
    }
  }
  printf("-------------------------------------------------------------------------------\n");
  printf(" j  i  k      alpha      g0         g1         g2\n");
  printf("------------------------------------------------------\n");
  for(i = 0; i < 2; i++){
    for(j = 0; j < 2; j++){
      printf("%s %s %s | %10.6f %10.6f %10.6f %10.6f\n",
             charac[i].name, charac[i].name,
             charac[j].name,
             wb3[i][i][j].alpha,wb3[i][i][j].g[0],wb3[i][i][j].g[1],
             wb3[i][i][j].g[2]);
      if(i != j){
        printf("%s %s %s | %10.6f %10.6f %10.6f %10.6f\n",
               charac[i].name, charac[j].name,
               charac[i].name,
               wb3[i][j][i].alpha,wb3[i][j][i].g[0],wb3[i][j][i].g[1],
               wb3[i][j][i].g[2]);
      }
    }
  }
  printf("------------------------------------------------------\n\n");
#endif

  wbcutsio.R   = 1.2;   /* [sig] */
  wbcutsio.D   = 0.08;   /* [sig] */
  wbcutsio.RpD = wbcutsio.R + wbcutsio.D;

  wbcutsisi.R   = 1.6;   /* [sig] */
  wbcutsisi.D   = 0.05;   /* [sig] */
  wbcutsisi.RpD = wbcutsisi.R + wbcutsisi.D;

  wbsat.na = 0.712;
  wbsat.nb = 0.522;
  wbsat.nc = 5.274;
  wbsat.nd = -0.0372;
  wbsat.ne = 4.52;

  wbsat.sia = 0.303672895420570121;
  wbsat.sib = 3.93232680731067283;
  wbsat.sic = 0.253450103965147155;
  wbsat.sid = 0.856344826609515408;
  wbsat.ha  = 0.0;

#ifndef SILENT
  printf("<< Si-O Bond-Softening Function >>\n");
  printf("    m1         m2         m3         m4         m5\n");
  printf("-------------------------------------------------------\n");
  printf("%10.6f %10.6f %10.6f %10.6f %10.6f\n",
         wbsat.nc, wbsat.na, wbsat.nb, wbsat.nd, wbsat.ne);
  printf("-------------------------------------------------------\n\n");
  printf("<< Si-O Cutoff Function >>\n");
  printf("     R          D\n");
  printf("------------------------\n");
  printf("%10.6f %10.6f\n", wbcutsio.R, wbcutsio.D);
  printf("------------------------\n\n");
#endif

  double cut, a;

  for(i = 0; i < 2; i++) {
    for(j = 0; j < 2; j++) {
      cut = wb2[i][j].a;
      for(k = 0; k < 2; k++) {
        a = wb3[i][j][k].a[0];
        if(cut < a) cut = a;
        a = wb3[j][i][k].a[0];
        if(cut < a) cut = a;
        a = wb3[i][j][k].b[0];
        if(cut < a) cut = a;
        a = wb3[j][i][k].b[0];
        if(cut < a) cut = a;
      }
      cutoff[i][j] = cut;
    }
  }
}

void set_WATANABE3_parameter(void)
{
  int i,j,k;

  for(i=0;i<2;i++){
    for(j=0;j<2;j++){
      wb2[i][j].A = 0.0;
      wb2[i][j].B = 0.0;
      wb2[i][j].p = 0.0;
      wb2[i][j].q = 0.0;
      wb2[i][j].a = 1.0;
    }
  }

  wb2[Si][Si].A = 7.04955;
  wb2[Si][Si].B = 0.60222;
  wb2[Si][Si].p = 4.0;
  wb2[Si][Si].q = 0.0;
  wb2[Si][Si].a = 1.80;

  wb2[Si][O].A = 21.0;
  wb2[Si][O].B = 0.038;
  wb2[Si][O].p = 5.3;
  wb2[Si][O].q = -1.1;
  wb2[Si][O].a = 1.30;

  wb2[O][O].A = -12.29;
  wb2[O][O].B = 0.00;
  wb2[O][O].p = 0.00;
  wb2[O][O].q = 2.24;
  wb2[O][O].a = 1.25;

  wb2[O][Si].A = wb2[Si][O].A;
  wb2[O][Si].B = wb2[Si][O].B;
  wb2[O][Si].p = wb2[Si][O].p;
  wb2[O][Si].q = wb2[Si][O].q;
  wb2[O][Si].a = wb2[Si][O].a;

  for(i=0;i<2;i++){
    for(j=0;j<2;j++){
      for(k=0;k<2;k++){
        wb3[i][j][k].lam    = 0.0;
        wb3[i][j][k].gam[0] = 0.0;
        wb3[i][j][k].gam[1] = 0.0;
        wb3[i][j][k].a[0]   = 0.1;
        wb3[i][j][k].a[1]   = 0.1;
        wb3[i][j][k].mu     = 0.0;
        wb3[i][j][k].nu[0]  = 0.0;
        wb3[i][j][k].nu[1]  = 0.0;
        wb3[i][j][k].b[0]   = 0.1;
        wb3[i][j][k].b[1]   = 0.1;
        wb3[i][j][k].cos0   = 0.0;
        wb3[i][j][k].alpha  = 0.0;
        wb3[i][j][k].beta   = 0.0;
        wb3[i][j][k].g[0]   = 0.0;
        wb3[i][j][k].g[1]   = 0.0;
        wb3[i][j][k].g[2]   = 0.0;
      }
    }
  }
  wb3[Si][Si][Si].lam    = 5.878;
  wb3[Si][Si][Si].gam[0] = 1.01;
  wb3[Si][Si][Si].gam[1] = 1.01;
  wb3[Si][Si][Si].a[0]   = 1.80;
  wb3[Si][Si][Si].a[1]   = 1.80;
  wb3[Si][Si][Si].mu     = 4.0;
  wb3[Si][Si][Si].nu[0]  = 0.088;
  wb3[Si][Si][Si].nu[1]  = 0.088;
  wb3[Si][Si][Si].b[0]   = 1.2;
  wb3[Si][Si][Si].b[1]   = 1.2;
  wb3[Si][Si][Si].cos0   = -0.3333333;
  wb3[Si][Si][Si].alpha  = 0.3;
  wb3[Si][Si][Si].beta   = -1.35;
  wb3[Si][Si][Si].g[0]   = 0.6;
  wb3[Si][Si][Si].g[1]   = 6.0;
  wb3[Si][Si][Si].g[2]   = 0.0;

  wb3[Si][Si][O].lam    = 3.0;
  wb3[Si][Si][O].gam[0] = 0.032;
  wb3[Si][Si][O].gam[1] = 0.124;
  wb3[Si][Si][O].a[0]   = 1.15;
  wb3[Si][Si][O].a[1]   = 1.10;
  wb3[Si][Si][O].mu     = 5.878;
  wb3[Si][Si][O].nu[0]  = 1.01;
  wb3[Si][Si][O].nu[1]  = 1.082079;
  wb3[Si][Si][O].b[0]   = 1.8;
  wb3[Si][Si][O].b[1]   = 1.5;
  wb3[Si][Si][O].cos0   = -0.3333333;
  wb3[Si][Si][O].alpha  = 0.3;
  wb3[Si][Si][O].beta   = 0.3;
  wb3[Si][Si][O].g[0]   = 2.33;
  wb3[Si][Si][O].g[1]   = 1.9;
  wb3[Si][Si][O].g[2]   = 3.0;

  wb3[Si][O][Si].lam    = 2.50;
  wb3[Si][O][Si].gam[0] = 0.15;
  wb3[Si][O][Si].gam[1] = 0.15;
  wb3[Si][O][Si].a[0]   = 1.20;
  wb3[Si][O][Si].a[1]   = 1.20;
  wb3[Si][O][Si].mu     = 0.0;
  wb3[Si][O][Si].nu[0]  = 0.0;
  wb3[Si][O][Si].nu[1]  = 0.0;
  wb3[Si][O][Si].b[0]   = 0.1;
  wb3[Si][O][Si].b[1]   = 0.1;
  wb3[Si][O][Si].cos0   = -0.812;
  wb3[Si][O][Si].alpha  = 1.6;
  wb3[Si][O][Si].beta   = 1.6;
  wb3[Si][O][Si].g[0]   = 0.0;
  wb3[Si][O][Si].g[1]   = 0.0;
  wb3[Si][O][Si].g[2]   = 0.0;

  wb3[O][Si][Si].lam    = wb3[Si][Si][O].lam;
  wb3[O][Si][Si].gam[0] = wb3[Si][Si][O].gam[1];
  wb3[O][Si][Si].gam[1] = wb3[Si][Si][O].gam[0];
  wb3[O][Si][Si].a[0]   = wb3[Si][Si][O].a[1];
  wb3[O][Si][Si].a[1]   = wb3[Si][Si][O].a[0];
  wb3[O][Si][Si].mu     = wb3[Si][Si][O].mu;
  wb3[O][Si][Si].nu[0]  = wb3[Si][Si][O].nu[1];
  wb3[O][Si][Si].nu[1]  = wb3[Si][Si][O].nu[0];
  wb3[O][Si][Si].b[0]   = wb3[Si][Si][O].b[1];
  wb3[O][Si][Si].b[1]   = wb3[Si][Si][O].b[0];
  wb3[O][Si][Si].cos0   = wb3[Si][Si][O].cos0;
  wb3[O][Si][Si].alpha  = wb3[Si][Si][O].alpha;
  wb3[O][Si][Si].beta   = wb3[Si][Si][O].beta;
  wb3[O][Si][Si].g[0]   = wb3[Si][Si][O].g[0];
  wb3[O][Si][Si].g[1]   = wb3[Si][Si][O].g[1];
  wb3[O][Si][Si].g[2]   = wb3[Si][Si][O].g[2];

  wb3[O][Si][O].lam    = 10.5;
  wb3[O][Si][O].gam[0] = 0.310;
  wb3[O][Si][O].gam[1] = 0.310;
  wb3[O][Si][O].a[0]   = 1.10;
  wb3[O][Si][O].a[1]   = 1.10;
  wb3[O][Si][O].mu     = 0.0;
  wb3[O][Si][O].nu[0]  = 0.0;
  wb3[O][Si][O].nu[1]  = 0.0;
  wb3[O][Si][O].b[0]   = 0.5;
  wb3[O][Si][O].b[1]   = 0.5;
  wb3[O][Si][O].cos0   = -0.333333333;
  wb3[O][Si][O].alpha  = 0.0;
  wb3[O][Si][O].beta   = 0.0;
  wb3[O][Si][O].g[0]   = 0.0;
  wb3[O][Si][O].g[1]   = 0.0;
  wb3[O][Si][O].g[2]   = 0.0;

#ifndef SILENT
  printf("###### Constants in the potential function ######\n");
  printf("<< Tow-Body Term >>\n");
  printf(" i  j        Aij            Bij           pij        qij        aij\n");
  printf("----------------------------------------------------------------------\n");
  for(i = 0; i < 2; i++){
    for(j = i; j < 2; j++){
      printf("%s %s | %14.10f %14.10f %10.6f %10.6f %10.6f\n",
             charac[i].name, charac[j].name,
             wb2[i][j].A,wb2[i][j].B,wb2[i][j].p,
             wb2[i][j].q,wb2[i][j].a);
    }
  }
  printf("----------------------------------------------------------------------\n\n");
  printf("<< Three-Body Term >>\n");
  printf(" j  i  k      lam        gam^ij     gam^ik     a^ij       a^ik     cos(theta^0)\n");
  printf("-------------------------------------------------------------------------------\n");
  for(i = 0; i < 2; i++){
    for(j = 0; j < 2; j++){
      printf("%s %s %s | %10.6f %10.6f %10.6f %10.6f %10.6f %13.7f\n",
             charac[i].name, charac[i].name,
             charac[j].name,
             wb3[i][i][j].lam,wb3[i][i][j].gam[0],wb3[i][i][j].gam[1],
             wb3[i][i][j].a[0],wb3[i][i][j].a[1],wb3[i][i][j].cos0);
      if(i != j){
        printf("%s %s %s | %10.6f %10.6f %10.6f %10.6f %10.6f %13.7f\n",
               charac[i].name, charac[j].name,
               charac[i].name,
               wb3[i][j][i].lam,wb3[i][j][i].gam[0],wb3[i][j][i].gam[1],
               wb3[i][j][i].a[0],wb3[i][j][i].a[1],wb3[i][j][i].cos0);
      }
    }
  }
  printf("-------------------------------------------------------------------------------\n");
  printf(" j  i  k      mu         nu^ij      nu^ik      b^ij       b^ik  \n");
  printf("----------------------------------------------------------------\n");
  for(i = 0; i < 2; i++){
    for(j = 0; j < 2; j++){
      printf("%s %s %s | %10.6f %10.6f %10.6f %10.6f %10.6f\n",
             charac[i].name, charac[i].name,
             charac[j].name,
             wb3[i][i][j].mu,wb3[i][i][j].nu[0],wb3[i][i][j].nu[1],
             wb3[i][i][j].b[0],wb3[i][i][j].b[1]);
      if(i != j){
        printf("%s %s %s | %10.6f %10.6f %10.6f %10.6f %10.6f\n",
               charac[i].name, charac[j].name,
               charac[i].name,
               wb3[i][j][i].mu,wb3[i][j][i].nu[0],wb3[i][j][i].nu[1],
               wb3[i][j][i].b[0],wb3[i][j][i].b[1]);
      }
    }
  }

  printf("----------------------------------------------------------------\n");
  printf(" j  i  k      alpha      g0         g1         g2\n");
  printf("------------------------------------------------------\n");
  for(i = 0; i < 2; i++){
    for(j = 0; j < 2; j++){
      printf("%s %s %s | %10.6f %10.6f %10.6f %10.6f\n",
             charac[i].name, charac[i].name,
             charac[j].name,
             wb3[i][i][j].alpha,wb3[i][i][j].g[0],wb3[i][i][j].g[1],
             wb3[i][i][j].g[2]);
      if(i != j){
        printf("%s %s %s | %10.6f %10.6f %10.6f %10.6f\n",
               charac[i].name, charac[j].name,
               charac[i].name,
               wb3[i][j][i].alpha,wb3[i][j][i].g[0],wb3[i][j][i].g[1],
               wb3[i][j][i].g[2]);
      }
    }
  }
  printf("------------------------------------------------------\n\n");
#endif

  wbcutsio.R   = 1.20;
  wbcutsio.D   = 0.05;
  wbcutsio.RpD = wbcutsio.R + wbcutsio.D;

  wbcutsisi.R   = 1.55;
  wbcutsisi.D   = 0.05;
  wbcutsisi.RpD = wbcutsisi.R + wbcutsisi.D;

  wbsat.na =  0.712;
  wbsat.nb =  0.522;
  wbsat.nc =  5.274;
  wbsat.nd = -0.0372;
  wbsat.ne =  4.52;

  wbsat.sia = 0.303672895420570121;
  wbsat.sib = 3.93232680731067283;
  wbsat.sic = 0.253450103965147155;
  wbsat.sid = 0.856344826609515408;
  wbsat.ha  = 0.0;

#ifndef SILENT
  printf("<< Si-O Bond-Softening Function >>\n");
  printf("    m1         m2         m3         m4         m5\n");
  printf("-------------------------------------------------------\n");
  printf("%10.6f %10.6f %10.6f %10.6f %10.6f\n",
         wbsat.nc, wbsat.na, wbsat.nb, wbsat.nd, wbsat.ne);
  printf("-------------------------------------------------------\n\n");
  printf("<< Si-O Cutoff Function >>\n");
  printf("     R          D\n");
  printf("------------------------\n");
  printf("%10.6f %10.6f\n", wbcutsio.R, wbcutsio.D);
  printf("------------------------\n\n");
#endif

  double cut, a;

  for(i = 0; i < 2; i++) {
    for(j = 0; j < 2; j++) {
      cut = wb2[i][j].a;
      for(k = 0; k < 2; k++) {
        a = wb3[i][j][k].a[0];
        if(cut < a) cut = a;
        a = wb3[j][i][k].a[0];
        if(cut < a) cut = a;
        a = wb3[i][j][k].b[0];
        if(cut < a) cut = a;
        a = wb3[j][i][k].b[0];
        if(cut < a) cut = a;
      }
      cutoff[i][j] = cut;
    }
  }
}

void calc_WATANABE1_force(void)
{
   count_coordinate();
   two_body();
   three_body();
}

void calc_WATANABE2_force(void)
{
   count_coordinate();
   two_body2();
   three_body();
}

void calc_WATANABE3_force(void)
{
   count_coordinate();
   two_body3();
   three_body();
}

static void count_coordinate(void)
{
  int i,j,i1,neighbor;
  double tmp_co;

  for(i=0;i<N;i++){
    pl[i]->co=0.0;  /* initialize */
    pl[i]->self_co = 0.0;
  }
  for(i=0;i<N;i++){
    if(pl[i]->atmtype==O){
      neighbor=pl[i]->neighbor;
      for(i1=1;i1<=neighbor;i1++){
        j=pl[i]->j[i1];
        if(pl[j]->atmtype==Si){
          tmp_co = fc(r1(i,j), wbcutsio);
          pl[i]->co += tmp_co;
          pl[j]->co += tmp_co;
        }
      }
    }
    if(pl[i]->atmtype==Si){
      neighbor=pl[i]->neighbor;
      for(i1=1;i1<=neighbor;i1++){
        j=pl[i]->j[i1];
        if(pl[j]->atmtype==Si){
          tmp_co = fc(r1(i,j), wbcutsisi);
          pl[i]->self_co += tmp_co;
        }
      }
/*
      printf("%d selfco = %f\n",i,pl[i]->self_co);
*/
    }
  }
}

static void two_body(void)
{
  int i,j,i1;
  double dx,dy,dz;
  int neighbor;

  for(i=0;i<N;i++){
    neighbor = pl[i]->neighbor;
    for(i1=1;i1<=neighbor;i1++){ /* pl[i]->j[0] is "i" itself */
      j=pl[i]->j[i1];
      dx=pl[i]->x-pl[j]->x;
      dy=pl[i]->y-pl[j]->y;
      dz=pl[i]->z-pl[j]->z;
      mirror(&dx,&dy,&dz);
      scale_to_real(&dx,&dy,&dz);
      switch(pl[i]->atmtype){
        case(Si):{
          switch(pl[j]->atmtype){
            case(Si):{
              force_Si_Si(i,j,dx,dy,dz);
            } break;
            case(O):{
              force_Si_O(i,j,dx,dy,dz);
            } break;
	    case(Omol):{
              force_simple_simple(i,j,dx,dy,dz);
  	    } break;
          }
        } break;
        case(O):{
          switch(pl[j]->atmtype){
            case(Si):{
              force_O_Si(i,j,dx,dy,dz);
            } break;
            case(O):{
              force_O_O(i,j,dx,dy,dz);
            } break;
            case(Omol):{
              force_simple_simple(i,j,dx,dy,dz);
            } break;
          }
        } break;
        case(Omol):{
          force_simple_simple(i,j,dx,dy,dz);
        } break;
      }
    }
  }
}

static void two_body2(void)
{
  int i,j,i1;
  double dx,dy,dz;
  int neighbor;

  for(i=0;i<N;i++){
    neighbor = pl[i]->neighbor;
    for(i1=1;i1<=neighbor;i1++){ /* pl[i]->j[0] is "i" itself */
      j=pl[i]->j[i1];
      dx=pl[i]->x-pl[j]->x;
      dy=pl[i]->y-pl[j]->y;
      dz=pl[i]->z-pl[j]->z;
      mirror(&dx,&dy,&dz);
      scale_to_real(&dx,&dy,&dz);
      switch(pl[i]->atmtype){
        case(Si):{
          switch(pl[j]->atmtype){
            case(Si):{
              force_Si_Si(i,j,dx,dy,dz);
            } break;
            case(O):{
              force_Si_O2(i,j,dx,dy,dz);
            } break;
          }
        } break;
        case(O):{
          switch(pl[j]->atmtype){
            case(Si):{
              force_O_Si2(i,j,dx,dy,dz);
            } break;
            case(O):{
             force_O_O(i,j,dx,dy,dz);
            } break;
          }
        } break;
      }
    }
    switch(pl[i]->atmtype){
      case(Si):{
        force_Si_other(i);
      } break;
      case(O):{
        force_O_other(i);
      } break;
    }
  }
}

static void two_body3(void)
{
  int i,j,i1;
  double dx,dy,dz;
  int neighbor;

  for(i=0;i<N;i++){
    neighbor = pl[i]->neighbor;
    for(i1=1;i1<=neighbor;i1++){ /* pl[i]->j[0] is "i" itself */
      j=pl[i]->j[i1];
      dx=pl[i]->x-pl[j]->x;
      dy=pl[i]->y-pl[j]->y;
      dz=pl[i]->z-pl[j]->z;
      mirror(&dx,&dy,&dz);
      scale_to_real(&dx,&dy,&dz);
      switch(pl[i]->atmtype){
        case(Si):{
          switch(pl[j]->atmtype){
            case(Si):{
              force_Si_Si3(i,j,dx,dy,dz);
            } break;
            case(O):{
              force_Si_O(i,j,dx,dy,dz);
            } break;
          }
        } break;
        case(O):{
          switch(pl[j]->atmtype){
            case(Si):{
              force_O_Si(i,j,dx,dy,dz);
            } break;
            case(O):{
             force_O_O(i,j,dx,dy,dz);
            } break;
          }
        } break;
      }
    }
  }
}


static void three_body(void)
{
  int i,j,k,l,i1,i2,i3,neighbor,neighborj,trio;
  double dxij,dyij,dzij,dxik,dyik,dzik,dxjl,dyjl,dzjl,dxjk,dyjk,dzjk,dxil,dyil,
         dzil,rij,rik,rjl,rjk,inv_rij,inv_rik,inv_rjl,cosjik,cosijl;
  double ril, lamijl, eijl, thaijl, dthaijl, thbijl, dthbijl, dfc ,coskjl;
  double lamjik, ejik, thajik, dthajik, thbjik, dthbjik, deij, deik, dcosij, dcosik, dlamij;
  double *a,*gam,cos,alpha,beta;
  Watanabe3Body *parameter;
  double fij, fik, fil, f;
  double dlam_h, hkjl, c_c, inv_rija, inv_rika, inv_rjka, inv_rjla, inv_rila;
  double mu, *nu, *b, e2jik, e2ijl, de2ij, de2ik, inv_rijb, inv_rikb, inv_rjlb;
 
  for(i=0;i<N;i++){
    if(pl[i]->atmtype != Dummy){
      neighbor=pl[i]->neighbor;
      for(i1=1;i1<=neighbor;i1++){
        j=pl[i]->j[i1];
        if(pl[j]->atmtype != Dummy){
          rij=r1(i,j);
          for(i2=i1+1;i2<=neighbor;i2++){
            k=pl[i]->j[i2];
            if(pl[k]->atmtype != Dummy){
              rik=r1(i,k);
              parameter=&(wb3[pl[j]->atmtype][pl[i]->atmtype][pl[k]->atmtype]);
              a = parameter->a;
              b = parameter->b;
              if(((rij<a[0])&&(rik<a[1]))||((rij<b[0])&&(rik<b[1]))){
                dxij = pl[i]->x-pl[j]->x;
                dyij = pl[i]->y-pl[j]->y;
                dzij = pl[i]->z-pl[j]->z;
                mirror(&dxij,&dyij,&dzij);
                scale_to_real(&dxij,&dyij,&dzij);
                dxik = pl[i]->x-pl[k]->x;
                dyik = pl[i]->y-pl[k]->y;
                dzik = pl[i]->z-pl[k]->z;
                mirror(&dxik,&dyik,&dzik);
                scale_to_real(&dxik,&dyik,&dzik);
                gam = parameter->gam;
                cos = parameter->cos0;
                alpha = parameter->alpha;
                beta = parameter->beta;
                mu = parameter->mu;
                nu = parameter->nu;
                inv_rij = 1.0/rij;
                inv_rik = 1.0/rik;
                cosjik = (dxij*dxik+dyij*dyik+dzij*dzik)*inv_rij*inv_rik;
                c_c = cosjik-cos;
                lamjik = lambda(pl[i]->co, parameter);
                if((rij<a[0])&&(rik<a[1])){
                  inv_rija = 1.0/(rij-a[0]);
                  inv_rika = 1.0/(rik-a[1]);
                  ejik  = exp(gam[0]*inv_rija+gam[1]*inv_rika);
                } else {
                  inv_rija = 0.0;
                  inv_rika = 0.0;
                  ejik = 0.0;
                }
                if((rij<b[0])&&(rik<b[1])){
                  inv_rijb = 1.0/(rij-b[0]);
                  inv_rikb = 1.0/(rik-b[1]);
                  e2jik = mu*exp(nu[0]*inv_rijb+nu[1]*inv_rikb);
                } else {
                  inv_rijb = 0.0;
                  inv_rikb = 0.0;
                  e2jik = 0.0;
                }
                deij   = -gam[0]*inv_rija*inv_rija*ejik;
                deik   = -gam[1]*inv_rika*inv_rika*ejik;
                de2ij  = -nu[0]*inv_rijb*inv_rijb*e2jik;
                de2ik  = -nu[1]*inv_rikb*inv_rikb*e2jik;
                thajik = c_c*c_c+alpha*c_c*c_c*c_c; 
                thbjik = c_c*c_c+beta*c_c*c_c*c_c; 
                dthajik= 2.0*c_c+3.0*alpha*c_c*c_c;
                dthbjik= 2.0*c_c+3.0*beta*c_c*c_c;
                dcosij= (rij-rik*cosjik)*inv_rik*inv_rij;
                dcosik= (rik-rij*cosjik)*inv_rik*inv_rij;
                fij = -(lamjik*deij*thajik
                       + de2ij*thbjik+lamjik*ejik*dthajik*dcosij
                       + e2jik*dthbjik*dcosij)*inv_rij;
                fik = -(lamjik*deik*thajik
                       + de2ik*thbjik+lamjik*ejik*dthajik*dcosik
                       + e2jik*dthbjik*dcosik)*inv_rik;
                (pl[i]->pe3) += lamjik*ejik*thajik+e2jik*thbjik;
                (pl[i]->fx)  += fij*dxij+fik*dxik;
                (pl[i]->fy)  += fij*dyij+fik*dyik;
                (pl[i]->fz)  += fij*dzij+fik*dzik;
#ifdef DEBUG
                printf("(%d-%d-%d 3body : %d <- %d) [%f %f %f]\n",
                        j, i, k, i, j, fij*dxij, fij*dyij, fij*dzij);
                printf("(%d-%d-%d 3body : %d <- %d) [%f %f %f]\n",
                        j, i, k, i, k, fik*dxik, fik*dyik, fik*dzik);
#endif
                if(i < j){
                  sys.nfx[0][0] += fij*dxij*dxij;
                  sys.nfx[0][1] += fij*dxij*dyij;
                  sys.nfx[0][2] += fij*dxij*dzij;
                  sys.nfx[1][0] += fij*dyij*dxij;
                  sys.nfx[1][1] += fij*dyij*dyij;
                  sys.nfx[1][2] += fij*dyij*dzij;
                  sys.nfx[2][0] += fij*dzij*dxij;
                  sys.nfx[2][1] += fij*dzij*dyij;
                  sys.nfx[2][2] += fij*dzij*dzij;
                }
                if(i < k){
                  sys.nfx[0][0] += fik*dxik*dxik;
                  sys.nfx[0][1] += fik*dxik*dyik;
                  sys.nfx[0][2] += fik*dxik*dzik;
                  sys.nfx[1][0] += fik*dyik*dxik;
                  sys.nfx[1][1] += fik*dyik*dyik;
                  sys.nfx[1][2] += fik*dyik*dzik;
                  sys.nfx[2][0] += fik*dzik*dxik;
                  sys.nfx[2][1] += fik*dzik*dyik;
                  sys.nfx[2][2] += fik*dzik*dzik;
                }
                pl[i]->str[0][0] += 0.5*fij*dxij*dxij;
                pl[i]->str[0][1] += 0.5*fij*dxij*dyij;
                pl[i]->str[0][2] += 0.5*fij*dxij*dzij;
                pl[i]->str[1][0] += 0.5*fij*dyij*dxij;
                pl[i]->str[1][1] += 0.5*fij*dyij*dyij;
                pl[i]->str[1][2] += 0.5*fij*dyij*dzij;
                pl[i]->str[2][0] += 0.5*fij*dzij*dxij;
                pl[i]->str[2][1] += 0.5*fij*dzij*dyij;
                pl[i]->str[2][2] += 0.5*fij*dzij*dzij;
                pl[i]->str[0][0] += 0.5*fik*dxik*dxik;
                pl[i]->str[0][1] += 0.5*fik*dxik*dyik;
                pl[i]->str[0][2] += 0.5*fik*dxik*dzik;
                pl[i]->str[1][0] += 0.5*fik*dyik*dxik;
                pl[i]->str[1][1] += 0.5*fik*dyik*dyik;
                pl[i]->str[1][2] += 0.5*fik*dyik*dzik;
                pl[i]->str[2][0] += 0.5*fik*dzik*dxik;
                pl[i]->str[2][1] += 0.5*fik*dzik*dyik;
                pl[i]->str[2][2] += 0.5*fik*dzik*dzik;
                dlam_h = -1.0*d_lambda(pl[i]->co,parameter)*ejik*thajik;
                for(i3=1;i3<=neighbor;i3++){
                  l = pl[i]->j[i3];
                  if(pl[i]->atmtype != pl[l]->atmtype){
                    ril=r1(i,l);
                    dxil=pl[i]->x-pl[l]->x;
                    dyil=pl[i]->y-pl[l]->y;
                    dzil=pl[i]->z-pl[l]->z;
                    mirror(&dxil,&dyil,&dzil);
                    scale_to_real(&dxil,&dyil,&dzil);
                    f = dlam_h * d_fc(ril, wbcutsio) / ril;
                    (pl[i]->fx) += f*dxil;
                    (pl[i]->fy) += f*dyil;
                    (pl[i]->fz) += f*dzil;
                    if(i < l){
                      sys.nfx[0][0] += f*dxil*dxil;
                      sys.nfx[0][1] += f*dxil*dyil;
                      sys.nfx[0][2] += f*dxil*dzil;
                      sys.nfx[1][0] += f*dyil*dxil;
                      sys.nfx[1][1] += f*dyil*dyil;
                      sys.nfx[1][2] += f*dyil*dzil;
                      sys.nfx[2][0] += f*dzil*dxil;
                      sys.nfx[2][1] += f*dzil*dyil;
                      sys.nfx[2][2] += f*dzil*dzil;
                    }
                    pl[i]->str[0][0] += 0.5*f*dxil*dxil;
	            pl[i]->str[0][1] += 0.5*f*dxil*dyil;
	            pl[i]->str[0][2] += 0.5*f*dxil*dzil;
	            pl[i]->str[1][0] += 0.5*f*dyil*dxil;
	            pl[i]->str[1][1] += 0.5*f*dyil*dyil;
	            pl[i]->str[1][2] += 0.5*f*dyil*dzil;
                    pl[i]->str[2][0] += 0.5*f*dzil*dxil;
	            pl[i]->str[2][1] += 0.5*f*dzil*dyil;
	            pl[i]->str[2][2] += 0.5*f*dzil*dzil;
		  }
		}
	      }
            }
          }
        } 
      }
    }
  }
 
  for(i=0;i<N;i++){
    if(pl[i]->atmtype != Dummy){
      trio=pl[i]->trio;
      for(i3=1;i3<=trio;i3++){
        j=pl[i]->jl[i3][0];
        l=pl[i]->jl[i3][1];
        if(pl[j]->atmtype != Dummy && pl[l]->atmtype != Dummy){
          rij=r1(i,j);
          rjl=r1(j,l);
          parameter=&(wb3[pl[i]->atmtype][pl[j]->atmtype][pl[l]->atmtype]);
          a=parameter->a;
          b=parameter->b;
          if(((rij<a[0])&&(rjl<a[1]))||((rij<b[0])&&(rjl<b[1]))){
            dxij=pl[i]->x-pl[j]->x;
            dyij=pl[i]->y-pl[j]->y;
            dzij=pl[i]->z-pl[j]->z;
            mirror(&dxij,&dyij,&dzij);
            scale_to_real(&dxij,&dyij,&dzij);
            dxjl=pl[j]->x-pl[l]->x;
            dyjl=pl[j]->y-pl[l]->y;
            dzjl=pl[j]->z-pl[l]->z;
            mirror(&dxjl,&dyjl,&dzjl);
            scale_to_real(&dxjl,&dyjl,&dzjl);
            dxil=pl[i]->x-pl[l]->x;
            dyil=pl[i]->y-pl[l]->y;
            dzil=pl[i]->z-pl[l]->z;
            mirror(&dxil,&dyil,&dzil);
            scale_to_real(&dxil,&dyil,&dzil);
            gam=parameter->gam;
            cos=parameter->cos0;
            alpha=parameter->alpha;
            beta=parameter->beta;
            mu = parameter->mu;
            nu = parameter->nu;
            inv_rij=1.0/rij;
            inv_rjl=1.0/rjl;
            cosijl =0.0-(dxij*dxjl+dyij*dyjl+dzij*dzjl)*inv_rij*inv_rjl;
            c_c = cosijl-cos;
            lamijl = lambda(pl[j]->co, parameter);
            if((rij<a[0])&&(rjl<a[1])){
              inv_rija = 1.0/(rij-a[0]);
              inv_rjla = 1.0/(rjl-a[1]);
              eijl   = exp(gam[0]*inv_rija+gam[1]*inv_rjla);
            } else {
              inv_rija = 0.0;
              inv_rjla = 0.0;
              eijl = 0.0;
            }
            if((rij<b[0])&&(rjl<b[1])){
              inv_rijb = 1.0/(rij-b[0]);
              inv_rjlb = 1.0/(rjl-b[1]);
              e2ijl  = mu*exp(nu[0]*inv_rijb+nu[1]*inv_rjlb);
            } else {
              inv_rijb = 0.0;
              inv_rjlb = 0.0;
              e2ijl = 0.0;
            }
            thaijl  = c_c*c_c+alpha*c_c*c_c*c_c;
            dthaijl = 2.0*c_c+3.0*alpha*c_c*c_c;
            thbijl  = c_c*c_c+beta*c_c*c_c*c_c;
            dthbijl = 2.0*c_c+3.0*beta*c_c*c_c;
            deij   = -gam[0]*inv_rija*inv_rija*eijl;
            de2ij  = -nu[0]*inv_rijb*inv_rijb*e2ijl;
            dcosij = (rij-rjl*cosijl)*inv_rij*inv_rjl;
            if(pl[i]->atmtype != pl[j]->atmtype){
              dlamij = d_lambda(pl[j]->co, parameter)*d_fc(rij, wbcutsio);
            } else {
              dlamij = 0.0;
            }
            fij = -(dlamij*eijl*thaijl
                   + lamijl*deij*thaijl+de2ij*thbijl
                   + lamijl*eijl*dthaijl*dcosij+e2ijl*dthbijl*dcosij)*inv_rij;
            fil = (lamijl*eijl*dthaijl+e2ijl*dthbijl)*inv_rij*inv_rjl;
            (pl[i]->fx)+=fij*dxij+fil*dxil;
            (pl[i]->fy)+=fij*dyij+fil*dyil;
            (pl[i]->fz)+=fij*dzij+fil*dzil;
#ifdef DEBUG
            printf("(%d-%d-%d 3body : %d <- %d) [%f %f %f]\n",
                   i, j, l, i, j, fij*dxij, fij*dyij, fij*dzij);
            printf("(%d-%d-%d 3body : %d <- %d) [%f %f %f]\n",
                   i, j, l, i, l, fik*dxil, fik*dyil, fik*dzil);
#endif
            if(i < j){
              sys.nfx[0][0] += fij*dxij*dxij;
              sys.nfx[0][1] += fij*dxij*dyij;
              sys.nfx[0][2] += fij*dxij*dzij;
              sys.nfx[1][0] += fij*dyij*dxij;
              sys.nfx[1][1] += fij*dyij*dyij;
              sys.nfx[1][2] += fij*dyij*dzij;
              sys.nfx[2][0] += fij*dzij*dxij;
              sys.nfx[2][1] += fij*dzij*dyij;
              sys.nfx[2][2] += fij*dzij*dzij;
            } 
            if(i < l){
              sys.nfx[0][0] += fil*dxil*dxil;
              sys.nfx[0][1] += fil*dxil*dyil;
              sys.nfx[0][2] += fil*dxil*dzil;
              sys.nfx[1][0] += fil*dyil*dxil;
              sys.nfx[1][1] += fil*dyil*dyil;
              sys.nfx[1][2] += fil*dyil*dzil;
              sys.nfx[2][0] += fil*dzil*dxil;
              sys.nfx[2][1] += fil*dzil*dyil;
              sys.nfx[2][2] += fil*dzil*dzil;
            }
            pl[i]->str[0][0] += 0.5*fij*dxij*dxij;
            pl[i]->str[0][1] += 0.5*fij*dxij*dyij;
            pl[i]->str[0][2] += 0.5*fij*dxij*dzij;
            pl[i]->str[1][0] += 0.5*fij*dyij*dxij;
            pl[i]->str[1][1] += 0.5*fij*dyij*dyij;
            pl[i]->str[1][2] += 0.5*fij*dyij*dzij;
            pl[i]->str[2][0] += 0.5*fij*dzij*dxij;
            pl[i]->str[2][1] += 0.5*fij*dzij*dyij;
            pl[i]->str[2][2] += 0.5*fij*dzij*dzij;
            pl[i]->str[0][0] += 0.5*fil*dxil*dxil;
            pl[i]->str[0][1] += 0.5*fil*dxil*dyil;
            pl[i]->str[0][2] += 0.5*fil*dxil*dzil;
            pl[i]->str[1][0] += 0.5*fil*dyil*dxil;
            pl[i]->str[1][1] += 0.5*fil*dyil*dyil;
            pl[i]->str[1][2] += 0.5*fil*dyil*dzil;
            pl[i]->str[2][0] += 0.5*fil*dzil*dxil;
            pl[i]->str[2][1] += 0.5*fil*dzil*dyil;
            pl[i]->str[2][2] += 0.5*fil*dzil*dzil;
	  }
	}
      }
    }
  }

  for(i=0;i<N;i++){
    if(pl[i]->atmtype != Dummy){
      neighbor = pl[i]->neighbor;
      for(i1=1;i1<=neighbor;i1++){
        j = pl[i]->j[i1];
        if(pl[i]->atmtype != pl[j]->atmtype && pl[j]->atmtype != Dummy){
          rij = r1(i,j);
          dfc = d_fc(rij, wbcutsio);
          if(not_zero(dfc)){
            dxij=pl[i]->x-pl[j]->x;
            dyij=pl[i]->y-pl[j]->y;
            dzij=pl[i]->z-pl[j]->z;
            mirror(&dxij,&dyij,&dzij);
            scale_to_real(&dxij,&dyij,&dzij);
            neighborj = pl[j]->neighbor;
            for(i2=1;i2<neighborj;i2++){
              k=pl[j]->j[i2];
              if(k!=i){
                rjk=r1(j,k);
                for(i3=i2+1;i3<=neighborj;i3++){
                  l=pl[j]->j[i3];
                  if(l!=i){
                    rjl=r1(j,l);
                    parameter=
                        &(wb3[pl[k]->atmtype][pl[j]->atmtype][pl[l]->atmtype]);
                    a=parameter->a;
                    if((rjk<a[0])&&(rjl<a[1])){
                      dxjk=pl[j]->x-pl[k]->x;
                      dyjk=pl[j]->y-pl[k]->y;
                      dzjk=pl[j]->z-pl[k]->z;
                      mirror(&dxjk,&dyjk,&dzjk);
                      scale_to_real(&dxjk,&dyjk,&dzjk);
                      dxjl=pl[j]->x-pl[l]->x;
                      dyjl=pl[j]->y-pl[l]->y;
                      dzjl=pl[j]->z-pl[l]->z;
                      mirror(&dxjl,&dyjl,&dzjl);
                      scale_to_real(&dxjl,&dyjl,&dzjl);
                      gam=parameter->gam;
                      cos=parameter->cos0;
                      alpha=parameter->alpha;
                      coskjl=(dxjk*dxjl+dyjk*dyjl+dzjk*dzjl)/rjk/rjl;
                      inv_rjka= 1.0/(rjk-a[0]);
                      inv_rjla= 1.0/(rjl-a[1]);
                      c_c= coskjl-cos;
                      hkjl= exp(gam[0]*inv_rjka+gam[1]*inv_rjla)
                            *(c_c*c_c+alpha*c_c*c_c*c_c);
                      f= -1.0*d_lambda(pl[j]->co, parameter)*dfc*hkjl/rij;
                      (pl[i]->fx)+=f*dxij;
                      (pl[i]->fy)+=f*dyij;
                      (pl[i]->fz)+=f*dzij;
                      if(i < j){
                        sys.nfx[0][0] += f*dxij*dxij;
                        sys.nfx[0][1] += f*dxij*dyij;
                        sys.nfx[0][2] += f*dxij*dzij;
                        sys.nfx[1][0] += f*dyij*dxij;
                        sys.nfx[1][1] += f*dyij*dyij;
                        sys.nfx[1][2] += f*dyij*dzij;
                        sys.nfx[2][0] += f*dzij*dxij;
                        sys.nfx[2][1] += f*dzij*dyij;
                        sys.nfx[2][2] += f*dzij*dzij;
                      }
                      pl[i]->str[0][0] += 0.5*f*dxij*dxij;
                      pl[i]->str[0][1] += 0.5*f*dxij*dyij;
                      pl[i]->str[0][2] += 0.5*f*dxij*dzij;
                      pl[i]->str[1][0] += 0.5*f*dyij*dxij;
                      pl[i]->str[1][1] += 0.5*f*dyij*dyij;
                      pl[i]->str[1][2] += 0.5*f*dyij*dzij;
                      pl[i]->str[2][0] += 0.5*f*dzij*dxij;
                      pl[i]->str[2][1] += 0.5*f*dzij*dyij;
                      pl[i]->str[2][2] += 0.5*f*dzij*dzij;
		    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

static void force_Si_Si(int i, int j,double dxij, double dyij, double dzij)
{
  force_simple_simple(i,j,dxij,dyij,dzij);
}

static void force_O_O(int i, int j,double dxij, double dyij, double dzij)
{
  force_simple_simple(i,j,dxij,dyij,dzij);
}

static void force_simple_simple(int i, int j,
                                double dxij, double dyij, double dzij)
{
  double A,B,p,q,a,r,f,tmp[7];
  AtomIonType ai, aj;

  r=r1(i,j);
  ai = pl[i]->atmtype;
  aj = pl[j]->atmtype;
  if(ai == Omol) ai = Ar;
  if(aj == Omol) aj = Ar;
  if(ai != Dummy && aj != Dummy){
    a=wb2[ai][aj].a;
    if(r<a){
      A=wb2[ai][aj].A;
      B=wb2[ai][aj].B;
      p=wb2[ai][aj].p;
      q=wb2[ai][aj].q;
      tmp[0]=1.0/(r-a);
      tmp[1]=exp(tmp[0]);
      tmp[2]=A*(B*pow(r,-p)-pow(r,-q));
      tmp[3]=tmp[0]*tmp[0];
      tmp[4]=A*(B*p*pow(r,-p-1.0)-q*pow(r,-q-1.0));
      tmp[5]=tmp[1]*tmp[2];
      tmp[6]= tmp[5]*0.5;
      (pl[i]->pe2)+=tmp[6];
      f=(tmp[4]+tmp[2]*tmp[3])*tmp[1]/r;
      (pl[i]->fx)+=f*dxij;
      (pl[i]->fy)+=f*dyij;
      (pl[i]->fz)+=f*dzij;
#ifdef DEBUG
      printf("(%d-%d 2body : %d <- %d) [%f %f %f]\n",
              i, j, i, j, f*dxij, f*dyij, f*dzij);
#endif
      if(i < j){
        sys.nfx[0][0] += f * dxij * dxij;
        sys.nfx[0][1] += f * dxij * dyij;
        sys.nfx[0][2] += f * dxij * dzij;
        sys.nfx[1][0] += f * dyij * dxij;
        sys.nfx[1][1] += f * dyij * dyij;
        sys.nfx[1][2] += f * dyij * dzij;
        sys.nfx[2][0] += f * dzij * dxij;
        sys.nfx[2][1] += f * dzij * dyij;
        sys.nfx[2][2] += f * dzij * dzij;
      }
      pl[i]->str[0][0] += 0.5*f*dxij*dxij;
      pl[i]->str[0][1] += 0.5*f*dxij*dyij;
      pl[i]->str[0][2] += 0.5*f*dxij*dzij;
      pl[i]->str[1][0] += 0.5*f*dyij*dxij;
      pl[i]->str[1][1] += 0.5*f*dyij*dyij;
      pl[i]->str[1][2] += 0.5*f*dyij*dzij;
      pl[i]->str[2][0] += 0.5*f*dzij*dxij;
      pl[i]->str[2][1] += 0.5*f*dzij*dyij;
      pl[i]->str[2][2] += 0.5*f*dzij*dzij;
    }
  }
}

static void force_Si_Si3(int i, int j, double dxij, double dyij, double dzij)
{
  double A,B,p,q,a,rij,rjk,f,gi,gj,gk,dgj,dfc,dgi,rik,dxik,dyik,dzik;
  int i2, neighbor,k;
  double invrija, ex, sw, dex, dsw;

  f = 0.0;
  rij=r1(i,j);
  a=wb2[Si][Si].a;
  A=wb2[Si][Si].A;
  B=wb2[Si][Si].B;
  p=wb2[Si][Si].p;
  q=wb2[Si][Si].q;
  if(rij<a){
    dxij=pl[i]->x-pl[j]->x;
    dyij=pl[i]->y-pl[j]->y;
    dzij=pl[i]->z-pl[j]->z;
    mirror(&dxij,&dyij,&dzij);
    scale_to_real(&dxij,&dyij,&dzij);
    gi = bsSi(pl[i]->self_co);
    gj = bsSi(pl[j]->self_co);
    invrija=1.0/(rij-a);
    ex=exp(invrija);
    sw=A*(B*pow(rij,-p)-pow(rij,-q));
    dex= -invrija*invrija*ex;
    dsw= -A*(B*p*pow(rij,-p-1.0)-q*pow(rij,-q-1.0));
    (pl[i]->pe2)+= gi*gj*sw*ex*0.5;
    f -= gi*gj*(dsw*ex + sw*dex);
    dgj = d_bsSi(pl[j]->self_co);
    dfc = d_fc(rij, wbcutsisi);
    f -= gi*dgj*dfc*sw*ex;
    f /= rij;
    (pl[i]->fx)+=f*dxij;
    (pl[i]->fy)+=f*dyij;
    (pl[i]->fz)+=f*dzij;
    if(i < j){
      sys.nfx[0][0] += f * dxij * dxij;
      sys.nfx[0][1] += f * dxij * dyij;
      sys.nfx[0][2] += f * dxij * dzij;
      sys.nfx[1][0] += f * dyij * dxij;
      sys.nfx[1][1] += f * dyij * dyij;
      sys.nfx[1][2] += f * dyij * dzij;
      sys.nfx[2][0] += f * dzij * dxij;
      sys.nfx[2][1] += f * dzij * dyij;
      sys.nfx[2][2] += f * dzij * dzij;
    }
    pl[i]->str[0][0] += 0.5*f*dxij*dxij;
    pl[i]->str[0][1] += 0.5*f*dxij*dyij;
    pl[i]->str[0][2] += 0.5*f*dxij*dzij;
    pl[i]->str[1][0] += 0.5*f*dyij*dxij;
    pl[i]->str[1][1] += 0.5*f*dyij*dyij;
    pl[i]->str[1][2] += 0.5*f*dyij*dzij;
    pl[i]->str[2][0] += 0.5*f*dzij*dxij;
    pl[i]->str[2][1] += 0.5*f*dzij*dyij;
    pl[i]->str[2][2] += 0.5*f*dzij*dzij;

    dgi = d_bsSi(pl[i]->self_co);
    neighbor=pl[i]->neighbor;
    for(i2=1;i2<=neighbor;i2++){
      k=pl[i]->j[i2];
      if(pl[k]->atmtype==Si){
        rik=r1(i,k);
        dfc=d_fc(rik, wbcutsisi);
        dxik=pl[i]->x-pl[k]->x;
        dyik=pl[i]->y-pl[k]->y;
        dzik=pl[i]->z-pl[k]->z;
        mirror(&dxik,&dyik,&dzik);
        scale_to_real(&dxik,&dyik,&dzik);
        f = -dgi*dfc*gj*sw*ex/rik;
        (pl[i]->fx)+=f*dxik;
        (pl[i]->fy)+=f*dyik;
        (pl[i]->fz)+=f*dzik;
        if(i < k){
          sys.nfx[0][0] += f * dxik * dxik;
          sys.nfx[0][1] += f * dxik * dyik;
          sys.nfx[0][2] += f * dxik * dzik;
          sys.nfx[1][0] += f * dyik * dxik;
          sys.nfx[1][1] += f * dyik * dyik;
          sys.nfx[1][2] += f * dyik * dzik;
          sys.nfx[2][0] += f * dzik * dxik;
          sys.nfx[2][1] += f * dzik * dyik;
          sys.nfx[2][2] += f * dzik * dzik;
        }
        pl[i]->str[0][0] += 0.5*f*dxik*dxik;
        pl[i]->str[0][1] += 0.5*f*dxik*dyik;
        pl[i]->str[0][2] += 0.5*f*dxik*dzik;
        pl[i]->str[1][0] += 0.5*f*dyik*dxik;
        pl[i]->str[1][1] += 0.5*f*dyik*dyik;
        pl[i]->str[1][2] += 0.5*f*dyik*dzik;
        pl[i]->str[2][0] += 0.5*f*dzik*dxik;
        pl[i]->str[2][1] += 0.5*f*dzik*dyik;
        pl[i]->str[2][2] += 0.5*f*dzik*dzik;
      }
    }
  }
  if(fc(rij, wbcutsisi) > 0.0){
    neighbor=pl[j]->neighbor;
    dfc = d_fc(rij, wbcutsisi);
    for(i2=1;i2<=neighbor;i2++){
      k=pl[j]->j[i2];
      if((pl[k]->atmtype == Si)&&(i != k)){
        rjk = r1(j,k);
        if(rjk < a){
          dgj = d_bsSi(pl[j]->self_co);
          gk = bsSi(pl[k]->self_co);
          sw = A*(B*pow(rjk,-p)-pow(rjk,-q))*exp(1.0/(rjk - a));
          f = -(dgj*dfc*gk*sw)/rij;
          (pl[i]->fx)+=f*dxij;
          (pl[i]->fy)+=f*dyij;
          (pl[i]->fz)+=f*dzij;
          if(i < j){
            sys.nfx[0][0] += f * dxij * dxij;
            sys.nfx[0][1] += f * dxij * dyij;
            sys.nfx[0][2] += f * dxij * dzij;
            sys.nfx[1][0] += f * dyij * dxij;
            sys.nfx[1][1] += f * dyij * dyij;
            sys.nfx[1][2] += f * dyij * dzij;
            sys.nfx[2][0] += f * dzij * dxij;
            sys.nfx[2][1] += f * dzij * dyij;
            sys.nfx[2][2] += f * dzij * dzij;
          }
          pl[i]->str[0][0] += 0.5*f*dxij*dxij;
          pl[i]->str[0][1] += 0.5*f*dxij*dyij;
          pl[i]->str[0][2] += 0.5*f*dxij*dzij;
          pl[i]->str[1][0] += 0.5*f*dyij*dxij;
          pl[i]->str[1][1] += 0.5*f*dyij*dyij;
          pl[i]->str[1][2] += 0.5*f*dyij*dzij;
          pl[i]->str[2][0] += 0.5*f*dzij*dxij;
          pl[i]->str[2][1] += 0.5*f*dzij*dyij;
          pl[i]->str[2][2] += 0.5*f*dzij*dzij;
        }
      }
    }
  }
}

static void force_Si_O(int i, int j, double dxij, double dyij, double dzij)
{
  double A,B,p,q,a,rij,rjk,rik,g,g2,dg,dg2,dfc,f,dxik,dyik,dzik,tmp[8];
  int i2,j1,k,neighbor;

  f=0.0; /* initialize */
  rij=r1(i,j);
  a=wb2[Si][O].a;
  A=wb2[Si][O].A;
  B=wb2[Si][O].B;
  p=wb2[Si][O].p;
  q=wb2[Si][O].q;
  if(rij<a){
    tmp[0]=1.0/(rij-a);
    tmp[1]=exp(tmp[0]);
    tmp[2]=A*(B*pow(rij,-p)-pow(rij,-q));
    tmp[3]=tmp[0]*tmp[0];
    tmp[4]=A*(B*p*pow(rij,-p-1.0)-q*pow(rij,-q-1.0));
    tmp[5]=tmp[1]*tmp[2];
    g=gO(pl[j]->co);
    g2 = gSi(pl[i]->co);
    tmp[6]= g*g2*tmp[5]*0.5;
    (pl[i]->pe2)+=tmp[6];
    f=g*g2*(tmp[4]+tmp[2]*tmp[3])*tmp[1]; /* be devided by rij before use */
  }
  dg=d_gO(pl[j]->co);
  dfc=d_fc(rij, wbcutsio);
  tmp[0]=dg*dfc;
  neighbor=pl[j]->neighbor;
  for(j1=1;j1<=neighbor;j1++){
    k=pl[j]->j[j1];
    if(pl[k]->atmtype==Si){
      rjk=r1(j,k);
      if(rjk<a){
        g2 = gSi(pl[k]->co);
        tmp[1]=A*(B*pow(rjk,-p)-pow(rjk,-q));
        tmp[2]=exp(1.0/(rjk-a));
        f-=g2*tmp[0]*tmp[1]*tmp[2];
      }
    }
  }
  f=f/rij;
  (pl[i]->fx)+=f*dxij;
  (pl[i]->fy)+=f*dyij;
  (pl[i]->fz)+=f*dzij;

  if(i < j){
    sys.nfx[0][0] += f * dxij * dxij;
    sys.nfx[0][1] += f * dxij * dyij;
    sys.nfx[0][2] += f * dxij * dzij;
    sys.nfx[1][0] += f * dyij * dxij;
    sys.nfx[1][1] += f * dyij * dyij;
    sys.nfx[1][2] += f * dyij * dzij;
    sys.nfx[2][0] += f * dzij * dxij;
    sys.nfx[2][1] += f * dzij * dyij;
    sys.nfx[2][2] += f * dzij * dzij;
  }
  pl[i]->str[0][0] += 0.5*f*dxij*dxij;
  pl[i]->str[0][1] += 0.5*f*dxij*dyij;
  pl[i]->str[0][2] += 0.5*f*dxij*dzij;
  pl[i]->str[1][0] += 0.5*f*dyij*dxij;
  pl[i]->str[1][1] += 0.5*f*dyij*dyij;
  pl[i]->str[1][2] += 0.5*f*dyij*dzij;
  pl[i]->str[2][0] += 0.5*f*dzij*dxij;
  pl[i]->str[2][1] += 0.5*f*dzij*dyij;
  pl[i]->str[2][2] += 0.5*f*dzij*dzij;

  g = gO(pl[j]->co);
  if(rij<a){
  } else {
    tmp[5] = 0.0;
  }
  dg2 = d_gSi(pl[i]->co);
  tmp[0] = g * dg2 * tmp[5];
  neighbor = pl[i]->neighbor;
  for(i2=1;i2<=neighbor;i2++){
    k=pl[i]->j[i2];
    if(pl[k]->atmtype==O){
      rik=r1(i,k);
      dfc=d_fc(rik, wbcutsio);
      dxik=pl[i]->x-pl[k]->x;
      dyik=pl[i]->y-pl[k]->y;
      dzik=pl[i]->z-pl[k]->z;
      mirror(&dxik,&dyik,&dzik);
      scale_to_real(&dxik,&dyik,&dzik);
      f=0.0-tmp[0]*dfc/rik;
      (pl[i]->fx)+=f*dxik;
      (pl[i]->fy)+=f*dyik;
      (pl[i]->fz)+=f*dzik;
      if(i < k){
        sys.nfx[0][0] += f * dxik * dxik;
        sys.nfx[0][1] += f * dxik * dyik;
        sys.nfx[0][2] += f * dxik * dzik;
        sys.nfx[1][0] += f * dyik * dxik;
        sys.nfx[1][1] += f * dyik * dyik;
        sys.nfx[1][2] += f * dyik * dzik;
        sys.nfx[2][0] += f * dzik * dxik;
        sys.nfx[2][1] += f * dzik * dyik;
        sys.nfx[2][2] += f * dzik * dzik;
      }
      pl[i]->str[0][0] += 0.5*f*dxik*dxik;
      pl[i]->str[0][1] += 0.5*f*dxik*dyik;
      pl[i]->str[0][2] += 0.5*f*dxik*dzik;
      pl[i]->str[1][0] += 0.5*f*dyik*dxik;
      pl[i]->str[1][1] += 0.5*f*dyik*dyik;
      pl[i]->str[1][2] += 0.5*f*dyik*dzik;
      pl[i]->str[2][0] += 0.5*f*dzik*dxik;
      pl[i]->str[2][1] += 0.5*f*dzik*dyik;
      pl[i]->str[2][2] += 0.5*f*dzik*dzik;
    }
  }
}

static void force_Si_O2(int i, int j, double dxij, double dyij, double dzij)
{
  double A,B,p,q,a,rij,rjk,rik,go,gsi,dgo,dgsi,dfc,f,dxik,dyik,dzik,R,F,dFr,dFR,dR,fcore;
  int i2,j1,k,neighbor;
  double invrija, ex, sw, dex, dsw;

  f=0.0; /* initialize */
  rij=r1(i,j);
  a=wb2[Si][O].a;
  A=wb2[Si][O].A;
  B=wb2[Si][O].B;
  p=wb2[Si][O].p;
  q=wb2[Si][O].q;
  R = Rc(pl[j]->co);
  if((rij<a)&&(rij < R+0.05)){
    dxij=pl[i]->x-pl[j]->x;
    dyij=pl[i]->y-pl[j]->y;
    dzij=pl[i]->z-pl[j]->z;
    mirror(&dxij,&dyij,&dzij);
    scale_to_real(&dxij,&dyij,&dzij);
    invrija = 1.0/(rij-a);
    ex = exp(invrija);
    dex = -invrija*invrija*ex; 
    sw = A*(B*pow(rij,-p)-pow(rij,-q));
    dsw = -A*(B*p*pow(rij,-p-1.0)-q*pow(rij,-q-1.0));
    go = gO(pl[j]->co);
    gsi = gSi(pl[i]->co);
    F = Fc(rij, R);
    dFr = d_Fc_r(rij, R);

    (pl[i]->pe2)+= go*gsi*F*sw*ex*0.5;
    f = -go*gsi*F*(sw*dex + dsw*ex); /* be devided by rij before use */

    dgo = d_gO(pl[j]->co);
    dfc = d_fc(rij, wbcutsio);

    f -= gsi*dgo*dfc*F*ex*sw;

    f -= gsi*go*dFr*ex*sw;

    dFR = d_Fc_R(rij, R);
    dR  = d_Rc(pl[j]->co);
   
    f -= gsi*go*sw*ex*dFR*dR*dfc;

    f /= rij;
    (pl[i]->fx)+=f*dxij;
    (pl[i]->fy)+=f*dyij;
    (pl[i]->fz)+=f*dzij;

    if(i < j){
      sys.nfx[0][0] += f * dxij * dxij;
      sys.nfx[0][1] += f * dxij * dyij;
      sys.nfx[0][2] += f * dxij * dzij;
      sys.nfx[1][0] += f * dyij * dxij;
      sys.nfx[1][1] += f * dyij * dyij;
      sys.nfx[1][2] += f * dyij * dzij;
      sys.nfx[2][0] += f * dzij * dxij;
      sys.nfx[2][1] += f * dzij * dyij;
      sys.nfx[2][2] += f * dzij * dzij;
    }
    pl[i]->str[0][0] += 0.5*f*dxij*dxij;
    pl[i]->str[0][1] += 0.5*f*dxij*dyij;
    pl[i]->str[0][2] += 0.5*f*dxij*dzij;
    pl[i]->str[1][0] += 0.5*f*dyij*dxij;
    pl[i]->str[1][1] += 0.5*f*dyij*dyij;
    pl[i]->str[1][2] += 0.5*f*dyij*dzij;
    pl[i]->str[2][0] += 0.5*f*dzij*dxij;
    pl[i]->str[2][1] += 0.5*f*dzij*dyij;
    pl[i]->str[2][2] += 0.5*f*dzij*dzij;

    neighbor = pl[i]->neighbor;
    dgsi = d_gSi(pl[i]->co);
    fcore = -go*F*sw*ex*dgsi;
    for(i2=1;i2<=neighbor;i2++){
      k=pl[i]->j[i2];
      if(pl[k]->atmtype==O){
        rik=r1(i,k);
        dfc=d_fc(rik, wbcutsio);
        dxik=pl[i]->x-pl[k]->x;
        dyik=pl[i]->y-pl[k]->y;
        dzik=pl[i]->z-pl[k]->z;
        mirror(&dxik,&dyik,&dzik);
        scale_to_real(&dxik,&dyik,&dzik);
        f = fcore*dfc/rik;
        (pl[i]->fx)+=f*dxik;
        (pl[i]->fy)+=f*dyik;
        (pl[i]->fz)+=f*dzik;
        if(i < k){
          sys.nfx[0][0] += f * dxik * dxik;
          sys.nfx[0][1] += f * dxik * dyik;
          sys.nfx[0][2] += f * dxik * dzik;
          sys.nfx[1][0] += f * dyik * dxik;
          sys.nfx[1][1] += f * dyik * dyik;
          sys.nfx[1][2] += f * dyik * dzik;
          sys.nfx[2][0] += f * dzik * dxik;
          sys.nfx[2][1] += f * dzik * dyik;
          sys.nfx[2][2] += f * dzik * dzik;
        }
        pl[i]->str[0][0] += 0.5*f*dxik*dxik;
        pl[i]->str[0][1] += 0.5*f*dxik*dyik;
        pl[i]->str[0][2] += 0.5*f*dxik*dzik;
        pl[i]->str[1][0] += 0.5*f*dyik*dxik;
        pl[i]->str[1][1] += 0.5*f*dyik*dyik;
        pl[i]->str[1][2] += 0.5*f*dyik*dzik;
        pl[i]->str[2][0] += 0.5*f*dzik*dxik;
        pl[i]->str[2][1] += 0.5*f*dzik*dyik;
        pl[i]->str[2][2] += 0.5*f*dzik*dzik;
      }
    }
  }
}

static void force_O_Si(int i, int j,double dxij, double dyij, double dzij)
{
  double A,B,p,q,a,rij,rik,rjk,dxik,dyik,dzik,f,g,g2,dg,dg2,dfc,tmp[7];
  int i2,k,neighbor;
  double tmp_fx, tmp_fy, tmp_fz;

  f=0.0;  /* initialize */
  tmp_fx = 0.0;
  tmp_fy = 0.0;
  tmp_fz = 0.0;
  rij=r1(i,j);
  a=wb2[O][Si].a;
  A=wb2[O][Si].A;
  B=wb2[O][Si].B;
  p=wb2[O][Si].p;
  q=wb2[O][Si].q;
  if(rij<a){
    g=gO(pl[i]->co);
    g2 = gSi(pl[j]->co);
    tmp[0]=1.0/(rij-a);
    tmp[1]=exp(tmp[0]);
    tmp[2]=A*(B*pow(rij,-p)-pow(rij,-q));
    tmp[3]=tmp[0]*tmp[0];
    tmp[4]=A*(B*p*pow(rij,-p-1.0)-q*pow(rij,-q-1.0));
    tmp[5]=tmp[1]*tmp[2];
    tmp[6]= g * g2 * tmp[5]*0.5;
    (pl[i]->pe2)+=tmp[6];
    f = g * g2 * (tmp[4]+tmp[2]*tmp[3])*tmp[1]/rij;
    (pl[i]->fx)+=f*dxij;
    (pl[i]->fy)+=f*dyij;
    (pl[i]->fz)+=f*dzij;

    if(i < j){
      sys.nfx[0][0] += f * dxij * dxij;
      sys.nfx[0][1] += f * dxij * dyij;
      sys.nfx[0][2] += f * dxij * dzij;
      sys.nfx[1][0] += f * dyij * dxij;
      sys.nfx[1][1] += f * dyij * dyij;
      sys.nfx[1][2] += f * dyij * dzij;
      sys.nfx[2][0] += f * dzij * dxij;
      sys.nfx[2][1] += f * dzij * dyij;
      sys.nfx[2][2] += f * dzij * dzij;
    }
    pl[i]->str[0][0] += 0.5*f*dxij*dxij;
    pl[i]->str[0][1] += 0.5*f*dxij*dyij;
    pl[i]->str[0][2] += 0.5*f*dxij*dzij;
    pl[i]->str[1][0] += 0.5*f*dyij*dxij;
    pl[i]->str[1][1] += 0.5*f*dyij*dyij;
    pl[i]->str[1][2] += 0.5*f*dyij*dzij;
    pl[i]->str[2][0] += 0.5*f*dzij*dxij;
    pl[i]->str[2][1] += 0.5*f*dzij*dyij;
    pl[i]->str[2][2] += 0.5*f*dzij*dzij;

    g2 = gSi(pl[j]->co);
    dg=d_gO(pl[i]->co);
    neighbor=pl[i]->neighbor;
    for(i2=1;i2<=neighbor;i2++){
      k=pl[i]->j[i2];
      if(pl[k]->atmtype==Si){
        rik=r1(i,k);
        dfc=d_fc(rik, wbcutsio);
        dxik=pl[i]->x-pl[k]->x;
        dyik=pl[i]->y-pl[k]->y;
        dzik=pl[i]->z-pl[k]->z;
        mirror(&dxik,&dyik,&dzik);
        scale_to_real(&dxik,&dyik,&dzik);
        f=0.0-g2*dg*dfc*tmp[5]/rik;
        (pl[i]->fx)+=f*dxik;
        (pl[i]->fy)+=f*dyik;
        (pl[i]->fz)+=f*dzik;
        if(i < k){
          sys.nfx[0][0] += f * dxik * dxik;
          sys.nfx[0][1] += f * dxik * dyik;
          sys.nfx[0][2] += f * dxik * dzik;
          sys.nfx[1][0] += f * dyik * dxik;
          sys.nfx[1][1] += f * dyik * dyik;
          sys.nfx[1][2] += f * dyik * dzik;
          sys.nfx[2][0] += f * dzik * dxik;
          sys.nfx[2][1] += f * dzik * dyik;
          sys.nfx[2][2] += f * dzik * dzik;
        }
        pl[i]->str[0][0] += 0.5*f*dxik*dxik;
        pl[i]->str[0][1] += 0.5*f*dxik*dyik;
        pl[i]->str[0][2] += 0.5*f*dxik*dzik;
        pl[i]->str[1][0] += 0.5*f*dyik*dxik;
        pl[i]->str[1][1] += 0.5*f*dyik*dyik;
        pl[i]->str[1][2] += 0.5*f*dyik*dzik;
        pl[i]->str[2][0] += 0.5*f*dzik*dxik;
        pl[i]->str[2][1] += 0.5*f*dzik*dyik;
        pl[i]->str[2][2] += 0.5*f*dzik*dzik;
      }
    }
  }

  f = 0.0;
  dg2 = d_gSi(pl[j]->co);
  dfc = d_fc(rij, wbcutsio);
  neighbor =pl[j]->neighbor;
  for(i2=1;i2<=neighbor;i2++){
    k=pl[j]->j[i2];
    if(pl[k]->atmtype==O){
      rjk=r1(j,k);
      if(rjk<a){
        g = gO(pl[k]->co);
        tmp[1]=A*(B*pow(rjk,-p)-pow(rjk,-q));
        tmp[2]=exp(1.0/(rjk-a));
        f-=g*tmp[1]*tmp[2];
      }
    }
  }
  f=f*dg2*dfc/rij;
  (pl[i]->fx)+=f*dxij;
  (pl[i]->fy)+=f*dyij;
  (pl[i]->fz)+=f*dzij;
  if(i < j){
    sys.nfx[0][0] += f * dxij * dxij;
    sys.nfx[0][1] += f * dxij * dyij;
    sys.nfx[0][2] += f * dxij * dzij;
    sys.nfx[1][0] += f * dyij * dxij;
    sys.nfx[1][1] += f * dyij * dyij;
    sys.nfx[1][2] += f * dyij * dzij;
    sys.nfx[2][0] += f * dzij * dxij;
    sys.nfx[2][1] += f * dzij * dyij;
    sys.nfx[2][2] += f * dzij * dzij;
  }
  pl[i]->str[0][0] += 0.5*f*dxij*dxij;
  pl[i]->str[0][1] += 0.5*f*dxij*dyij;
  pl[i]->str[0][2] += 0.5*f*dxij*dzij;
  pl[i]->str[1][0] += 0.5*f*dyij*dxij;
  pl[i]->str[1][1] += 0.5*f*dyij*dyij;
  pl[i]->str[1][2] += 0.5*f*dyij*dzij;
  pl[i]->str[2][0] += 0.5*f*dzij*dxij;
  pl[i]->str[2][1] += 0.5*f*dzij*dyij;
  pl[i]->str[2][2] += 0.5*f*dzij*dzij;
}

static void force_O_Si2(int i, int j,double dxij, double dyij, double dzij)
{
  double A,B,p,q,a,rij,rjk,rik,go,gsi,dgo,dgsi,dfc,f,dxik,dyik,dzik,R,F,dFr,dFR,dR,fcore;
  int i2,j1,k,neighbor;
  double invrija, ex, sw, dex, dsw;

  f=0.0;  /* initialize */
  rij=r1(i,j);
  a=wb2[O][Si].a;
  A=wb2[O][Si].A;
  B=wb2[O][Si].B;
  p=wb2[O][Si].p;
  q=wb2[O][Si].q;
  R = Rc(pl[i]->co);
  if((rij<a)&&(rij<R+0.05)){
    dxij=pl[i]->x-pl[j]->x;
    dyij=pl[i]->y-pl[j]->y;
    dzij=pl[i]->z-pl[j]->z;
    mirror(&dxij,&dyij,&dzij);
    scale_to_real(&dxij,&dyij,&dzij);

    go  = gO(pl[i]->co);
    gsi = gSi(pl[j]->co);
    F   = Fc(rij, R);
    invrija=1.0/(rij-a);
    ex=exp(invrija);
    sw=A*(B*pow(rij,-p)-pow(rij,-q));
    dex= -invrija*invrija*ex;
    dsw= -A*(B*p*pow(rij,-p-1.0)-q*pow(rij,-q-1.0));
    (pl[i]->pe2)+= go*gsi*F*sw*ex*0.5;
    f -= go*gsi*F*(dsw*ex + sw*dex);

    dgsi = d_gSi(pl[j]->co);
    dfc = d_fc(rij, wbcutsio);
    f -= go*dgsi*dfc*F*sw*ex;
   
    dFr = d_Fc_r(rij, R);
    f -= go*gsi*sw*ex*dFr;
    
    f /= rij;

    (pl[i]->fx)+=f*dxij;
    (pl[i]->fy)+=f*dyij;
    (pl[i]->fz)+=f*dzij;

    if(i < j){
      sys.nfx[0][0] += f * dxij * dxij;
      sys.nfx[0][1] += f * dxij * dyij;
      sys.nfx[0][2] += f * dxij * dzij;
      sys.nfx[1][0] += f * dyij * dxij;
      sys.nfx[1][1] += f * dyij * dyij;
      sys.nfx[1][2] += f * dyij * dzij;
      sys.nfx[2][0] += f * dzij * dxij;
      sys.nfx[2][1] += f * dzij * dyij;
      sys.nfx[2][2] += f * dzij * dzij;
    }
    pl[i]->str[0][0] += 0.5*f*dxij*dxij;
    pl[i]->str[0][1] += 0.5*f*dxij*dyij;
    pl[i]->str[0][2] += 0.5*f*dxij*dzij;
    pl[i]->str[1][0] += 0.5*f*dyij*dxij;
    pl[i]->str[1][1] += 0.5*f*dyij*dyij;
    pl[i]->str[1][2] += 0.5*f*dyij*dzij;
    pl[i]->str[2][0] += 0.5*f*dzij*dxij;
    pl[i]->str[2][1] += 0.5*f*dzij*dyij;
    pl[i]->str[2][2] += 0.5*f*dzij*dzij;

    dgo = d_gO(pl[i]->co);
    neighbor=pl[i]->neighbor;
    for(i2=1;i2<=neighbor;i2++){
      k=pl[i]->j[i2];
      if(pl[k]->atmtype==Si){
        rik=r1(i,k);
        dfc=d_fc(rik, wbcutsio);
        dxik=pl[i]->x-pl[k]->x;
        dyik=pl[i]->y-pl[k]->y;
        dzik=pl[i]->z-pl[k]->z;
        mirror(&dxik,&dyik,&dzik);
        scale_to_real(&dxik,&dyik,&dzik);
        f = -gsi*dgo*dfc*sw*ex/rik;
        (pl[i]->fx)+=f*dxik;
        (pl[i]->fy)+=f*dyik;
        (pl[i]->fz)+=f*dzik;
        if(i < k){
          sys.nfx[0][0] += f * dxik * dxik;
          sys.nfx[0][1] += f * dxik * dyik;
          sys.nfx[0][2] += f * dxik * dzik;
          sys.nfx[1][0] += f * dyik * dxik;
          sys.nfx[1][1] += f * dyik * dyik;
          sys.nfx[1][2] += f * dyik * dzik;
          sys.nfx[2][0] += f * dzik * dxik;
          sys.nfx[2][1] += f * dzik * dyik;
          sys.nfx[2][2] += f * dzik * dzik;
        }
        pl[i]->str[0][0] += 0.5*f*dxik*dxik;
        pl[i]->str[0][1] += 0.5*f*dxik*dyik;
        pl[i]->str[0][2] += 0.5*f*dxik*dzik;
        pl[i]->str[1][0] += 0.5*f*dyik*dxik;
        pl[i]->str[1][1] += 0.5*f*dyik*dyik;
        pl[i]->str[1][2] += 0.5*f*dyik*dzik;
        pl[i]->str[2][0] += 0.5*f*dzik*dxik;
        pl[i]->str[2][1] += 0.5*f*dzik*dyik;
        pl[i]->str[2][2] += 0.5*f*dzik*dzik;
      }
    }

    dFR = d_Fc_R(rij, R);
    dR  = d_Rc(pl[i]->co);
    neighbor=pl[i]->neighbor;
    for(i2=1;i2<=neighbor;i2++){
      k=pl[i]->j[i2];
      if(pl[k]->atmtype==Si){
        rik=r1(i,k);
        dfc=d_fc(rik, wbcutsio);
        dxik=pl[i]->x-pl[k]->x;
        dyik=pl[i]->y-pl[k]->y;
        dzik=pl[i]->z-pl[k]->z;
        mirror(&dxik,&dyik,&dzik);
        scale_to_real(&dxik,&dyik,&dzik);
        f = -gsi*go*dFR*dR*dfc*sw*ex/rik;
        (pl[i]->fx)+=f*dxik;
        (pl[i]->fy)+=f*dyik;
        (pl[i]->fz)+=f*dzik;
        if(i < k){
          sys.nfx[0][0] += f * dxik * dxik;
          sys.nfx[0][1] += f * dxik * dyik;
          sys.nfx[0][2] += f * dxik * dzik;
          sys.nfx[1][0] += f * dyik * dxik;
          sys.nfx[1][1] += f * dyik * dyik;
          sys.nfx[1][2] += f * dyik * dzik;
          sys.nfx[2][0] += f * dzik * dxik;
          sys.nfx[2][1] += f * dzik * dyik;
          sys.nfx[2][2] += f * dzik * dzik;
        }
        pl[i]->str[0][0] += 0.5*f*dxik*dxik;
        pl[i]->str[0][1] += 0.5*f*dxik*dyik;
        pl[i]->str[0][2] += 0.5*f*dxik*dzik;
        pl[i]->str[1][0] += 0.5*f*dyik*dxik;
        pl[i]->str[1][1] += 0.5*f*dyik*dyik;
        pl[i]->str[1][2] += 0.5*f*dyik*dzik;
        pl[i]->str[2][0] += 0.5*f*dzik*dxik;
        pl[i]->str[2][1] += 0.5*f*dzik*dyik;
        pl[i]->str[2][2] += 0.5*f*dzik*dzik;
      }
    }
  }
}

static void force_Si_other(int i)
{
  double A,B,p,q,a,rij,rjk,rik,go,gsi,dgo,dgsi,dfc,f,dxij,dyij,dzij,dxik,dyik,dzik,R,F,dFR,dR,fcore;
  int j,i1,j1,k,neighbori, neighborj;
  double invrija, ex, sw;

  a=wb2[Si][O].a;
  A=wb2[Si][O].A;
  B=wb2[Si][O].B;
  p=wb2[Si][O].p;
  q=wb2[Si][O].q;
  neighbori=pl[i]->neighbor;
  for(i1=1;i1<=neighbori;i1++){
    j=pl[i]->j[i1];
    rij = r1(i,j);
    neighborj = pl[j]->neighbor;
    for(j1=1;j1<=neighborj;j1++){
      k=pl[j]->j[j1];
      if(j<k){
        rjk = r1(j,k);
        rik = r1(i,k);
        if((pl[j]->atmtype == Si)&&(pl[k]->atmtype == O)){
          R = Rc(pl[k]->co);
          if((rjk < a)&&(rjk < R+0.05)){
            dxij=pl[i]->x-pl[j]->x;
            dyij=pl[i]->y-pl[j]->y;
            dzij=pl[i]->z-pl[j]->z;
            mirror(&dxij,&dyij,&dzij);
            scale_to_real(&dxij,&dyij,&dzij);
            dxik=pl[i]->x-pl[k]->x;
            dyik=pl[i]->y-pl[k]->y;
            dzik=pl[i]->z-pl[k]->z;
            mirror(&dxik,&dyik,&dzik);
            scale_to_real(&dxik,&dyik,&dzik);
            gsi = gSi(pl[j]->co);
            F = Fc(rjk,R);
            dFR = d_Fc_R(rjk,R);
            dR = d_Rc(pl[k]->co);
            go = gO(pl[k]->co);
            sw = A*(B*pow(rjk,-p)-pow(rjk,-q))*exp(1.0/(rjk - a));
            if(rik < a){
              dgo = d_gO(pl[k]->co);
              dfc = d_fc(rik, wbcutsio);
              f = -(gsi*F*sw*dgo*dfc+gsi*go*dFR*dR*dfc*sw)/rik;
              (pl[i]->fx)+=f*dxik;
              (pl[i]->fy)+=f*dyik;
              (pl[i]->fz)+=f*dzik;
              if(i < k){
                sys.nfx[0][0] += f*dxik*dxik;
                sys.nfx[0][1] += f*dxik*dyik;
                sys.nfx[0][2] += f*dxik*dzik;
                sys.nfx[1][0] += f*dyik*dxik;
                sys.nfx[1][1] += f*dyik*dyik;
                sys.nfx[1][2] += f*dyik*dzik;
                sys.nfx[2][0] += f*dzik*dxik;
                sys.nfx[2][1] += f*dzik*dyik;
                sys.nfx[2][2] += f*dzik*dzik;
              }
              pl[i]->str[0][0] += 0.5*f*dxik*dxik;
              pl[i]->str[0][1] += 0.5*f*dxik*dyik;
              pl[i]->str[0][2] += 0.5*f*dxik*dzik;
              pl[i]->str[1][0] += 0.5*f*dyik*dxik;
              pl[i]->str[1][1] += 0.5*f*dyik*dyik;
              pl[i]->str[1][2] += 0.5*f*dyik*dzik;
              pl[i]->str[2][0] += 0.5*f*dzik*dxik;
              pl[i]->str[2][1] += 0.5*f*dzik*dyik;
              pl[i]->str[2][2] += 0.5*f*dzik*dzik;
            }
          }
        }
        if((pl[j]->atmtype == O)&&(pl[k]->atmtype == Si)){
          R = Rc(pl[i]->co);
          if((rjk < a)&&(rjk < R+0.05)){
            dxij=pl[i]->x-pl[j]->x;
            dyij=pl[i]->y-pl[j]->y;
            dzij=pl[i]->z-pl[j]->z;
            mirror(&dxij,&dyij,&dzij);
            scale_to_real(&dxij,&dyij,&dzij);
            dxik=pl[i]->x-pl[k]->x;
            dyik=pl[i]->y-pl[k]->y;
            dzik=pl[i]->z-pl[k]->z;
            mirror(&dxik,&dyik,&dzik);
            scale_to_real(&dxik,&dyik,&dzik);
            gsi = gSi(pl[k]->co);
            F = Fc(rjk,R);
            dFR = d_Fc_R(rjk,R);
            dR = d_Rc(pl[k]->co);
            go = gO(pl[j]->co);
            sw = A*(B*pow(rjk,-p)-pow(rjk,-q))*exp(1.0/(rjk - a));
            if(rij < a){
              dgo = d_gO(pl[j]->co);
              dfc = d_fc(rij, wbcutsio);
              f = -(gsi*F*sw*dgo*dfc+gsi*go*dFR*dR*dfc*sw)/rij;
              (pl[i]->fx)+=f*dxij;
              (pl[i]->fy)+=f*dyij;
              (pl[i]->fz)+=f*dzij;
              if(i < j){
                sys.nfx[0][0] += f*dxij*dxij;
                sys.nfx[0][1] += f*dxij*dyij;
                sys.nfx[0][2] += f*dxij*dzij;
                sys.nfx[1][0] += f*dyij*dxij;
                sys.nfx[1][1] += f*dyij*dyij;
                sys.nfx[1][2] += f*dyij*dzij;
                sys.nfx[2][0] += f*dzij*dxij;
                sys.nfx[2][1] += f*dzij*dyij;
                sys.nfx[2][2] += f*dzij*dzij;
              }
              pl[i]->str[0][0] += 0.5*f*dxij*dxij;
              pl[i]->str[0][1] += 0.5*f*dxij*dyij;
              pl[i]->str[0][2] += 0.5*f*dxij*dzij;
              pl[i]->str[1][0] += 0.5*f*dyij*dxij;
              pl[i]->str[1][1] += 0.5*f*dyij*dyij;
              pl[i]->str[1][2] += 0.5*f*dyij*dzij;
              pl[i]->str[2][0] += 0.5*f*dzij*dxij;
              pl[i]->str[2][1] += 0.5*f*dzij*dyij;
              pl[i]->str[2][2] += 0.5*f*dzij*dzij;
            }
          }
        }
      }
    }
  }
}

static void force_O_other(int i)
{
  double A,B,p,q,a,rij,rjk,rik,go,gsi,dgo,dgsi,dfc,f,dxij,dyij,dzij,dxik,dyik,dzik,R,F,dF,dR,fcore;
  int j,i1,j1,k,neighbori, neighborj;
  double invrija, ex, sw;

  a=wb2[Si][O].a;
  A=wb2[Si][O].A;
  B=wb2[Si][O].B;
  p=wb2[Si][O].p;
  q=wb2[Si][O].q;
  neighbori=pl[i]->neighbor;
  for(i1=1;i1<=neighbori;i1++){
    j=pl[i]->j[i1];
    rij = r1(i,j);
    neighborj = pl[j]->neighbor;
    for(j1=1;j1<=neighborj;j1++){
      k=pl[j]->j[j1];
      if(j<k){
        rjk = r1(j,k);
        rik = r1(i,k);
        if((pl[j]->atmtype == O)&&(pl[k]->atmtype == Si)){
          R = Rc(pl[j]->co);
          if((rjk < a)&&(rjk < R+0.05)){
            dxij=pl[i]->x-pl[j]->x;
            dyij=pl[i]->y-pl[j]->y;
            dzij=pl[i]->z-pl[j]->z;
            mirror(&dxij,&dyij,&dzij);
            scale_to_real(&dxij,&dyij,&dzij);
            dxik=pl[i]->x-pl[k]->x;
            dyik=pl[i]->y-pl[k]->y;
            dzik=pl[i]->z-pl[k]->z;
            mirror(&dxik,&dyik,&dzik);
            scale_to_real(&dxik,&dyik,&dzik);
            gsi = gSi(pl[k]->co);
            F = Fc(rjk,R);
            go = gO(pl[j]->co);
            sw = A*(B*pow(rjk,-p)-pow(rjk,-q))*exp(1.0/(rjk - a));
            if(rik < a){
              dgsi = d_gSi(pl[k]->co);
              dfc = d_fc(rik, wbcutsio);
              f = -go*dgsi*F*sw*dfc/rik;
              (pl[i]->fx)+=f*dxik;
              (pl[i]->fy)+=f*dyik;
              (pl[i]->fz)+=f*dzik;
              if(i < k){
                sys.nfx[0][0] += f*dxik*dxik;
                sys.nfx[0][1] += f*dxik*dyik;
                sys.nfx[0][2] += f*dxik*dzik;
                sys.nfx[1][0] += f*dyik*dxik;
                sys.nfx[1][1] += f*dyik*dyik;
                sys.nfx[1][2] += f*dyik*dzik;
                sys.nfx[2][0] += f*dzik*dxik;
                sys.nfx[2][1] += f*dzik*dyik;
                sys.nfx[2][2] += f*dzik*dzik;
              }
              pl[i]->str[0][0] += 0.5*f*dxik*dxik;
              pl[i]->str[0][1] += 0.5*f*dxik*dyik;
              pl[i]->str[0][2] += 0.5*f*dxik*dzik;
              pl[i]->str[1][0] += 0.5*f*dyik*dxik;
              pl[i]->str[1][1] += 0.5*f*dyik*dyik;
              pl[i]->str[1][2] += 0.5*f*dyik*dzik;
              pl[i]->str[2][0] += 0.5*f*dzik*dxik;
              pl[i]->str[2][1] += 0.5*f*dzik*dyik;
              pl[i]->str[2][2] += 0.5*f*dzik*dzik;
            }
          }
        }
        if((pl[j]->atmtype == Si)&&(pl[k]->atmtype == O)){
          R = Rc(pl[k]->co);
          if((rjk < a)&&(rjk < R+0.05)){
            dxij=pl[i]->x-pl[j]->x;
            dyij=pl[i]->y-pl[j]->y;
            dzij=pl[i]->z-pl[j]->z;
            mirror(&dxij,&dyij,&dzij);
            scale_to_real(&dxij,&dyij,&dzij);
            dxik=pl[i]->x-pl[k]->x;
            dyik=pl[i]->y-pl[k]->y;
            dzik=pl[i]->z-pl[k]->z;
            mirror(&dxik,&dyik,&dzik);
            scale_to_real(&dxik,&dyik,&dzik);
            gsi = gSi(pl[j]->co);
            F = Fc(rjk,R);
            go = gO(pl[k]->co);
            sw = A*(B*pow(rjk,-p)-pow(rjk,-q))*exp(1.0/(rjk - a));
            if(rij < a){
              dgsi = d_gSi(pl[j]->co);
              dfc = d_fc(rij, wbcutsio);
              f = -go*dgsi*F*sw*dfc/rij;
              (pl[i]->fx)+=f*dxij;
              (pl[i]->fy)+=f*dyij;
              (pl[i]->fz)+=f*dzij;
              if(i < j){
                sys.nfx[0][0] += f*dxij*dxij;
                sys.nfx[0][1] += f*dxij*dyij;
                sys.nfx[0][2] += f*dxij*dzij;
                sys.nfx[1][0] += f*dyij*dxij;
                sys.nfx[1][1] += f*dyij*dyij;
                sys.nfx[1][2] += f*dyij*dzij;
                sys.nfx[2][0] += f*dzij*dxij;
                sys.nfx[2][1] += f*dzij*dyij;
                sys.nfx[2][2] += f*dzij*dzij;
              }
              pl[i]->str[0][0] += 0.5*f*dxij*dxij;
              pl[i]->str[0][1] += 0.5*f*dxij*dyij;
              pl[i]->str[0][2] += 0.5*f*dxij*dzij;
              pl[i]->str[1][0] += 0.5*f*dyij*dxij;
              pl[i]->str[1][1] += 0.5*f*dyij*dyij;
              pl[i]->str[1][2] += 0.5*f*dyij*dzij;
              pl[i]->str[2][0] += 0.5*f*dzij*dxij;
              pl[i]->str[2][1] += 0.5*f*dzij*dyij;
              pl[i]->str[2][2] += 0.5*f*dzij*dzij;
            }
          }
        }
      }
    }
  }
}


static double fc(double r, WatanabeCutoff cutoff_param)
  /* r : interatomic distance */
{
  double value,disl,disr;

  disl = cutoff_param.R - cutoff_param.D;
  disr = cutoff_param.R + cutoff_param.D;

  if(r<=disl){
    value=1.0;
  }
  else if(r>=disr){
    value=0.0;
  }
  else{
    value=0.5/PI*(sin(PI*(r-disl)/cutoff_param.D))
          -0.5*(r-disl)/cutoff_param.D+1.0;
  }
  return value;
}

static double d_fc(double r, WatanabeCutoff cutoff_param)
  /* r : interatomic distance */
{
  double value;

  if((r>(cutoff_param.R-cutoff_param.D))&&(r<(cutoff_param.R+cutoff_param.D))){
    value=0.5/cutoff_param.D
          *cos(PI/cutoff_param.D*(r-cutoff_param.R+cutoff_param.D))
          -0.5/cutoff_param.D;
/*
    printf("trangent:dfc=%f\n",value);
*/
  }
  else{
    value=0.0;
  }
  return value;
}

static double bsSi(double co)
{
  double a = 4.0, b = 0.5;

  if(co <= 4.0){
    return 1.0;
  } else {
    return (a/co - 1.0)*exp(b/(4.0-co))+1.0;
  }
}

static double d_bsSi(double co)
{
  double a = 4.0, b = 0.5;

  if(co <= 4.0){
    return 0.0;
  } else {
    return exp(b/(4.0-co))*(-1.0*b*co*co - a*(co*co-co*(8.0+b)+16.0))/((co-4.0)*(co-4.0)*co*co);
  }
}

static double gO(double co)
{
  double tmp;

  tmp = co + wbsat.ne;
  return (wbsat.nc*exp(wbsat.nd*tmp*tmp)/(exp((wbsat.na-co)/wbsat.nb)+1.0));
}

static double d_gO(double co)
{
  double tmp[7];

  tmp[0] = exp(wbsat.na/wbsat.nb);
  tmp[1] = exp(co/wbsat.nb);
  tmp[2] = co + wbsat.ne;
  tmp[3] = tmp[0] + tmp[1];
  tmp[4] = wbsat.nc*tmp[1]*exp(wbsat.nd*tmp[2]*tmp[2]);
  tmp[5] = tmp[0]+2.0*wbsat.nb*wbsat.nd*tmp[3]*tmp[2];
  tmp[6] = wbsat.nb*tmp[3]*tmp[3];
  return (tmp[4]*tmp[5]/tmp[6]);
}

static double gSi(double co)
{
  if(co < 4.0){
    return (wbsat.sia*sqrt(co+wbsat.sib)-wbsat.sid)
            *exp(wbsat.sic/(co-4.0))+wbsat.sid;
  } else {
    return wbsat.sid;
  }
}

static double d_gSi(double co)
{
  double tmp[2];

  if(co < 4.0){
    tmp[0] = sqrt(co + wbsat.sib);
    tmp[1] = exp(wbsat.sic/(co-4.0));
    return 0.5*tmp[1]*(wbsat.sia/tmp[0]-2.0*wbsat.sic
           *(wbsat.sia*tmp[0]-wbsat.sid)/(co-4.0)/(co-4.0));
  } else {
    return 0.0;
  }
}

static double lambda(double z, Watanabe3Body * parameter){
  return (parameter->lam)*(1.0+(parameter->g[0])
         *exp(-(parameter->g[1])*pow(z-(parameter->g[2]),2.0)));
}

static double d_lambda(double z, Watanabe3Body *  parameter){
  return -2.0*(parameter->lam)*(parameter->g[0])*(parameter->g[1])
         *(z-(parameter->g[2]))*exp(-(parameter->g[1])
         *pow(z-(parameter->g[2]),2.0));
}

static double Rc(double z)
{
  double value;

  return 1.35;
  return 0.8;

  if(z<=2.0){
    value=1.25;
  }
  else if(z>=3.0){
    value=0.8;
  }
  else{
    value=0.225/PI*(sin(PI*(z-2.0)/0.5))
          -0.225*(z-2.0)/0.5 + 1.25;
  }
  return value;
}

static double d_Rc(double z)
{
  double value;

  return 0.0;

  if((z>2.0)&&(z<3.0)){
    value=0.225/0.5
          *cos(PI/0.5*(z-2.0))
          -0.225/0.5;
/*
    printf("gradient %f\n",value);
*/
  }
  else{
    value=0.0;
  }
  return value;
}

static double Fc(double r, double R)
{
  double value,disl,disr,D;

  return 1.0;

  D = 0.05;

  disl = R - D;
  disr = R + D;

  if(r<=disl){
    value=1.0;
  }
  else if(r>=disr){
    value=0.0;
  }
  else{
    value=0.5/PI*(sin(PI*(r-disl)/D))
          -0.5*(r-disl)/D+1.0;
  }
  return value;
}

static double d_Fc_r(double r, double R)
{
  double value,D = 0.05;

  return 0.0;

  if((r>(R-D))&&(r<(R+D))){
    value=0.5/D*cos(PI/D*(r-R+D))-0.5/D;
  }
  else{
    value=0.0;
  }
  return value;
}

static double d_Fc_R(double r, double R)
{
  double value,D=0.05;

  return 0.0;

  if((r>(R-D))&&(r<(R+D))){
    value= -0.5/D*cos(PI/D*(r-R+D))+0.5/D;
  }
  else{
    value=0.0;
  }
  return value;
}

void set_SIMPLECOULOMB_parameter(void)
{
  cutoff[O][O] = 1.6;
  cutoff[Si][O] = 1.8/sig;
  cutoff[O][Si] = 1.8/sig;
  cutoff[Si][Si] = 3.5/sig;
}

void calc_SIMPLECOULOMB_force(void)
{
  int i,j,i1,neighbor;
  double total_charge = 0.0;
  double A,B,p,q,a;
  double rij,dx,dy,dz,f;
  double tmp[7];

  for(i=0;i<N;i++){
    pl[i]->co=0.0;  /* initialize */
    pl[i]->self_co = 0.0;
  }
  for(i=0;i<N;i++){
    if(pl[i]->atmtype==O){
      neighbor=pl[i]->neighbor;
      for(i1=1;i1<=neighbor;i1++){
        j=pl[i]->j[i1];
        if(pl[j]->atmtype==Si){
          if(r1(i,j) < 1.8/sig){
            pl[i]->co += 1.0;
            pl[j]->co += 1.0;
          }
        }
      }
    }
    if(pl[i]->atmtype==Si){
      neighbor=pl[i]->neighbor;
      for(i1=1;i1<=neighbor;i1++){
        j=pl[i]->j[i1];
        if(pl[j]->atmtype==Si){
          if(r1(i,j)<3.4/sig){
            pl[i]->self_co += 1.0;
          }
        }
      }
    }
  }
  for(i=0;i<N;i++){
    if(pl[i]->atmtype == Si){
      pl[i]->charge = +0.4*(double)(pl[i]->co);
    } else if (pl[i]->atmtype == O){
      pl[i]->charge = -0.4*(double)(pl[i]->co);
    } else {
      pl[i]->charge = 0.0;
    }
/*
    printf("%d %d %f %f\n",pl[i]->atmtype, i, pl[i]->co, pl[i]->self_co);
*/
  }

  for(i=0;i<N;i++){
    if(pl[i]->atmtype == Si){
      pl[i]->ion = (1.18946 * pl[i]->charge + 3.239735 * pl[i]->charge * pl[i]->charge);
    }
    if(pl[i]->atmtype == O){
      pl[i]->ion = exp(0.525/pl[i]->charge)*(3.27 * pl[i]->charge+2.71 * pl[i]->charge * pl[i]->charge);
    }
  }

  for(i=0;i<N;i++){
    if(pl[i]->atmtype == Si){
      pl[i]->pe2 +=  0.0*(-3.9744*pl[i]->co - (pl[i]->self_co - pl[i]->co)*0.5*2.48853971)*eV2eps;
/** Hamman 
      if(fabs(pl[i]->co - 1.0) < 0.01){
        pl[i]->pe2 += 0.47  * eV2eps;
      }
      if(fabs(pl[i]->co - 2.0) < 0.01){
        pl[i]->pe2 += 0.51 * eV2eps;
      }
      if(fabs(pl[i]->co - 3.0) < 0.01){
        pl[i]->pe2 += 0.24 * eV2eps;
      }
**/
/*
      if(fabs(pl[i]->co - 1.0) < 0.01){
        pl[i]->pe2 += 0.611275 * eV2eps;
      }
      if(fabs(pl[i]->co - 2.0) < 0.01){
        pl[i]->pe2 += 0.616009 * eV2eps;
      }
      if(fabs(pl[i]->co - 3.0) < 0.01){
        pl[i]->pe2 += 0.290535 * eV2eps;
      }
*/
    }
  }


/*
  for(i=0;i<N;i++){
    if(pl[i]->atmtype == Si){
      pl[i]->charge = +1.6;
       total_charge += 1.6;
    } else if (pl[i]->atmtype == O){
      pl[i]->charge = -0.8;
       total_charge -= 0.8;
    } else {
      pl[i]->charge = 0.0;
    }
  }
  if(fabs(total_charge) > 0.00000001){
    fprintf(stderr,"Error : Total charge must be zero.\n");
    exit(1);
  }
*/

  execute_pme();

  for(i=0;i<N;i++){
    pl[i]->fx = 0.0;
    pl[i]->fy = 0.0;
    pl[i]->fz = 0.0;
  }

/* CORE
  a = 1.30;
  A = 21.0;
  B = 0.038;
  p = 5.3;
  q = -1.1;
  for(i=0;i<N;i++){
    neighbor = pl[i]->neighbor;
    for(i1=1;i1<=neighbor;i1++){
      j = pl[i]->j[i1];
      rij = r1(i,j);
      if(rij < a){
        dx=pl[i]->x-pl[j]->x;
        dy=pl[i]->y-pl[j]->y;
        dz=pl[i]->z-pl[j]->z;
        mirror(&dx,&dy,&dz);
        scale_to_real(&dx,&dy,&dz);
        tmp[0]=1.0/(rij-a);
        tmp[1]=exp(tmp[0]);
        tmp[2]=A*(B*pow(rij,-p)-pow(rij,-q));
        tmp[3]=tmp[0]*tmp[0];
        tmp[4]=A*(B*p*pow(rij,-p-1.0)-q*pow(rij,-q-1.0));
        tmp[5]=tmp[1]*tmp[2];
        tmp[6]= tmp[5]*0.5;
        (pl[i]->pe2)+=tmp[6];
        f=(tmp[4]+tmp[2]*tmp[3])*tmp[1]/rij;
        (pl[i]->fx)+=f*dx;
        (pl[i]->fy)+=f*dy;
        (pl[i]->fz)+=f*dz;
      }
    }
  }
*/
}

} // end of namespace WASEDA_MD_LABO

/*** end of watanabe.c ***/
