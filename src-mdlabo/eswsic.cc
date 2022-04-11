/******************************************************************************
                                  eswsic.cc

    --- Extended Stillinger-Weber Potential for Si, C, O Mixed Systems ---

      Copyright (c) N. Aoki, T. Zushi, T.Watanabe 2012  All Rights Reserved
******************************************************************************/
#include "mddef.h"

namespace WASEDA_MD_LABO {

typedef struct {
  double A;
  double B;
  double p;
  double q;
  double a;
  double d;
} ESWSIC2Body;

typedef struct {
  double lam1, lam2;
  double gam1ij, gam1ik, gam2ij, gam2ik;
  double a1ij, a1ik, a2ij, a2ik;
  double cos01, cos02;
  double alpha1, alpha2;
} ESWSIC3Body;

static ESWSIC2Body * get_ESWSIC_2body(AtomIonType, AtomIonType);
static ESWSIC3Body * get_ESWSIC_3body(AtomIonType, AtomIonType, AtomIonType);




///////////////////////////////////////////////////////////////////
static double fc_SiO(double);
static double d_fc_SiO(double);
static double gO_SiO(double);
static double d_gO_SiO(double);
static double cutR_SiO, cutD_SiO, M1_SiO, M2_SiO, M3_SiO, M4_SiO, M5_SiO;

static double fc_CO(double);
static double d_fc_CO(double);
static double gO_CO(double);
static double d_gO_CO(double);
static double cutR_CO, cutD_CO, M1_CO, M2_CO, M3_CO, M4_CO, M5_CO;
///////////////////////////////////////////////////////////////////




void set_ESWSIC_parameter(void)   // to be modified
{
  printf("set_ESWSIC_parameter\n");

  double sigma = 2.0951;

  cutR_SiO = 1.30*sigma;//Dr [A]
  cutD_SiO = 0.1*sigma;//Dr
  M1_SiO = 0.097;//Dr
  M2_SiO = 1.6;//Dr
  M3_SiO = 0.3654;//Dr
  M4_SiO = 0.1344;//Dr
  M5_SiO = 6.4176;//Dr

  cutR_CO = 1.7;//to be modified [A]
  cutD_CO = 0.2;//to be modified
  M1_CO = 1.24246;//self
  M2_CO = -36.1839;//self
  M3_CO = 0.0295862;//self
  M4_CO = -0.287896;//self
  M5_CO = 1.13163;//self
}

double set_ESWSIC_scut(void)
{
  double cut, tmpcut;
  int i,j;

  cut = 0.0;
  for(i=0;i<MAXELEMENT;i++){
    for(j=i;j<MAXELEMENT;j++){
      tmpcut = get_ESWSIC_scutoff(AtomIonType(i), AtomIonType(j));
      if(cut < tmpcut) cut = tmpcut;
    }
  }

  return cut;
}

double get_ESWSIC_scutoff(AtomIonType a1, AtomIonType a2)
{
  double cut = 0.0;
  int i;
  ESWSIC2Body *param2;
  ESWSIC3Body *param3;

  double sigma = 2.0951;

  if(a1 == Si && a2 == Si) return 1.80*sigma;//Dr [A]
  if(a1 == Si && a2 == O ) return 1.40*sigma;//Dr
  if(a1 == O  && a2 == Si) return 1.40*sigma;//Dr
  if(a1 == O  && a2 == O ) return 1.25*sigma;//Dr
  if(a1 == C  && a2 == C) return 3.0;  	// to be modified completed [A]
  if(a1 == C  && a2 == O) return 2.4;  // to be modified completed
  if(a1 == O  && a2 == C) return 2.4;  // to be modified completed
  if(a1 == C  && a2 == Si) return 3.9;  // to be modified completed
  if(a1 == Si  && a2 == C) return 3.9;  // to be modified completed

/*

  if(a1 == Omol && a2 == Omol) return 1.25*sigma;//md-4.19
  if(a1 == Omol && a2 == O ) return 1.43*sigma;//md-4.19
  if(a1 == O  && a2 == Omol) return 1.43*sigma;//md-4.19
  if(a1 == Si && a2 == Omol) return 1.81*sigma;//md-4.19
  if(a1 == Omol && a2 == Si) return 1.81*sigma;//md-4.19
  if(a1 == C  && a2 == Omol) return 3.30;  // to be modified
  if(a1 == Omol  && a2 == C) return 3.30;  // to be modified

*/

  param2 = get_ESWSIC_2body(a1, a2);
  if(cut < param2->a) cut = param2->a;
  for(i=0;i<MAXELEMENT;i++){
    param3 = get_ESWSIC_3body(a1, a2, (AtomIonType)i);
    if(cut < param3->a1ij) cut = param3->a1ij;
    if(cut < param3->a2ik) cut = param3->a2ik;
    param3 = get_ESWSIC_3body(a2, a1, (AtomIonType)i);
    if(cut < param3->a1ij) cut = param3->a1ij;
    if(cut < param3->a2ik) cut = param3->a2ik;
  }

  return cut;
}

static ESWSIC2Body * get_ESWSIC_2body(AtomIonType ai, AtomIonType aj)
{
  static ESWSIC2Body prm;

  prm.A = 0.0;
  prm.B = 0.0;
  prm.p = 0.0;
  prm.q = 0.0;
  prm.a = 0.1;
  prm.d = 1.0;

  double sigma = 2.0951;//dimension change [no-dimension -> A]
  double epsilon = 50.0;//[kcal/mol]
  double dim_change=0.04337;//eV/(kcal/mol)
  double norm1;//dimension change (no-dimension -> eV)

  if(ai == Si && aj == Si){//Dr
    prm.p = 4.0;//[no-dimension]
    prm.q = 0.0;//[no-dimension]
    norm1=epsilon*dim_change*pow(sigma,prm.q);
    prm.A = 7.049556277*norm1;//[eV*A^(q)]
    prm.B = 0.6022245584*pow(sigma,prm.p-prm.q);//[A^(p-q)] (old-B*sigma^(p-q))
    prm.a = 1.80*sigma;//[A]
    prm.d = 1.0*sigma;//[A]
    return &prm;
  }

  if((ai == Si && aj == O) || (ai == O && aj == Si)){//Dr
    prm.p = 2.58759;
    prm.q = 2.39370;
    norm1=epsilon*dim_change*pow(sigma,prm.q);
    prm.A = 115.364065913*norm1;
    prm.B = 0.9094442793*pow(sigma,prm.p-prm.q);
    prm.a = 1.40*sigma;
    prm.d = 1.0*sigma;
    return &prm;
  }

  if(ai == O && aj == O){//Dr
    prm.p = 0.0;
    prm.q = 2.24432;
    norm1=epsilon*dim_change*pow(sigma,prm.q);
    prm.A = -12.292427744*norm1;
    prm.B = 0.0*pow(sigma,prm.p-prm.q);
    prm.a = 1.25*sigma;
    prm.d = 1.0*sigma;
    return &prm;
  }

  if(ai == C && aj == C){	// to be modified completed
    prm.p = 3.5;//[no-dimension]
    prm.q = 2.0;//[no-dimension]
    prm.A = 44.28928;//[eV*A^(q)] (old-A * norm * sigma^(q))
    prm.B = 1.42489;//[A^(p-q)] (old-B * sigma^(p-q))
    prm.a = 3.0;//[A]
    prm.d = 0.87096;//[A]
    return &prm;
  }

  if((ai == C && aj == O) || (ai == O && aj == C)){	// to be modified completed
    prm.p = 4.0;
    prm.q = 2.0;
    prm.A = 31.6872;
    prm.B = 1.23825;
    prm.a = 2.4;
    prm.d = 0.704113;
    return &prm;
  }

  if((ai == C && aj == Si) || (ai == Si && aj == C)){	// to be modified completed
    prm.p = 3.5;
    prm.q = 1.8;
    prm.A = 70.1942;
    prm.B = 1.85335;
    prm.a = 3.9;
    prm.d = 1.63943;
    return &prm;
  }

/*

  if(ai == Omol && aj == Omol){//md-4.19
    prm.p = 3.6;
    prm.q = 2.63;
    norm1=epsilon*dim_change*pow(sigma,prm.q);
    prm.A = 11.57*norm1;
    prm.B = 0.4598*pow(sigma,prm.p-prm.q);
    prm.a = 1.25*sigma;
    prm.d = 1.0*sigma;
    return &prm;
  }

  if((ai == Omol && aj == O) || (ai == O && aj == Omol)){//md-4.19
    prm.p = 0.0;
    prm.q = 3.0;
    norm1=epsilon*dim_change*pow(sigma,prm.q);
    prm.A = -38.0*norm1;
    prm.B = 0.0*pow(sigma,prm.p-prm.q);
    prm.a = 1.43*sigma;
    prm.d = 1.0*sigma;
    return &prm;
  }

  if((ai == Omol && aj == Si) || (ai == Si && aj == Omol)){//md-4.19
    prm.p = 0.0;
    prm.q = 3.6;
    norm1=epsilon*dim_change*pow(sigma,prm.q);
    prm.A = -46.0*norm1;
    prm.B = 0.0*pow(sigma,prm.p-prm.q);
    prm.a = 1.81*sigma;
    prm.d = 1.0*sigma;
    return &prm;
  }

  if((ai == Omol && aj == C) || (ai == C && aj == Omol)){//to be modified
    prm.p = 0.0;
    prm.q = 0.0;
    prm.A = 0.0;
    prm.B = 0.0;
    prm.a = 0.1;
    prm.d = 1.0;
    return &prm;
  }

*/

  return &prm;
}

static ESWSIC3Body * get_ESWSIC_3body(AtomIonType aj,
                                      AtomIonType ai,
                                      AtomIonType ak)
{
  static ESWSIC3Body prm;

  prm.lam1 = 0.0;
  prm.gam1ij = 1.0;
  prm.gam1ik = 1.0;
  prm.a1ij = 0.1;
  prm.a1ik = 0.1;
  prm.cos01 = 0.0;
  prm.alpha1 = 0.0;

  prm.lam2 = 0.0;
  prm.gam2ij = 1.0;
  prm.gam2ik = 1.0;
  prm.a2ij = 0.1;
  prm.a2ik = 0.1;
  prm.cos02 = 0.0;
  prm.alpha2 = 0.0;

  double sigma = 2.0951;
  double epsilon = 50.0;//[kcal/mol]
  double dim_change=0.04337;//eV/(kcal/mol)
  double norm2 = epsilon*dim_change;//dimension change (epsilon -> eV)

//{Si}
  if(aj == Si && ai == Si && ak == Si){//Dr

    prm.lam1 = 16.404*norm2;//[eV]
    prm.gam1ij = 1.0473*sigma;//[A]
    prm.gam1ik = 1.0473*sigma;//[A]
    prm.a1ij = 1.80*sigma;//[A]
    prm.a1ik = 1.80*sigma;//[A]
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0;

    return &prm;
  }

//{Si,O}
  if(aj == Si && ai == O && ak == Si){//Dr

    prm.lam1 = 2.9572*norm2;
    prm.gam1ij = 0.71773*sigma;
    prm.gam1ik = 0.71773*sigma;
    prm.a1ij = 1.4*sigma;
    prm.a1ik = 1.4*sigma;
    prm.cos01 = -0.6155238;
    prm.alpha1 = 0.0;

    prm.lam2 = 2.9572*norm2;
    prm.gam2ij = 0.71773*sigma;
    prm.gam2ik = 0.71773*sigma;
    prm.a2ij = 1.4*sigma;
    prm.a2ik = 1.4*sigma;
    prm.cos02 = -0.6155238;
    prm.alpha2 = 0.0;

    return &prm;
  }

  if(aj == O && ai == Si && ak == O){//Dr
    prm.lam1 = 3.1892*norm2;
    prm.gam1ij = 0.3220*sigma;
    prm.gam1ik = 0.3220*sigma;
    prm.a1ij = 1.65*sigma;
    prm.a1ik = 1.65*sigma;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0;

    prm.lam2 = 3.1892*norm2;
    prm.gam2ij = 0.3220*sigma;
    prm.gam2ik = 0.3220*sigma;
    prm.a2ij = 1.65*sigma;
    prm.a2ik = 1.65*sigma;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.0;

    return &prm;
  }

  if(aj == Si && ai == Si && ak == O){//Dr
    prm.lam1 = 10.667*norm2;
    prm.gam1ij = 1.93973*sigma;
    prm.gam1ik = 0.25*sigma;
    prm.a1ij = 1.9*sigma;
    prm.a1ik = 1.4*sigma;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0;

    prm.lam2 = 10.667*norm2;
    prm.gam2ij = 1.93973*sigma;
    prm.gam2ik = 0.25*sigma;
    prm.a2ij = 1.9*sigma;
    prm.a2ik = 1.4*sigma;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.0;

    return &prm;
  }

  if(aj == O && ai == Si && ak == Si){//Dr
    prm.lam1 = 10.667*norm2;
    prm.gam1ij = 1.93973*sigma;
    prm.gam1ik = 0.25*sigma;
    prm.a1ij = 1.9*sigma;
    prm.a1ik = 1.4*sigma;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0;

    prm.lam2 = 10.667*norm2;
    prm.gam2ij = 1.93973*sigma;
    prm.gam2ik = 0.25*sigma;
    prm.a2ij = 1.9*sigma;
    prm.a2ik = 1.4*sigma;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.0;

    return &prm;
  }

//{C}
  if(aj == C && ai == C && ak == C){	// to be modified completed
    prm.lam1 = 579.691;//[eV]
    prm.gam1ij = 4.13236;//[A]
    prm.gam1ik = 4.13236;//[A]
    prm.a1ij = 3.35398;//[A]
    prm.a1ik = 3.35398;//[A]
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0001;
    return &prm;
  }

//{C,O}
  if(aj == C && ai == O && ak == C){	// to be modified completed
    prm.lam1 = 109.567;
    prm.gam1ij = 2.16786;
    prm.gam1ik = 2.16786;
    prm.a1ij = 2.54447;
    prm.a1ik = 2.54447;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0001;
    prm.lam2 = 109.567;
    prm.gam2ij = 2.16786;
    prm.gam2ik = 2.16786;
    prm.a2ij = 2.54447;
    prm.a2ik = 2.54447;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.0001;
    return &prm;
  }

  if(aj == O && ai == C && ak == O){	// to be modified
    prm.lam1 = 18.7243;//18.7243
    prm.gam1ij = 1.17075;
    prm.gam1ik = 1.17075;
    prm.a1ij = 2.72792;
    prm.a1ik = 2.72792;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0001;
    prm.lam2 = 18.7243;//18.7243
    prm.gam2ij = 1.17075;
    prm.gam2ik = 1.17075;
    prm.a2ij = 2.72792;
    prm.a2ik = 2.72792;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.0001;
    return &prm;
  }

  if(aj == C && ai == C && ak == O){	// to be modified completed
    prm.lam1 = 339.409;
    prm.gam1ij = 35.8556;
    prm.gam1ik = -4.36577;
    prm.a1ij = 4.96774;
    prm.a1ik = 2.05934;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0001;
    prm.lam2 = 339.409;
    prm.gam2ij = 35.8556;
    prm.gam2ik = -4.36577;
    prm.a2ij = 4.96774;
    prm.a2ik = 2.05934;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.0001;
    return &prm;
  }

  if(aj == O && ai == C && ak == C){	// to be modified completed
    prm.lam1 = 339.409;
    prm.gam1ij = -4.36577;
    prm.gam1ik = 35.8556;
    prm.a1ij = 2.05934;
    prm.a1ik = 4.96774;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0001;
    prm.lam2 = 339.409;
    prm.gam2ij = -4.36577;
    prm.gam2ik = 35.8556;
    prm.a2ij = 2.05934;
    prm.a2ik = 4.96774;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.0001;
    return &prm;
  }

//{Si,C}
  if(aj == Si && ai == C && ak == Si){	// to be modified completed
    prm.lam1 = 112.794;
    prm.gam1ij = 3.48977;
    prm.gam1ik = 3.48977;
    prm.a1ij = 3.9;
    prm.a1ik = 3.9;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0001;
    prm.lam2 = 112.794;
    prm.gam2ij = 3.48977;
    prm.gam2ik = 3.48977;
    prm.a2ij = 3.9;
    prm.a2ik = 3.9;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.0001;
    return &prm;
  }

  if(aj == C && ai == Si && ak == C){	// to be modified completed
    prm.lam1 = 46.6202;
    prm.gam1ij = 2.72235;
    prm.gam1ik = 2.72235;
    prm.a1ij = 3.72789;
    prm.a1ik = 3.72789;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0001;
    prm.lam2 = 46.6202;
    prm.gam2ij = 2.72235;
    prm.gam2ik = 2.72235;
    prm.a2ij = 3.72789;
    prm.a2ik = 3.72789;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.0001;
    return &prm;
  }

  if(aj == Si && ai == Si && ak == C){	// to be modified completed
    prm.lam1 = 34.2004;
    prm.gam1ij = 4.06712;
    prm.gam1ik = 3.1453;
    prm.a1ij = 3.9258;
    prm.a1ik = 3.9;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0001;
    prm.lam2 = 34.2004;
    prm.gam2ij = 4.06712;
    prm.gam2ik = 3.1453;
    prm.a2ij = 3.9258;
    prm.a2ik = 3.9;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.0001;
    return &prm;
  }

  if(aj == C && ai == Si && ak == Si){	// to be modified completed
    prm.lam1 = 34.2004;
    prm.gam1ij = 3.1453;
    prm.gam1ik = 4.06712;
    prm.a1ij = 3.9;
    prm.a1ik = 3.9258;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0001;
    prm.lam2 = 34.2004;
    prm.gam2ij = 3.1453;
    prm.gam2ik = 4.06712;
    prm.a2ij = 3.9;
    prm.a2ik = 3.9258;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.0001;
    return &prm;
  }

  if(aj == C && ai == C && ak == Si){	// to be modified completed
    prm.lam1 = 19.3261;
    prm.gam1ij = 1.62965;
    prm.gam1ik = 2.5356;
    prm.a1ij = 3.11521;
    prm.a1ik = 4.47205;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0001;
    prm.lam2 = 19.3261;
    prm.gam2ij = 1.62965;
    prm.gam2ik = 2.5356;
    prm.a2ij = 3.11521;
    prm.a2ik = 4.47205;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.0001;
    return &prm;
  }

  if(aj == Si && ai == C && ak == C){	// to be modified completed
    prm.lam1 = 19.3261;
    prm.gam1ij = 2.5356;
    prm.gam1ik = 1.62965;
    prm.a1ij = 4.47205;
    prm.a1ik = 3.11521;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0001;
    prm.lam2 = 19.3261;
    prm.gam2ij = 2.5356;
    prm.gam2ik = 1.62965;
    prm.a2ij = 4.47205;
    prm.a2ik = 3.11521;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.0001;
    return &prm;
  }

//{Si,C,O}
  if(aj == Si && ai == C && ak == O){	// to be modified completed
    prm.lam1 = 207.581;
    prm.gam1ij = 11.0585;
    prm.gam1ik = -1.73465;
    prm.a1ij = 3.9;
    prm.a1ik = 2.4;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0001;
    prm.lam2 = 207.581;
    prm.gam2ij = 11.0585;
    prm.gam2ik = -1.73465;
    prm.a2ij = 3.9;
    prm.a2ik = 2.4;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.0001;
    return &prm;
  }

  if(aj == Si && ai == O && ak == C){	// to be modified completed
    prm.lam1 = 51.945;
    prm.gam1ij = 2.63341;
    prm.gam1ik = 1.898;
    prm.a1ij = 2.75;
    prm.a1ik = 2.4;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0001;
    prm.lam2 = 51.945;
    prm.gam2ij = 2.63341;
    prm.gam2ik = 1.898;
    prm.a2ij = 2.75;
    prm.a2ik = 2.4;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.0001;
    return &prm;
  }

  if(aj == C && ai == Si && ak == O){	// to be modified completed
    prm.lam1 = 37.0115;
    prm.gam1ij = 1.69221;
    prm.gam1ik = 2.33945;
    prm.a1ij = 3.9;
    prm.a1ik = 2.75;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0001;
    prm.lam2 = 37.0115;
    prm.gam2ij = 1.69221;
    prm.gam2ik = 2.33945;
    prm.a2ij = 3.9;
    prm.a2ik = 2.75;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.0001;
    return &prm;
  }

  if(aj == C && ai == O && ak == Si){	// to be modified completed
    prm.lam1 = 51.945;
    prm.gam1ij = 1.898;
    prm.gam1ik = 2.63341;
    prm.a1ij = 2.4;
    prm.a1ik = 2.75;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0001;
    prm.lam2 = 51.945;
    prm.gam2ij = 1.898;
    prm.gam2ik = 2.63341;
    prm.a2ij = 2.4;
    prm.a2ik = 2.75;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.0001;
    return &prm;
  }

  if(aj == O && ai == Si && ak == C){	// to be modified completed
    prm.lam1 = 37.0115;
    prm.gam1ij = 2.33945;
    prm.gam1ik = 1.69221;
    prm.a1ij = 2.75;
    prm.a1ik = 3.9;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0001;
    prm.lam2 = 37.0115;
    prm.gam2ij = 2.33945;
    prm.gam2ik = 1.69221;
    prm.a2ij = 2.75;
    prm.a2ik = 3.9;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.0001;
    return &prm;
  }

  if(aj == O && ai == C && ak == Si){	// to be modified completed
    prm.lam1 = 207.581;
    prm.gam1ij = -1.73465;
    prm.gam1ik = 11.0585;
    prm.a1ij = 2.4;
    prm.a1ik = 3.9;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0001;
    prm.lam2 = 207.581;
    prm.gam2ij = -1.73465;
    prm.gam2ik = 11.0585;
    prm.a2ij = 2.4;
    prm.a2ik = 3.9;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.0001;
    return &prm;
  }

  return &prm;
}

void calc_ESWSIC_force(void)
{
  int i, j, k, i1, i2;
  int neighbor, neighborj;
  double rij, rik, rjk, dxij, dyij, dzij, dxik, dyik, dzik, dxjk, dyjk, dzjk;
  double f, fij, fik, fjk, tmp_co;




/////////////////////////////////////////////////////////////////////////////
  double A, B, p, q, a, d, g, dg, dfc_SiO, dfc_CO;
/////////////////////////////////////////////////////////////////////////////




  double lam, gamij, gamik, aij, aik, cosjik, alpha, inv_rij, inv_rik, cos0;
  AtomIonType ai, aj, ak;
  ESWSIC2Body * prm2; ESWSIC3Body * prm3;
  double tmp[11], pe;


  // count_coordinate
  for(i=0;i<N;i++){
    pl[i]->co=0.0;  /* initialize */
    pl[i]->self_co = 0.0;
  }
  for(i=0;i<N;i++){
    if(pl[i]->atmtype==O){
      neighbor=pl[i]->neighbor;
      for(i1=1;i1<=neighbor;i1++){
        j=pl[i]->j[i1];




/////////////////////////////////////////////////////////
        if(pl[j]->atmtype==Si){
          tmp_co = fc_SiO(r1(i,j));
          pl[i]->co += tmp_co;
          pl[j]->co += tmp_co;

        if(pl[j]->atmtype==C){
          tmp_co = fc_CO(r1(i,j));
          pl[i]->co += tmp_co;
          pl[j]->co += tmp_co;
        }
/////////////////////////////////////////////////////////




        }
      }
    }
  }

  // two-body term
  for(i=0;i<N;i++){
    neighbor = pl[i]->neighbor;   // pl[i]->neighbor : number of neighbor atoms
    ai = pl[i]->atmtype;
    for(i1=1;i1<=neighbor;++i1){  // neighbor list is from 1 to "neighbor"
      j = pl[i]->j[i1];
      if(i < j){
        rij = pl[i]->d[i1];
        aj = pl[j]->atmtype;
        dxij = pl[i]->x - pl[j]->x;
        dyij = pl[i]->y - pl[j]->y;
        dzij = pl[i]->z - pl[j]->z; 
        mirror(&dxij,&dyij,&dzij);
        scale_to_real(&dxij,&dyij,&dzij);
        if((ai == Si && aj == Si) ||(ai == C && aj == C) || (ai == O  && aj == O)  ||
           (ai == Si && aj == C) || (ai == C && aj == Si) || ai == Omol || aj == Omol){
          prm2 = get_ESWSIC_2body(ai, aj);
          a = prm2->a;
          if(rij < a){
            A = prm2->A;
            B = prm2->B;
            p = prm2->p;
            q = prm2->q;
            d = prm2->d;
	    tmp[0] = 1.0/(rij-a);
	    tmp[1] = exp(d * tmp[0]);
	    tmp[2] = A*(B*pow(rij,-p)-pow(rij,-q));
	    tmp[3] = d * tmp[0]*tmp[0];
	    tmp[4] = A*(-B*p*pow(rij,-p-1.0)+q*pow(rij,-q-1.0));
	    pe = tmp[1]*tmp[2];
            f = (tmp[4] - tmp[2]*tmp[3])*tmp[1]/rij;
            pl[i]->pe2 += 0.5 * pe;
            pl[j]->pe2 += 0.5 * pe;
            pl[i]->fx -= f * dxij;
            pl[i]->fy -= f * dyij;
            pl[i]->fz -= f * dzij;
            pl[j]->fx += f * dxij;
            pl[j]->fy += f * dyij;
            pl[j]->fz += f * dzij;
            add_partial_stress(i, -f, dxij,dyij,dzij);
            add_partial_stress(j, -f, dxij,dyij,dzij);
	  } 
        }
        if(ai == Si && aj == O){
	  prm2 = get_ESWSIC_2body(ai, aj);
          a = prm2->a;
          if(rij < a){
            A = prm2->A;
            B = prm2->B;
            p = prm2->p;
            q = prm2->q;
            d = prm2->d;
            tmp[0] = 1.0/(rij-a);
            tmp[1] = exp(d * tmp[0]);
            tmp[2] = A*(B*pow(rij,-p)-pow(rij,-q));
            tmp[3] = d * tmp[0]*tmp[0];
            tmp[4] = A*(-B*p*pow(rij,-p-1.0)+q*pow(rij,-q-1.0));
            g = gO_SiO(pl[j]->co);
            pe = g * tmp[1]*tmp[2];
            pl[i]->pe2 += 0.5 * pe;
            pl[j]->pe2 += 0.5 * pe;

            f = g * (tmp[4] - tmp[2]*tmp[3])*tmp[1]/rij;
            pl[i]->fx -= f * dxij;
	    pl[i]->fy -= f * dyij;
	    pl[i]->fz -= f * dzij;
	    pl[j]->fx += f * dxij;
            pl[j]->fy += f * dyij;
            pl[j]->fz += f * dzij;
	    add_partial_stress(i, -f, dxij, dyij, dzij);
            add_partial_stress(j, -f, dxij, dyij, dzij);

            dg = d_gO_SiO(pl[j]->co);
            neighborj = pl[j]->neighbor;
            for(i2=1;i2<=neighborj;++i2){
              k = pl[j]->j[i2];
	      if(pl[k]->atmtype == Si){
	        rjk = pl[j]->d[i2];
	        dxjk = pl[j]->x - pl[k]->x;
	        dyjk = pl[j]->y - pl[k]->y;
	        dzjk = pl[j]->z - pl[k]->z;
                mirror(&dxjk, &dyjk, &dzjk);
	        scale_to_real(&dxjk, &dyjk, &dzjk);
                dfc_SiO = d_fc_SiO(rjk);
                f = dg * dfc_SiO * tmp[1]*tmp[2] / rjk;
                pl[j]->fx -= f * dxjk;
                pl[j]->fy -= f * dyjk;
                pl[j]->fz -= f * dzjk;
                pl[k]->fx += f * dxjk;
                pl[k]->fy += f * dyjk;
                pl[k]->fz += f * dzjk;
  		add_partial_stress(j, -f, dxjk,dyjk,dzjk);
  		add_partial_stress(k, -f, dxjk,dyjk,dzjk);
	      }
	    }
	  }
        }
        if(ai == O && aj == Si){
	  prm2 = get_ESWSIC_2body(ai, aj);
          a = prm2->a;
          if(rij < a){
            A = prm2->A;
            B = prm2->B;
            p = prm2->p;
            q = prm2->q;
            d = prm2->d;
            tmp[0] = 1.0/(rij-a);
            tmp[1] = exp(d * tmp[0]);
            tmp[2] = A*(B*pow(rij,-p)-pow(rij,-q));
            tmp[3] = d * tmp[0]*tmp[0];
            tmp[4] = A*(-B*p*pow(rij,-p-1.0)+q*pow(rij,-q-1.0));
            g = gO_SiO(pl[i]->co);
            pe = g * tmp[1]*tmp[2];
            pl[i]->pe2 += 0.5 * pe;
            pl[j]->pe2 += 0.5 * pe;

            f = g * (tmp[4] - tmp[2]*tmp[3])*tmp[1]/rij;
            pl[i]->fx -= f * dxij;
	    pl[i]->fy -= f * dyij;
	    pl[i]->fz -= f * dzij;
	    pl[j]->fx += f * dxij;
            pl[j]->fy += f * dyij;
            pl[j]->fz += f * dzij;
	    add_partial_stress(i, -f, dxij, dyij, dzij);
            add_partial_stress(j, -f, dxij, dyij, dzij);

            dg = d_gO_SiO(pl[i]->co);
            for(i2=1;i2<=neighbor;++i2){
              k = pl[i]->j[i2];
	      if(pl[k]->atmtype == Si){
	        rik = pl[i]->d[i2];
	        dxik = pl[i]->x - pl[k]->x;
	        dyik = pl[i]->y - pl[k]->y;
	        dzik = pl[i]->z - pl[k]->z;
                mirror(&dxik, &dyik, &dzik);
	        scale_to_real(&dxik, &dyik, &dzik);
                dfc_SiO = d_fc_SiO(rik);
                f = dg * dfc_SiO * tmp[1]*tmp[2] / rik;
                pl[i]->fx -= f * dxik;
                pl[i]->fy -= f * dyik;
                pl[i]->fz -= f * dzik;
                pl[k]->fx += f * dxik;
                pl[k]->fy += f * dyik;
                pl[k]->fz += f * dzik;
  		add_partial_stress(i, -f, dxik,dyik,dzik);
  		add_partial_stress(k, -f, dxik,dyik,dzik);
	      }
	    }
	  }
	}




//////////////////////////////////////////////////////////////
        if(ai == C && aj == O){
	  prm2 = get_ESWSIC_2body(ai, aj);
          a = prm2->a;
          if(rij < a){
            A = prm2->A;
            B = prm2->B;
            p = prm2->p;
            q = prm2->q;
            d = prm2->d;
            tmp[0] = 1.0/(rij-a);
            tmp[1] = exp(d * tmp[0]);
            tmp[2] = A*(B*pow(rij,-p)-pow(rij,-q));
            tmp[3] = d * tmp[0]*tmp[0];
            tmp[4] = A*(-B*p*pow(rij,-p-1.0)+q*pow(rij,-q-1.0));
            g = gO_CO(pl[i]->co);
            pe = g * tmp[1]*tmp[2];
            pl[i]->pe2 += 0.5 * pe;
            pl[j]->pe2 += 0.5 * pe;

            f = g * (tmp[4] - tmp[2]*tmp[3])*tmp[1]/rij;
            pl[i]->fx -= f * dxij;
	    pl[i]->fy -= f * dyij;
	    pl[i]->fz -= f * dzij;
	    pl[j]->fx += f * dxij;
            pl[j]->fy += f * dyij;
            pl[j]->fz += f * dzij;
	    add_partial_stress(i, -f, dxij, dyij, dzij);
            add_partial_stress(j, -f, dxij, dyij, dzij);

            dg = d_gO_CO(pl[i]->co);
            for(i2=1;i2<=neighbor;++i2){
              k = pl[i]->j[i2];
	      if(pl[k]->atmtype == C){
	        rik = pl[i]->d[i2];
	        dxik = pl[i]->x - pl[k]->x;
	        dyik = pl[i]->y - pl[k]->y;
	        dzik = pl[i]->z - pl[k]->z;
                mirror(&dxik, &dyik, &dzik);
	        scale_to_real(&dxik, &dyik, &dzik);
                dfc_SiO = d_fc_CO(rik);
                f = dg * dfc_CO * tmp[1]*tmp[2] / rik;
                pl[i]->fx -= f * dxik;
                pl[i]->fy -= f * dyik;
                pl[i]->fz -= f * dzik;
                pl[k]->fx += f * dxik;
                pl[k]->fy += f * dyik;
                pl[k]->fz += f * dzik;
  		add_partial_stress(i, -f, dxik,dyik,dzik);
  		add_partial_stress(k, -f, dxik,dyik,dzik);
	      }
	    }
	  }
	}

        if(ai == O && aj == C){
	  prm2 = get_ESWSIC_2body(ai, aj);
          a = prm2->a;
          if(rij < a){
            A = prm2->A;
            B = prm2->B;
            p = prm2->p;
            q = prm2->q;
            d = prm2->d;
            tmp[0] = 1.0/(rij-a);
            tmp[1] = exp(d * tmp[0]);
            tmp[2] = A*(B*pow(rij,-p)-pow(rij,-q));
            tmp[3] = d * tmp[0]*tmp[0];
            tmp[4] = A*(-B*p*pow(rij,-p-1.0)+q*pow(rij,-q-1.0));
            g = gO_CO(pl[i]->co);
            pe = g * tmp[1]*tmp[2];
            pl[i]->pe2 += 0.5 * pe;
            pl[j]->pe2 += 0.5 * pe;

            f = g * (tmp[4] - tmp[2]*tmp[3])*tmp[1]/rij;
            pl[i]->fx -= f * dxij;
	    pl[i]->fy -= f * dyij;
	    pl[i]->fz -= f * dzij;
	    pl[j]->fx += f * dxij;
            pl[j]->fy += f * dyij;
            pl[j]->fz += f * dzij;
	    add_partial_stress(i, -f, dxij, dyij, dzij);
            add_partial_stress(j, -f, dxij, dyij, dzij);

            dg = d_gO_CO(pl[i]->co);
            for(i2=1;i2<=neighbor;++i2){
              k = pl[i]->j[i2];
	      if(pl[k]->atmtype == C){
	        rik = pl[i]->d[i2];
	        dxik = pl[i]->x - pl[k]->x;
	        dyik = pl[i]->y - pl[k]->y;
	        dzik = pl[i]->z - pl[k]->z;
                mirror(&dxik, &dyik, &dzik);
	        scale_to_real(&dxik, &dyik, &dzik);
                dfc_CO = d_fc_CO(rik);
                f = dg * dfc_CO * tmp[1]*tmp[2] / rik;
                pl[i]->fx -= f * dxik;
                pl[i]->fy -= f * dyik;
                pl[i]->fz -= f * dzik;
                pl[k]->fx += f * dxik;
                pl[k]->fy += f * dyik;
                pl[k]->fz += f * dzik;
  		add_partial_stress(i, -f, dxik,dyik,dzik);
  		add_partial_stress(k, -f, dxik,dyik,dzik);
	      }
	    }
	  }
	}
//////////////////////////////////////////////////////////////




      }
    }
  }

  // three-body term
  for(i=0;i<N;i++){
    ai = pl[i]->atmtype;
    if(ai != Dummy){
      neighbor = pl[i]->neighbor;
      for(i1=1;i1<=neighbor;++i1){
        j = pl[i]->j[i1];
        aj = pl[j]->atmtype;
        if(aj != Dummy){
          rij = pl[i]->d[i1];
          for(i2=i1+1;i2<=neighbor;i2++){
            k = pl[i]->j[i2];
            ak = pl[k]->atmtype;
            if(ak != Dummy){
              rik = pl[i]->d[i2];
              dxij = pl[i]->x - pl[j]->x;
              dyij = pl[i]->y - pl[j]->y;
              dzij = pl[i]->z - pl[j]->z;
              mirror(&dxij, &dyij, &dzij);
              scale_to_real(&dxij, &dyij, &dzij);
              dxik = pl[i]->x - pl[k]->x;
              dyik = pl[i]->y - pl[k]->y;
              dzik = pl[i]->z - pl[k]->z;
              mirror(&dxik, &dyik, &dzik);
              scale_to_real(&dxik, &dyik, &dzik);
              dxjk = pl[j]->x - pl[k]->x;
              dyjk = pl[j]->y - pl[k]->y;
              dzjk = pl[j]->z - pl[k]->z;
              mirror(&dxjk, &dyjk, &dzjk);
              scale_to_real(&dxjk, &dyjk, &dzjk);
              rjk = sqrt(dxjk*dxjk + dyjk*dyjk + dzjk*dzjk);
              prm3 = get_ESWSIC_3body(aj, ai, ak);
              aij = prm3->a1ij;
              aik = prm3->a1ik;
              if((rij < aij) && (rik < aik)){
                lam   = prm3->lam1;
                gamij = prm3->gam1ij;
                gamik = prm3->gam1ik;
                cos0  = prm3->cos01;
                alpha = prm3->alpha1; 
                inv_rij = 1.0/rij;
                inv_rik = 1.0/rik;
                cosjik = (dxij*dxik+dyij*dyik+dzij*dzik)*inv_rij*inv_rik;
                tmp[0] = 1.0/(rij-aij);
                tmp[1] = 1.0/(rik-aik);
                tmp[2] = cosjik-cos0;
                tmp[3] = lam*exp(gamij*tmp[0]+gamik*tmp[1]);
                tmp[4] = tmp[2] * tmp[2] + alpha * pow(tmp[2],3.0);
                tmp[5] = 2.0 * tmp[2] + 3.0 * alpha * tmp[2] * tmp[2];
                tmp[6] = -gamij * tmp[0] * tmp[0] * tmp[3] * tmp[4];
                tmp[7] = -gamik * tmp[1] * tmp[1] * tmp[3] * tmp[4];
                tmp[8] = (rij-rik*cosjik)*inv_rij*inv_rik;
                tmp[9] = (rik-rij*cosjik)*inv_rij*inv_rik;
                tmp[10] = -rjk*inv_rij*inv_rik;
                fij = -(tmp[6]+tmp[3]*tmp[5]*tmp[8])*inv_rij;
                fik = -(tmp[7]+tmp[3]*tmp[5]*tmp[9])*inv_rik;
                fjk = -tmp[3]*tmp[5]*tmp[10]/rjk;
                pl[i]->pe3 += tmp[3]*tmp[4];
                pl[i]->fx +=  fij*dxij + fik*dxik;
                pl[i]->fy +=  fij*dyij + fik*dyik;
                pl[i]->fz +=  fij*dzij + fik*dzik;
                pl[j]->fx += -fij*dxij + fjk*dxjk;
                pl[j]->fy += -fij*dyij + fjk*dyjk;
                pl[j]->fz += -fij*dzij + fjk*dzjk;
                pl[k]->fx += -fik*dxik - fjk*dxjk;
                pl[k]->fy += -fik*dyik - fjk*dyjk;
                pl[k]->fz += -fik*dzik - fjk*dzjk;
		add_partial_stress(i, fij, dxij, dyij, dzij);
		add_partial_stress(j, fij, dxij, dyij, dzij);
                add_partial_stress(i, fik, dxik, dyik, dzik);
                add_partial_stress(k, fik, dxik, dyik, dzik);
                add_partial_stress(j, fjk, dxjk, dyjk, dzjk);
                add_partial_stress(k, fjk, dxjk, dyjk, dzjk);
	      }
              aij = prm3->a2ij;
              aik = prm3->a2ik;
              if((rij < aij) && (rik < aik)){
                lam = prm3->lam2;
                gamij = prm3->gam2ij;
                gamik = prm3->gam2ik;
                cos0 = prm3->cos02;
                alpha = prm3->alpha2; 
                inv_rij = 1.0/rij;
                inv_rik = 1.0/rik;
                cosjik = (dxij*dxik+dyij*dyik+dzij*dzik)*inv_rij*inv_rik;
                tmp[0] = 1.0/(rij-aij);
                tmp[1] = 1.0/(rik-aik);
                tmp[2] = cosjik-cos0;
                tmp[3] = lam*exp(gamij*tmp[0]+gamik*tmp[1]);
                tmp[4] = tmp[2] * tmp[2] + alpha * pow(tmp[2],3.0);
                tmp[5] = 2.0 * tmp[2] + 3.0 * alpha * tmp[2] * tmp[2];
                tmp[6] = -gamij * tmp[0] * tmp[0] * tmp[3] * tmp[4];
                tmp[7] = -gamik * tmp[1] * tmp[1] * tmp[3] * tmp[4];
                tmp[8] = (rij-rik*cosjik)*inv_rij*inv_rik;
                tmp[9] = (rik-rij*cosjik)*inv_rij*inv_rik;
                tmp[10] = -rjk*inv_rij*inv_rik;
                fij = -(tmp[6]+tmp[3]*tmp[5]*tmp[8])*inv_rij;
                fik = -(tmp[7]+tmp[3]*tmp[5]*tmp[9])*inv_rik;
                fjk = -tmp[3]*tmp[5]*tmp[10]/rjk;
                pl[i]->pe3 += tmp[3]*tmp[4];
                pl[i]->fx +=  fij*dxij + fik*dxik;
                pl[i]->fy +=  fij*dyij + fik*dyik;
                pl[i]->fz +=  fij*dzij + fik*dzik;
                pl[j]->fx += -fij*dxij + fjk*dxjk;
                pl[j]->fy += -fij*dyij + fjk*dyjk;
                pl[j]->fz += -fij*dzij + fjk*dzjk;
                pl[k]->fx += -fik*dxik - fjk*dxjk;
                pl[k]->fy += -fik*dyik - fjk*dyjk;
                pl[k]->fz += -fik*dzik - fjk*dzjk;
		add_partial_stress(i, fij, dxij, dyij, dzij);
		add_partial_stress(j, fij, dxij, dyij, dzij);
                add_partial_stress(i, fik, dxik, dyik, dzik);
                add_partial_stress(k, fik, dxik, dyik, dzik);
                add_partial_stress(j, fjk, dxjk, dyjk, dzjk);
                add_partial_stress(k, fjk, dxjk, dyjk, dzjk);
	      }
	    }
	  } 
        }
      }
    }
  }

  for(i=0;i<N;++i){
    for(j=0;j<3;++j){
      for(k=0;k<3;++k){
        sys.nfx[j][k] += pl[i]->str[j][k];
      }
    }
  }
}




///////////////////////////////////////////////////////////////////////////
static double fc_SiO(double r)
  // r : interatomic distance
{
  double value,disl,disr;

  disl = cutR_SiO - cutD_SiO;
  disr = cutR_SiO + cutD_SiO;

  if(r <= disl){
    value = 1.0;
  }
  else if(r >= disr){
    value = 0.0;
  }
  else{
    value = 0.5/PI*(sin(PI*(r-disl)/cutD_SiO))-0.5*(r-disl)/cutD_SiO+1.0;
  }
  return value;
}

static double fc_CO(double r)
  // r : interatomic distance
{
  double value,disl,disr;

  disl = cutR_CO - cutD_CO;
  disr = cutR_CO + cutD_CO;

  if(r <= disl){
    value = 1.0;
  }
  else if(r >= disr){
    value = 0.0;
  }
  else{
    value = 0.5/PI*(sin(PI*(r-disl)/cutD_CO))-0.5*(r-disl)/cutD_CO+1.0;
  }
  return value;
}
///////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////
static double d_fc_SiO(double r)
{
  double value,disl,disr;

  disl = cutR_SiO - cutD_SiO;
  disr = cutR_SiO + cutD_SiO;

  if((r > disl) && (r < disr)){
    value = 0.5/cutD_SiO * cos(PI*(r-disl)/cutD_SiO)-0.5/cutD_SiO;
  }
  else{
    value=0.0;
  }
  return value;
}

static double d_fc_CO(double r)
{
  double value,disl,disr;

  disl = cutR_CO - cutD_CO;
  disr = cutR_CO + cutD_CO;

  if((r > disl) && (r < disr)){
    value = 0.5/cutD_CO * cos(PI*(r-disl)/cutD_CO)-0.5/cutD_CO;
  }
  else{
    value=0.0;
  }
  return value;
}
///////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////
static double gO_SiO(double z)
{
  return M1_SiO/(exp((M2_SiO-z)/M3_SiO)+1.0)*exp(M4_SiO*(z-M5_SiO)*(z-M5_SiO));
}

static double d_gO_SiO(double z)
{
  return  M1_SiO/(exp((M2_SiO-z)/M3_SiO)+1.0)*exp(M4_SiO*(z-M5_SiO)*(z-M5_SiO))
          * (2.0*M4_SiO*(z-M5_SiO) + 1.0/M3_SiO/(1.0+exp((z-M2_SiO)/M3_SiO)));
}

static double gO_CO(double z)
{
  return M1_CO/(exp((M2_CO-z)/M3_CO)+1.0)*exp(M4_CO*(z-M5_CO)*(z-M5_CO));
}

static double d_gO_CO(double z)
{
  return  M1_CO/(exp((M2_CO-z)/M3_CO)+1.0)*exp(M4_CO*(z-M5_CO)*(z-M5_CO))
          * (2.0*M4_CO*(z-M5_CO) + 1.0/M3_CO/(1.0+exp((z-M2_CO)/M3_CO)));
}
///////////////////////////////////////////////////////////////////////////




} // end of namespace WASEDA_MD_LABO

/*** end of watanabe.c ***/

