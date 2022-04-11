/******************************************************************************
                                  eswge.cc

    --- Extended Stillinger-Weber Potential for Ge, O Mixed Systems ---

         Copyright (c) Takanobu Watanabe 2008  All Rights Reserved
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
} ESWGE2Body;

typedef struct {
  double lam1, lam2;
  double gam1ij, gam1ik, gam2ij, gam2ik;
  double a1ij, a1ik, a2ij, a2ik;
  double cos01, cos02;
  double alpha1, alpha2;
} ESWGE3Body;

static ESWGE2Body * get_ESWGE_2body(AtomIonType, AtomIonType);
static ESWGE3Body * get_ESWGE_3body(AtomIonType, AtomIonType, AtomIonType);
static double fc(double);
static double d_fc(double);
static double gO(double);
static double d_gO(double);

static double cutD, cutR, M1, M2, M3, M4, M5;

void set_ESWGE_parameter(void)
{
  printf("set_ESWGE_parameter\n");
  cutR = 2.75;
  cutD = 0.2;
  M1 = 0.01; M2 = 1.22; M3 = 0.34; M4 = 0.0315; M5 = 14.22;
}

double set_ESWGE_scut(void)
{
  double cut, tmpcut;
  int i,j;

  cut = 0.0;
  for(i=0;i<MAXELEMENT;i++){
    for(j=i;j<MAXELEMENT;j++){
      tmpcut = get_ESWGE_scutoff(AtomIonType(i), AtomIonType(j));
      if(cut < tmpcut) cut = tmpcut;
    }
  }

  return cut;
}

double get_ESWGE_scutoff(AtomIonType a1, AtomIonType a2)
{
  double cut = 0.0;
  int i;
  ESWGE2Body *param2;
  ESWGE3Body *param3;

  if(a1 == Ge && a2 == Ge) return 3.9258;
  if(a1 == Ge && a2 == O ) return 3.80;
  if(a1 == O  && a2 == Ge) return 3.80;
  if(a1 == O  && a2 == O ) return 2.30;
  if(a1 == Omol && a2 == Omol) return 3.80;
  if(a1 == Omol && a2 == O ) return 3.0;
  if(a1 == O  && a2 == Omol) return 3.0;
  if(a1 == Omol && a2 == Ge) return 3.30;
  if(a1 == Ge && a2 == Omol) return 3.30;

  param2 = get_ESWGE_2body(a1, a2);
  if(cut < param2->a) cut = param2->a;
  for(i=0;i<MAXELEMENT;i++){
    param3 = get_ESWGE_3body(a1, a2, (AtomIonType)i);
    if(cut < param3->a1ij) cut = param3->a1ij;
    if(cut < param3->a2ik) cut = param3->a2ik;
    param3 = get_ESWGE_3body(a2, a1, (AtomIonType)i);
    if(cut < param3->a1ij) cut = param3->a1ij;
    if(cut < param3->a2ik) cut = param3->a2ik;
  }

  return cut;
}

static ESWGE2Body * get_ESWGE_2body(AtomIonType ai, AtomIonType aj)
{
  static ESWGE2Body prm;

  prm.A = 0.0;
  prm.B = 0.0;
  prm.p = 0.0;
  prm.q = 0.0;
  prm.a = 0.1;
  prm.d = 1.0;

  if(ai == Ge && aj == Ge){
    prm.A = 13.6056;
    prm.B = 13.6264;
    prm.p = 4.0;
    prm.q = 0.0;
    prm.a = 3.9258;
    prm.d = 2.181;
    return &prm;
  }
  if(ai == O && aj == O){
    prm.A = -91.0;
    prm.B = 0.0;
    prm.p = 1.0;
    prm.q = 3.54;
    prm.a = 2.3;
    prm.d = 1.0;
    return &prm;
  }
  if((ai == Ge && aj == O) || (ai == O && aj == Ge)){
    prm.A = 45743.7;
    prm.B = 1.0016;
    prm.p = 2.8435;
    prm.q = 2.8388;
    prm.a = 2.75;
    prm.d = 1.0;
    return &prm;
  }
  if(ai == Omol && aj == Omol){
    prm.A = 175.392;
    prm.B = 0.9422;
    prm.p = 3.6;
    prm.q = 2.63;
    prm.a = 2.62;
    prm.d = 2.0951;
    return &prm;
  }
  if((ai == Omol && aj == O) || (ai == O && aj == Omol)){
    prm.A = -38.0;
    prm.B = 0.0;
    prm.p = 0.0;
    prm.q = 3.0;
    prm.a = 3.0;
    prm.d = 1.0;
    return &prm;
  }
  if((ai == Omol && aj == Ge) || (ai == Ge && aj == Omol)){
    prm.A = -82.0;
    prm.B = 0.0;
    prm.p = 0.0;
    prm.q = 3.95;
    prm.a = 3.3;
    prm.d = 1.0;
    return &prm;
  }

  return &prm;
}

static ESWGE3Body * get_ESWGE_3body(AtomIonType aj,
                                    AtomIonType ai,
                                    AtomIonType ak)
{
  static ESWGE3Body prm;

  prm.lam1 = 0.0;
  prm.lam2 = 0.0;
  prm.gam1ij = 1.0;
  prm.gam1ik = 1.0;
  prm.gam2ij = 1.0;
  prm.gam2ik = 1.0;
  prm.a1ij = 0.1;
  prm.a1ik = 0.1;
  prm.a2ij = 0.1;
  prm.a2ik = 0.1;
  prm.cos01 = 0.0;
  prm.cos02 = 0.0;
  prm.alpha1 = 0.0;
  prm.alpha2 = 0.0;

  if(aj == Ge && ai == O && ak == Ge){
    prm.lam1 = 41.0;
    prm.gam1ij = 1.5;
    prm.gam1ik = 1.5;
    prm.a1ij = 2.75;
    prm.a1ik = 2.75;
    prm.cos01 = -1.0;
    prm.alpha1 = 0.1;

    if(eswge_option == ESWGE_OPTION_133DEG){
      prm.lam2 = 60.0;
      prm.gam2ij = 1.5;
      prm.gam2ik = 1.5;
      prm.a2ij = 2.75;
      prm.a2ik = 2.75;
      prm.cos02 = -0.69;
      prm.alpha2 = 2.2;
    }
    if(eswge_option == ESWGE_OPTION_144DEG){
      prm.lam2 = 143.0;
      prm.gam2ij = 1.5;
      prm.gam2ik = 1.5;
      prm.a2ij = 2.75;
      prm.a2ik = 2.75;
      prm.cos02 = -0.809;
      prm.alpha2 = 3.40;
    }

    return &prm;
  }
  if(aj == O && ai == Ge && ak == O){
    prm.lam1 = 17.0;
    prm.gam1ij = 1.0;
    prm.gam1ik = 1.0;
    prm.a1ij = 2.75;
    prm.a1ik = 2.75;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.5;

    //    prm.lam2 = 1.5;
    prm.lam2 = 0.0;

    prm.gam2ij = 1.0;
    prm.gam2ik = 1.0;
    prm.a2ij = 3.8;
    prm.a2ik = 3.8;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 0.5;

    return &prm;
  }
  if(aj == Ge && ai == Ge && ak == Ge){
    prm.lam1 = 59.83;
    prm.gam1ij = 2.51;
    prm.gam1ik = 2.51;
    prm.a1ij = 3.9258;
    prm.a1ik = 3.9258;
    prm.cos01 = -1.0/3.0;
    prm.alpha1 = 0.0;

    return &prm;
  }
  if(aj == Ge && ai == Ge && ak == O){
    if(eswge_option == ESWGE_OPTION_133DEG){
      prm.lam1 = 125.0;
      prm.gam1ij = 1.0;
      prm.gam1ik = 1.0;
      prm.a1ij = 3.0;
      prm.a1ik = 2.2;
      prm.cos01 = -1.0/3.0;
      prm.alpha1 = 0.0;
    }
    if(eswge_option == ESWGE_OPTION_144DEG){
      prm.lam1 = 200.0;
      prm.gam1ij = 1.3;
      prm.gam1ik = 1.0;
      prm.a1ij = 3.0;
      prm.a1ik = 2.2;
      prm.cos01 = -1.0/3.0;
      prm.alpha1 = 0.0;
    }
    prm.lam2 = 3.7;
    prm.gam2ij = 1.9;
    prm.gam2ik = 1.0;
    prm.a2ij = 3.9;
    prm.a2ik = 2.2;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 1.2;

    return &prm;
  }
  if(aj == O && ai == Ge && ak == Ge){
    if(eswge_option == ESWGE_OPTION_133DEG){
      prm.lam1 = 125.0;
      prm.gam1ij = 1.0;
      prm.gam1ik = 1.0;
      prm.a1ij = 2.2;
      prm.a1ik = 3.0;
      prm.cos01 = -1.0/3.0;
      prm.alpha1 = 0.0;
    }
    if(eswge_option == ESWGE_OPTION_144DEG){
      prm.lam1 = 200.0;
      prm.gam1ij = 1.0;
      prm.gam1ik = 1.3;
      prm.a1ij = 2.2;
      prm.a1ik = 3.0;
      prm.cos01 = -1.0/3.0;
      prm.alpha1 = 0.0;
    }
    prm.lam2 = 3.7;
    prm.gam2ij = 1.0;
    prm.gam2ik = 1.9;
    prm.a2ij = 2.2;
    prm.a2ik = 3.9;
    prm.cos02 = -1.0/3.0;
    prm.alpha2 = 1.2;

    return &prm;
  }

  return &prm;
}

void calc_ESWGE_force(void)
{
  int i, j, k, i1, i2;
  int neighbor, neighborj;
  double rij, rik, rjk, dxij, dyij, dzij, dxik, dyik, dzik, dxjk, dyjk, dzjk;
  double f, fij, fik, fjk, tmp_co;
  double A, B, p, q, a, d, g, dg, dfc;
  double lam, gamij, gamik, aij, aik, cosjik, alpha, inv_rij, inv_rik, cos0;
  AtomIonType ai, aj, ak;
  ESWGE2Body * prm2; ESWGE3Body * prm3;
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
        if(pl[j]->atmtype==Ge){
          tmp_co = fc(r1(i,j));
          pl[i]->co += tmp_co;
          pl[j]->co += tmp_co;
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
        if((ai == Ge && aj == Ge) || (ai == O  && aj == O) || ai == Omol
           || aj == Omol){
          prm2 = get_ESWGE_2body(ai, aj);
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
        if(ai == Ge && aj == O){
	  prm2 = get_ESWGE_2body(ai, aj);
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
            g = gO(pl[j]->co);
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

            dg = d_gO(pl[j]->co);
            neighborj = pl[j]->neighbor;
            for(i2=1;i2<=neighborj;++i2){
              k = pl[j]->j[i2];
	      if(pl[k]->atmtype == Ge){
	        rjk = pl[j]->d[i2];
	        dxjk = pl[j]->x - pl[k]->x;
	        dyjk = pl[j]->y - pl[k]->y;
	        dzjk = pl[j]->z - pl[k]->z;
                mirror(&dxjk, &dyjk, &dzjk);
	        scale_to_real(&dxjk, &dyjk, &dzjk);
                dfc = d_fc(rjk);
                f = dg * dfc * tmp[1]*tmp[2] / rjk;
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
        if(ai == O && aj == Ge){
	  prm2 = get_ESWGE_2body(ai, aj);
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
            g = gO(pl[i]->co);
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

            dg = d_gO(pl[i]->co);
            for(i2=1;i2<=neighbor;++i2){
              k = pl[i]->j[i2];
	      if(pl[k]->atmtype == Ge){
	        rik = pl[i]->d[i2];
	        dxik = pl[i]->x - pl[k]->x;
	        dyik = pl[i]->y - pl[k]->y;
	        dzik = pl[i]->z - pl[k]->z;
                mirror(&dxik, &dyik, &dzik);
	        scale_to_real(&dxik, &dyik, &dzik);
                dfc = d_fc(rik);
                f = dg * dfc * tmp[1]*tmp[2] / rik;
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
              prm3 = get_ESWGE_3body(aj, ai, ak);
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

static double fc(double r)
  // r : interatomic distance
{
  double value,disl,disr;

  disl = cutR - cutD;
  disr = cutR + cutD;

  if(r <= disl){
    value = 1.0;
  }
  else if(r >= disr){
    value = 0.0;
  }
  else{
    value = 0.5/PI*(sin(PI*(r-disl)/cutD))-0.5*(r-disl)/cutD+1.0;
  }
  return value;
}

static double d_fc(double r)
{
  double value,disl,disr;

  disl = cutR - cutD;
  disr = cutR + cutD;

  if((r > disl) && (r < disr)){
    value = 0.5/cutD * cos(PI*(r-disl)/cutD)-0.5/cutD;
  }
  else{
    value=0.0;
  }
  return value;
}

static double gO(double z)
{
  return M1/(exp((M2-z)/M3)+1.0)*exp(M4*(z-M5)*(z-M5));
}

static double d_gO(double z)
{
  return  M1/(exp((M2-z)/M3)+1.0)*exp(M4*(z-M5)*(z-M5))
          * (2.0*M4*(z-M5) + 1.0/M3/(1.0+exp((z-M2)/M3)));
}

} // end of namespace WASEDA_MD_LABO

/*** end of watanabe.c ***/
