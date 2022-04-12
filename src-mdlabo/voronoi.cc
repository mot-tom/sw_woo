/******************************************************************************
                                  voronoi.cc

                      --- Voronoi Sphere Calculator ---

         Copyright (c) Takanobu Watanabe 2008  All Rights Reserved
******************************************************************************/
#include "mddef.h"

namespace WASEDA_MD_LABO {

typedef struct TmpVoronoiSphere{
  double x,y,z,d;
  bool included;
  struct TmpVoronoiSphere *next;
} VoronoiSphere;

void calc_voronoi(void)
{
  int i,j,k,l,m,a;
  int b1, b2, b3;
  double **mat, **invmat, **tmp;
  double *vec, *ans;
  double volume;
  double xi, yi, zi, xk, yk, zk, xl, yl, zl, xm, ym, zm, xv, yv, zv, dv;
  double sxv, syv, szv;
  double dx, dy, dz;
  int *neighbors, num_neighbors;
  VoronoiSphere *vs, *vshead, *vsj, *ld, *sd;
  int num_voronoi_spheres;

  printf("***** Start Voronoi Sphere Search *****\n");
  printf("  cut off distance : %f [angstrom]\n", voronoi_cut);

  matrix_invers(ul,iul,3);
  mat = malloc_for_matrix(4);
  invmat = malloc_for_matrix(4);
  tmp = malloc_for_matrix(3);
  vec = (double *)malloc(4*sizeof(double));
  ans = (double *)malloc(4*sizeof(double));
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      mat[i][j] = 0.0;
      invmat[i][j] = 0.0;
    }
    mat[i][3] = 1.0;
    vec[i] = 0.0;
    ans[i] = 0.0;
  }
  vshead = (VoronoiSphere *)malloc(sizeof(VoronoiSphere));
  vshead->next = NULL;
  vs = vshead;

  num_voronoi_spheres = 0;
  volume = 0.0;
  neighbors = NULL;
  for(i=0;i<N-3;i++){
    printf("Voronoi search : atom %d/%d\n",i, N);
    // make neighbor list
    num_neighbors = 0;
    for(j=i+1;j<N;j++){
      if(r1(i,j) < voronoi_cut/sig) num_neighbors++; 
    }
    if(neighbors) free(neighbors);
    neighbors = (int *)malloc(num_neighbors * sizeof(int));
    a = 0;
    for(j=i+1;j<N;j++){
      if(r1(i,j) < voronoi_cut/sig) neighbors[a++] = j; 
    }

    // search voronoi sphere
    for(b1=0;b1<num_neighbors-2;b1++){
      k = neighbors[b1];
      xi = pl[i]->x; yi = pl[i]->y; zi = pl[i]->z;
      scale_to_real(&xi,&yi,&zi);      
      mat[0][0] = xi; mat[0][1] = yi; mat[0][2] = zi;
      vec[0] = -1.0 *(xi*xi + yi*yi + zi*zi);
      dx = pl[k]->x - pl[i]->x;
      dy = pl[k]->y - pl[i]->y;
      dz = pl[k]->z - pl[i]->z;
      mirror(&dx,&dy,&dz);
      scale_to_real(&dx,&dy,&dz);
      xk = xi + dx; yk = yi + dy; zk = zi + dz;
      mat[1][0] = xk; mat[1][1] = yk; mat[1][2] = zk;
      vec[1] = -1.0 *(xk*xk + yk*yk + zk*zk);
      for(b2=b1+1;b2<num_neighbors-1;b2++){
        l = neighbors[b2];
        dx = pl[l]->x - pl[i]->x;
        dy = pl[l]->y - pl[i]->y;
        dz = pl[l]->z - pl[i]->z;
        mirror(&dx,&dy,&dz);
        scale_to_real(&dx,&dy,&dz);
        xl = xi + dx; yl = yi + dy; zl = zi + dz;
        mat[2][0] = xl; mat[2][1] = yl; mat[2][2] = zl;
        vec[2] = -1.0 *(xl*xl + yl*yl + zl*zl);
        for(b3=b2+1;b3<num_neighbors;b3++){
          m = neighbors[b3];
          dx = pl[m]->x - pl[i]->x;
          dy = pl[m]->y - pl[i]->y;
          dz = pl[m]->z - pl[i]->z;
          mirror(&dx,&dy,&dz);
          scale_to_real(&dx,&dy,&dz);
          xm = xi + dx; ym = yi + dy; zm = zi + dz;
          mat[3][0] = xm; mat[3][1] = ym; mat[3][2] = zm;
          vec[3] = -1.0 *(xm*xm + ym*ym + zm*zm);

          // calc center of Voronoi sphere
	  // print_matrix(mat, 4);

          tmp[0][0] = xk - xi; tmp[0][1] = yk - yi; tmp[0][2] = zk - zi;
          tmp[1][0] = xl - xi; tmp[1][1] = yl - yi; tmp[1][2] = zl - zi;
          tmp[2][0] = xm - xi; tmp[2][1] = ym - yi; tmp[2][2] = zm - zi;
          volume = fabs(cell_volume(tmp));

          if(fabs(volume) >  0.000001){
            matrix_invers(mat, invmat, 4);
            matrix_vector_product(invmat, vec, ans, 4);          
            xv = -0.5 * ans[0];
            yv = -0.5 * ans[1];
            zv = -0.5 * ans[2];
            dv = sqrt(xv*xv + yv*yv + zv*zv - ans[3]);
            sxv = xv*iul[0][0] + yv*iul[0][1] + zv*iul[0][2];
	    syv = xv*iul[1][0] + yv*iul[1][1] + zv*iul[1][2];
	    szv = xv*iul[2][0] + yv*iul[2][1] + zv*iul[2][2];

            // check
            for(j=0; j<N; j++){
              if(j != i && j != k && j != l && j != m){
                dx = pl[j]->x - sxv;
                dy = pl[j]->y - syv;
                dz = pl[j]->z - szv;
                mirror(&dx,&dy,&dz);
                scale_to_real(&dx,&dy,&dz);
                if(dv > sqrt(dx*dx + dy*dy + dz*dz)){
                  break;
	        }
	      }
	    }
            if(j == N){ // no atom inside the sphere
              // register voronoi sphere
              vs->x = xv;
              vs->y = yv;
              vs->z = zv;
              vs->d = dv;
              vs->included = false;
	      // printf("sv[%f, %f, %f]\n",sxv, syv, szv);
              // printf("v[%f, %f, %f] %f\n",
	      // xv*sig,yv*sig,zv*sig,4.0*PI/3.0*dv*dv*dv*sig*sig*sig);
              vs->next = (VoronoiSphere *)malloc(sizeof(VoronoiSphere));
              vs = vs->next;
              vs->next = NULL;
              num_voronoi_spheres++;
	    }

            // add the tetrahedron volume
	  }
	}
      }
    }
  }
  printf("Voronoi search completed.\n");
     
  // check included sphere
  vs = vshead;
  for(i=0;i<num_voronoi_spheres;i++){
    if(vs->included == false){
      vsj = vs->next;
      while(vsj->next){
        if(vsj->included == false){
          if(vs->d > vsj->d){
            ld = vs;
            sd = vsj;
	  } else {
            sd = vs;
            ld = vsj;
	  }
          dx = vsj->x - vs->x;
          dy = vsj->y - vs->y;
          dz = vsj->z - vs->z;
          if(ld->d > sqrt(dx*dx + dy*dy + dz*dz)){
            sd->included = true;
	  }
	}
        vsj = vsj->next;
      }
    }
    vs = vs->next;
  }

  printf("number of voronoi spheres : %d\n", num_voronoi_spheres);

  // print voronoi spheres
  printf("###############################################\n");
  printf(" Summary of Voronoi Analysis\n");
  printf("###############################################\n");
  printf("  xv  yv  xv  diameter [angstrom]  volume [angstrom^3]\n");
  vs = vshead;
  while(vs->next){
    if(vs->included == false){
      printf("%f %f %f %f %f\n", vs->x * sig, vs->y * sig, vs->z * sig,
	     vs->d * sig,
             4.0 * PI / 3.0 * pow(vs->d*sig, 3.0));
    }
    vs = vs->next;
  }
  printf("###############################################\n");
}

} // end of namespace WASEDA_MD_LABO

//*** end of voronoi.cc 
