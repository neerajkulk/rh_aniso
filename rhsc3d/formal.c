/* ------- file: -------------------------- formal.c ----------------
       Version:       rh2.0, 3-D Cartesian, short characteristics
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Jun  6 16:52:25 2018 --

       --------------------------                      ----------RH-- */

/* --- Formal solution with given source function and PRD emission
   profile. --                                     -------------- */


#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "error.h"
#include "background.h"
#include "inputs.h"
#include "constant.h"
#include "statistics.h"
#include "xdr.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Geometry geometry;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];

/* ------- begin -------------------------- Formal.c ---------------- */

double Formal(int nspect, bool_t eval_operator, bool_t redistribute)
{
  const char routineName[] = "Formal";
  register int k, l, mu, n, nact;

  bool_t   initialize, boundbound, polarized_as, polarized_c,
    PRD_angle_dep, to_obs, solveStokes, angle_dep;
  enum     FeautrierOrder F_order;
  int      Nspace = atmos.Nspace, Nrays = atmos.Nrays, nt,la;

  double  *phi, *I, *chi, *S, **Ipol, **Spol, *Psi, *Jdag, wmu, dJmax, dJ,
    *eta_Q, *eta_U, *eta_V, *eta_c_Q, *eta_c_U, *eta_c_V, *J20dag, musq, threemu1, threemu2, *J, *J20, *reJ21, *imJ21, *reJ22, *imJ22, inc, azi, wlambda, domg_dlam;

  ActiveSet *as;
  AtomicLine *line;
  Atom *atom;
  pthread_mutex_t *rate_lock;


  /* --- Retrieve active set as of transitions at wavelength nspect - */

  as = &spectrum.as[nspect];
  nt = nspect % input.Nthreads;
  alloc_as(nspect, eval_operator);

  //printf("formal.c is running again \n");

  /* --- Check whether current active set includes a bound-bound
	 and/or polarized transition and/or angle-dependent PRD
	 transition, and/or polarization through background scattering.
	 Otherwise, only angle-independent opacity and source functions
	 are needed --                                 -------------- */

  /* --- Check for bound-bound transition in active set -- ---------- */

  boundbound    = containsBoundBound(as);

  /* --- Check for line with angle-dependent PRD in set -- ---------- */

  PRD_angle_dep = (containsPRDline(as) && input.PRD_angle_dep);

  /* --- Check for polarized bound-bound transition in active set - - */

  polarized_as  = containsPolarized(as);

  /* --- Check for polarized bound-bound transition in background - - */

  polarized_c   = atmos.backgrflags[nspect].ispolarized;

  /* --- Determine if we solve for I, or for I, Q, U, V -- ---------- */

  solveStokes   = (input.StokesMode == FULL_STOKES &&
		   (polarized_as || polarized_c || input.backgr_pol));

  /* --- Determine if we have to do angle-dependent opacity and
     emissivity --                                 -------------- */

  angle_dep     = (polarized_as || polarized_c || PRD_angle_dep ||
		   (input.backgr_pol && input.StokesMode == FULL_STOKES) ||
		   (atmos.moving &&
		    (boundbound || atmos.backgrflags[nspect].hasline)));

  /* --- Allocate temporary space --                   -------------- */

  if (eval_operator)
    Psi = (double *) malloc(Nspace * sizeof(double));
  else
    Psi = NULL;

  if (solveStokes) {
    Spol = matrix_double(4, Nspace);
    S    = Spol[0];
    Ipol = matrix_double(4, Nspace);
    I    = Ipol[0];
  } else {
    S = (double *) malloc(Nspace * sizeof(double));
    I = (double *) malloc(Nspace * sizeof(double));
  }
  chi  = (double *) malloc(Nspace * sizeof(double));

  /* --- Store current mean intensity, initialize new one to zero - - */

  Jdag = (double *) malloc(Nspace * sizeof(double));
  if (input.limit_memory) {
    J = (double *) malloc(atmos.Nspace * sizeof(double));
    readJlambda(nspect, Jdag);
  } else {
    J = spectrum.J[nspect];
    for (k = 0;  k < Nspace;  k++) Jdag[k] = J[k];
  }
  if (spectrum.updateJ)
    for (k = 0;  k < Nspace;  k++) J[k] = 0.0;

  /* --- Store current anisotropy, initialize new one to zero ---- -- */

  if (input.backgr_pol) {
    J20dag = (double *) malloc(Nspace * sizeof(double));
    if (input.limit_memory) {
      J20 = (double *) malloc(Nspace * sizeof(double));
      reJ21 = (double *) malloc(Nspace * sizeof(double));
      imJ21 = (double *) malloc(Nspace * sizeof(double));
      reJ22 = (double *) malloc(Nspace * sizeof(double));
      imJ22 = (double *) malloc(Nspace * sizeof(double));

      readJ20lambda(nspect, J20dag);
    } else {
      J20 = spectrum.J20[nspect];
      reJ21 = spectrum.reJ21[nspect];
      imJ21 = spectrum.imJ21[nspect];
      reJ22 = spectrum.reJ22[nspect];
      imJ22 = spectrum.imJ22[nspect];

      for (k = 0;  k < Nspace;  k++)
	J20dag[k] = J20[k];
    }
    if (spectrum.updateJ)
      for (k = 0;  k < Nspace;  k++){
	J20[k]   = 0.0;
	reJ21[k] = 0.0;
	imJ21[k] = 0.0;
	reJ22[k] = 0.0;
	imJ22[k] = 0.0;
      }

  }
  /* --- Case of angle-dependent opacity and source function -- ----- */

  if (angle_dep) {
    for (mu = 0;  mu < Nrays;  mu++) {
      wmu = 0.5 * geometry.wmu[mu];
      if (input.backgr_pol) {

	azi   = atan2(geometry.muy[mu],geometry.mux[mu]);
	inc = geometry.muz[mu];                   /* inc is the same is mu*/

	musq = SQ(geometry.muz[mu]);
	threemu1 = TWOSQRTTWO * (3.0*musq - 1.0);
	threemu2 = (3.0 * TWOSQRTTWO) * (musq - 1.0);



      }
      for (to_obs = 0;  to_obs <= 1;  to_obs++) {
	initialize = (mu == 0 && to_obs == 0);

	if (initialize || atmos.backgrflags[nspect].hasline)
	  readBackground(nspect, mu, to_obs);

	if (initialize || boundbound)
	  Opacity(nspect, mu, to_obs, initialize);

	if (eval_operator) addtoCoupling(nspect);
	for (k = 0;  k < Nspace;  k++) {
	  chi[k] = as->chi[k] + as->chi_c[k];
	  S[k]   = as->eta[k] + as->eta_c[k] + as->sca_c[k]*Jdag[k];
	}
	if (solveStokes) {
	  for (k = Nspace;  k < 4*Nspace;  k++) Spol[0][k] = 0.0;

	  /* --- Add emissivity due to active set for Q, U, V -- ---- */

	  if (polarized_as) {
	    for (k = Nspace;  k < 4*Nspace;  k++)
	      Spol[0][k] += as->eta[k];
	  }
	  /* --- Add emissivity due to background lines -- ---------- */

	  if (polarized_c) {
	    for (k = Nspace;  k < 4*Nspace;  k++)
	      Spol[0][k] += as->eta_c[k];
	  }
	  /* --- Add emissivity due to background scattering -- ----- */

	  if (input.backgr_pol && input.StokesMode == FULL_STOKES) {
	    for (k = 0;  k < Nspace;  k++) {
	      Spol[0][k] += threemu1 * as->sca_c[k]*J20dag[k];
	      Spol[1][k] += threemu2 * as->sca_c[k]*J20dag[k];
	    }
	  }
	  for (n = 0;  n < 4;  n++) {
	    for (k = 0;  k < Nspace;  k++)
	      Spol[n][k] /= chi[k];
	  }
	  ShortChar_Stokes(&geometry, nspect, mu, to_obs,
			   chi, Spol, Ipol, Psi);
	} else {
	  for (k = 0;  k < Nspace;  k++)
	    S[k] /= chi[k];
	  ShortChar(&geometry, nspect, mu, to_obs, chi, S, I, Psi);
	}
	if (eval_operator) {
	  for (k = 0;  k < Nspace;  k++) Psi[k] /= chi[k];
	  addtoGamma(nspect, wmu, I, Psi);
	}

	if (spectrum.updateJ) {

	  /* --- Accumulate mean intensity and rates -- ----------- */

	  for (k = 0;  k < Nspace;  k++)
	    J[k] += wmu * I[k];
	  addtoRates(nspect, mu, to_obs, wmu, I, redistribute);

	  /* --- Accumulate anisotropy --            -------------- */

	  if (PRD_angle_dep) writeImu(nspect, mu, to_obs, I);
	}
      }
      /* --- Save emergent intensity --              -------------- */

      for (l = 0;  l < geometry.Nplane;  l++)
	spectrum.I[nspect*Nrays + mu][l] = I[l];

      if (solveStokes) {
	for (l = 0;  l < geometry.Nplane;  l++) {
	  spectrum.Stokes_Q[nspect*Nrays + mu][l] = Ipol[1][l];
	  spectrum.Stokes_U[nspect*Nrays + mu][l] = Ipol[2][l];
	  spectrum.Stokes_V[nspect*Nrays + mu][l] = Ipol[3][l];
	}
      }
    }
  } else {

    /* --- The angle-independent case --               -------------- */

    readBackground(nspect, 0, 0);
    Opacity(nspect, 0, 0, initialize=TRUE);
    if (eval_operator) addtoCoupling(nspect);

    for (k = 0;  k < Nspace;  k++) {
      chi[k] = as->chi[k] + as->chi_c[k];
      S[k]   = (as->eta[k] +
		as->eta_c[k] + as->sca_c[k]*Jdag[k]) / chi[k];
    }
    for (mu = 0;  mu < Nrays;  mu++) {
      wmu = 0.5 * geometry.wmu[mu];
      for (to_obs = 0;  to_obs <= 1;  to_obs++) {
	ShortChar(&geometry, nspect, mu, to_obs, chi, S, I, Psi);

	if (eval_operator) {
	  for (k = 0;  k < Nspace;  k++) Psi[k] /= chi[k];
	  addtoGamma(nspect, wmu, I, Psi);
	}
	if (spectrum.updateJ) {
	  for (k = 0;  k < Nspace;  k++)
	    J[k] += I[k] * wmu;
	  addtoRates(nspect, mu, to_obs, wmu, I, redistribute);
	}
      }
      for (l = 0;  l < geometry.Nplane;  l++)
	spectrum.I[nspect*Nrays + mu][l] = I[l];
    }
  }
  /* --- Write new J for current position in the spectrum -- -------- */

  dJmax = 0.0;
  if (spectrum.updateJ) {
    for (k = 0;  k < Nspace;  k++) {
      dJ = fabs(1.0 - Jdag[k]/J[k]);
      dJmax = MAX(dJmax, dJ);
    }
    if (input.limit_memory) writeJlambda(nspect, J);
  }

  /* --- Clean up --                                 ---------------- */

  free_as(nspect, eval_operator);
  if (eval_operator) free(Psi);

  free(chi);
  if (solveStokes) {
    freeMatrix((void **) Ipol);
    freeMatrix((void **) Spol);
  } else {
    free(I);
    free(S);
  }

  free(Jdag);
  if (input.limit_memory) free(J);
  if (input.backgr_pol) {
    free(J20dag);
    if (input.limit_memory) {
      free(J20);
      free(reJ21);
      free(imJ21);
      free(reJ22);
      free(imJ22);

    }
  }

  return dJmax;
}
/* ------- end ---------------------------- Formal.c ---------------- */

/* ------- begin -------------------------- Hydrostatic.c ----------- */

void Hydrostatic(int NmaxIter, double iterLimit)
{
  const char routineName[] = "Hydrostatic";

  if (atmos.hydrostatic) {
    sprintf(messageStr,
	    "Can only establish hydrostatic equilibrium in 1-D geometry");
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
}
/* ------- end ---------------------------- Hydrostatic.c ----------- */

void write_spherical_tensor()
{
  printf("\n \n \n Computing and writing spherical tensors \n \n \n");

  int nspect;
  bool_t eval_operator = FALSE;
  bool_t redistribute = FALSE;

  register int k, l, mu, n, nact;

  bool_t   initialize, boundbound, polarized_as, polarized_c,
    PRD_angle_dep, to_obs, solveStokes, angle_dep;
  enum     FeautrierOrder F_order;
  int      Nspace = atmos.Nspace, Nrays = atmos.Nrays, nt,la,ni,nf;

  double  *phi, *I, *chi, *S, **Ipol, **Spol, *Psi, *Jdag, wmu, dJmax, dJ,
    *eta_Q, *eta_U, *eta_V, *eta_c_Q, *eta_c_U, *eta_c_V, *J20dag, musq, threemu1, threemu2, *J, *J20, *reJ21, *imJ21, *reJ22, *imJ22, inc, azi, wlambda, domg_dlam;

  ActiveSet *as;
  AtomicLine *line;
  Atom *atom;
  pthread_mutex_t *rate_lock;

  /* Loop over all wavelengths to find indices corresponding to strontium 4607 wavelengths */

  for (nspect = 0; nspect <= spectrum.Nspect - 1; nspect++){
    as = &spectrum.as[nspect];
    for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
      atom = atmos.activeatoms[nact];
      for (n = 0;  n < as->Nactiveatomrt[nact];  n++) {

	line = as->art[nact][n].ptype.line;

	if (fabs(line->lambda0 - 460.861732) < 0.00001)
	  {
	    ni = line->Nblue; // first index of 4607 line
	    nf = ni + line->Nlambda - 1; // last index of 4607 line
	  }
      }
    }
  }

  /* Now start integration over just 4607 wavelengths */

  for (nspect =ni; nspect <=nf; nspect++){  // For loop over SR wavelengths

    printf("wavelength = %.20f \n", spectrum.lambda[nspect]/1.00029);

    /* --- Retrieve active set as of transitions at wavelength nspect - */

    as = &spectrum.as[nspect];
    nt = nspect % input.Nthreads;
    alloc_as(nspect, eval_operator);

    /* --- Check whether current active set includes a bound-bound
       and/or polarized transition and/or angle-dependent PRD
       transition, and/or polarization through background scattering.
       Otherwise, only angle-independent opacity and source functions
       are needed --                                 -------------- */

    /* --- Check for bound-bound transition in active set -- ---------- */

    boundbound    = containsBoundBound(as);

    /* --- Check for line with angle-dependent PRD in set -- ---------- */

    PRD_angle_dep = (containsPRDline(as) && input.PRD_angle_dep);

    /* --- Check for polarized bound-bound transition in active set - - */

    polarized_as  = containsPolarized(as);

    /* --- Check for polarized bound-bound transition in background - - */

    polarized_c   = atmos.backgrflags[nspect].ispolarized;

    /* --- Determine if we solve for I, or for I, Q, U, V -- ---------- */

    solveStokes   = (input.StokesMode == FULL_STOKES &&
		     (polarized_as || polarized_c || input.backgr_pol));

    /* --- Determine if we have to do angle-dependent opacity and
       emissivity --                                 -------------- */

    angle_dep     = (polarized_as || polarized_c || PRD_angle_dep ||
		     (input.backgr_pol && input.StokesMode == FULL_STOKES) ||
		     (atmos.moving &&
		      (boundbound || atmos.backgrflags[nspect].hasline)));

    /* --- Allocate temporary space --                   -------------- */

    if (eval_operator)
      Psi = (double *) malloc(Nspace * sizeof(double));
    else
      Psi = NULL;

    if (solveStokes) {
      Spol = matrix_double(4, Nspace);
      S    = Spol[0];
      Ipol = matrix_double(4, Nspace);
      I    = Ipol[0];
    } else {
      S = (double *) malloc(Nspace * sizeof(double));
      I = (double *) malloc(Nspace * sizeof(double));
    }
    chi  = (double *) malloc(Nspace * sizeof(double));

    /* --- Store current mean intensity, initialize new one to zero - - */

    Jdag = (double *) malloc(Nspace * sizeof(double));
    if (input.limit_memory) {
      J = (double *) malloc(atmos.Nspace * sizeof(double));
      readJlambda(nspect, Jdag);
    } else {
      J = spectrum.J[nspect];
      for (k = 0;  k < Nspace;  k++) Jdag[k] = J[k];
    }
    if (spectrum.updateJ)
      for (k = 0;  k < Nspace;  k++) J[k] = 0.0;

    /* --- Store current anisotropy, initialize new one to zero ---- -- */

    if (input.backgr_pol) {
      J20dag = (double *) malloc(Nspace * sizeof(double));
      if (input.limit_memory) {
	J20 = (double *) malloc(Nspace * sizeof(double));
	reJ21 = (double *) malloc(Nspace * sizeof(double));
	imJ21 = (double *) malloc(Nspace * sizeof(double));
	reJ22 = (double *) malloc(Nspace * sizeof(double));
	imJ22 = (double *) malloc(Nspace * sizeof(double));

	readJ20lambda(nspect, J20dag);
      } else {
	J20 = spectrum.J20[nspect];
	reJ21 = spectrum.reJ21[nspect];
	imJ21 = spectrum.imJ21[nspect];
	reJ22 = spectrum.reJ22[nspect];
	imJ22 = spectrum.imJ22[nspect];

	for (k = 0;  k < Nspace;  k++)
	  J20dag[k] = J20[k];
      }
      if (spectrum.updateJ)
	for (k = 0;  k < Nspace;  k++){
	  J20[k]   = 0.0;
	  reJ21[k] = 0.0;
	  imJ21[k] = 0.0;
	  reJ22[k] = 0.0;
	  imJ22[k] = 0.0;
	}

    }


    if (angle_dep) {
      for (mu = 0;  mu < Nrays;  mu++) {
	wmu = 0.5 * geometry.wmu[mu];
	if (input.backgr_pol) {

	  azi   = atan2(geometry.muy[mu],geometry.mux[mu]);
	  inc = geometry.muz[mu];                   /* inc is the same is mu*/

	  musq = SQ(geometry.muz[mu]);
	  threemu1 = TWOSQRTTWO * (3.0*musq - 1.0);
	  threemu2 = (3.0 * TWOSQRTTWO) * (musq - 1.0);



	}
	for (to_obs = 0;  to_obs <= 1;  to_obs++) {
	  initialize = (mu == 0 && to_obs == 0);

	  if (initialize || atmos.backgrflags[nspect].hasline)
	    readBackground(nspect, mu, to_obs);

	  if (initialize || boundbound)
	    Opacity(nspect, mu, to_obs, initialize);

	  if (eval_operator) addtoCoupling(nspect);
	  for (k = 0;  k < Nspace;  k++) {
	    chi[k] = as->chi[k] + as->chi_c[k];
	    S[k]   = as->eta[k] + as->eta_c[k] + as->sca_c[k]*Jdag[k];
	  }
	  if (solveStokes) {
	    for (k = Nspace;  k < 4*Nspace;  k++) Spol[0][k] = 0.0;

	    /* --- Add emissivity due to active set for Q, U, V -- ---- */

	    if (polarized_as) {
	      for (k = Nspace;  k < 4*Nspace;  k++)
		Spol[0][k] += as->eta[k];
	    }
	    /* --- Add emissivity due to background lines -- ---------- */

	    if (polarized_c) {
	      for (k = Nspace;  k < 4*Nspace;  k++)
		Spol[0][k] += as->eta_c[k];
	    }
	    /* --- Add emissivity due to background scattering -- ----- */

	    if (input.backgr_pol && input.StokesMode == FULL_STOKES) {
	      for (k = 0;  k < Nspace;  k++) {
		Spol[0][k] += threemu1 * as->sca_c[k]*J20dag[k];
		Spol[1][k] += threemu2 * as->sca_c[k]*J20dag[k];
	      }
	    }
	    for (n = 0;  n < 4;  n++) {
	      for (k = 0;  k < Nspace;  k++)
		Spol[n][k] /= chi[k];
	    }
	    ShortChar_Stokes(&geometry, nspect, mu, to_obs,
			     chi, Spol, Ipol, Psi);
	  } else {
	    for (k = 0;  k < Nspace;  k++)
	      S[k] /= chi[k];
	    ShortChar(&geometry, nspect, mu, to_obs, chi, S, I, Psi);
	  }
	  if (eval_operator) {
	    for (k = 0;  k < Nspace;  k++) Psi[k] /= chi[k];
	    addtoGamma(nspect, wmu, I, Psi);
	  }

	  if (spectrum.updateJ) {

	    /* --- Accumulate mean intensity and rates -- ----------- */


	    addtoRates(nspect, mu, to_obs, wmu, I, redistribute);

	    /* --- Accumulate anisotropy --            -------------- */

	    if (PRD_angle_dep) writeImu(nspect, mu, to_obs, I);
	  }	}
	/* --- start calculating spherical tensors --              -------------- */

	for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
	  atom = atmos.activeatoms[nact];

	  for (n = 0;  n < as->Nactiveatomrt[nact];  n++) {

	    line = as->art[nact][n].ptype.line;

	    if (!redistribute || line->PRD)
	      rate_lock = &line->rate_lock;

	    // Start calculating anisotropy here...

	    if (input.backgr_pol) {

	      if (input.Nthreads > 1)
		pthread_mutex_lock(rate_lock);

	      la = nspect - line->Nblue;

	      wlambda = getwlambda_line(line,la);

	      for (k = 0;  k < atmos.Nspace;  k++) {

		domg_dlam = wmu * line->phi[la][k] *line->wphi[k] * wlambda;

		J[k] += Ipol[0][k] * domg_dlam;

		J20[k] += (threemu1 * Ipol[0][k] + threemu2 * Ipol[1][k]) * domg_dlam;

		reJ21[k] += (sqrt(3.0)/2.0)*( sqrt(1.0 - musq) * ( -1.0 * inc * cos(azi) * (Ipol[0][k] + Ipol[1][k]) + sin(azi)*Ipol[2][k]) )   * domg_dlam;

		imJ21[k] += (sqrt(3.0)/2.0)*( sqrt(1.0 - musq) * ( -1.0 * inc * sin(azi) * (Ipol[0][k] + Ipol[1][k]) - cos(azi)*Ipol[2][k]) ) * domg_dlam;

		reJ22[k] += (sqrt(3.0)/4.0)*(cos(2.0*azi)*((1.0-musq) * Ipol[0][k] - (1.0 + musq)* Ipol[1][k]) + 2.0*sin(2.0*azi)*inc*Ipol[3][k])  * domg_dlam;

		imJ22[k] += (sqrt(3.0)/4.0)*(sin(2.0*azi)*((1.0-musq) * Ipol[0][k] - (1.0 + musq)* Ipol[1][k]) - 2.0*cos(2.0*azi)*inc*Ipol[3][k])  * domg_dlam;

	      } // k for loop
	    } //background polarization if statement
	  } //atom - radiative transition for loop
	} //atom for loop
      }
    }
  }  // nspect for loop

  FILE *fptr;

  fptr = fopen("J.txt","w");
  if(fptr == NULL)
    {
      printf("Error!");
      exit(1);
    }

  for (k = 0;  k < atmos.Nspace;  k++) {
    fprintf(fptr,"%0.20f \n",J[k]);
  }
  fclose(fptr);


  fptr = fopen("J20.txt","w");
  if(fptr == NULL)
    {
      printf("Error!");
      exit(1);
    }

  for (k = 0;  k < atmos.Nspace;  k++) {
    fprintf(fptr,"%0.20f \n",J20[k]);
  }
  fclose(fptr);


  fptr = fopen("reJ21.txt","w");
  if(fptr == NULL)
    {
      printf("Error!");
      exit(1);
    }

  for (k = 0;  k < atmos.Nspace;  k++) {
    fprintf(fptr,"%0.20f \n",reJ21[k]);
  }
  fclose(fptr);


  fptr = fopen("reJ22.txt","w");
  if(fptr == NULL)
    {
      printf("Error!");
      exit(1);
    }

  for (k = 0;  k < atmos.Nspace;  k++) {
    fprintf(fptr,"%0.40f \n",reJ22[k]);
  }
  fclose(fptr);

  fptr = fopen("imJ21.txt","w");
  if(fptr == NULL)
    {
      printf("Error!");
      exit(1);
    }

  for (k = 0;  k < atmos.Nspace;  k++) {
    fprintf(fptr,"%0.20f \n",imJ21[k]);
  }
  fclose(fptr);


  fptr = fopen("imJ22.txt","w");
  if(fptr == NULL)
    {
      printf("Error!");
      exit(1);
    }

  for (k = 0;  k < atmos.Nspace;  k++) {
    fprintf(fptr,"%0.20f \n",imJ22[k]);
  }
  fclose(fptr);
}
