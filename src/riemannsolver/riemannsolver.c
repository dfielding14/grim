#include "riemannsolver.h"

/* Returns a wavespeed */
REAL riemannSolver(const REAL fluxLeft[ARRAY_ARGS DOF],
                   const REAL fluxRight[ARRAY_ARGS DOF],
                   const REAL conservedVarsLeft[ARRAY_ARGS DOF],
                   const REAL conservedVarsRight[ARRAY_ARGS DOF],
                   const REAL primVarsLeft[ARRAY_ARGS DOF],
                   const REAL primVarsRight[ARRAY_ARGS DOF],
                   const struct geometry geom[ARRAY_ARGS DOF],
                   const int dir, REAL fluxes[ARRAY_ARGS DOF])
{
  REAL cMinLeft, cMaxLeft;
  struct fluidElement elem;
  setFluidElement(primVarsLeft, geom, &elem);
  waveSpeeds(&elem, geom, dir, &cMinLeft, &cMaxLeft);

  REAL cMinRight, cMaxRight;
  setFluidElement(primVarsRight, geom, &elem);
  waveSpeeds(&elem, geom, dir, &cMinRight, &cMaxRight);
  
	REAL cMax = fabs(fmax(fmax(0., cMaxLeft), cMaxRight));
	REAL cMin = fabs(fmax(fmax(0., -cMinLeft), -cMinRight));
	REAL cLaxFriedrichs = fmax(cMax, cMin);
    
  for (int var=0; var<DOF; var++) 
  {
    fluxes[var] = 0.5*(fluxLeft[var] + fluxRight[var]
                       - cLaxFriedrichs*(  conservedVarsRight[var]
                                         - conservedVarsLeft[var]
                                        )
                      );
                         
  }

  return cLaxFriedrichs;
}

void waveSpeeds(const struct fluidElement elem[ARRAY_ARGS 1],
                const struct geometry geom[ARRAY_ARGS 1],
                const int dir,
                REAL cMin[ARRAY_ARGS 1], REAL cMax[ARRAY_ARGS 1])
{
  REAL bCov[NDIM];
  conToCov(elem->bCon, geom, bCov);

  REAL bSqr = covDotCon(bCov, elem->bCon);

  REAL cAlvenSqr = bSqr/(bSqr + elem->primVars[RHO] 
                              + ADIABATIC_INDEX*elem->primVars[UU]);
  REAL csSqr = (ADIABATIC_INDEX)*(ADIABATIC_INDEX-1)*elem->primVars[UU]
              /(elem->primVars[RHO] + ADIABATIC_INDEX*elem->primVars[UU]);

  //Correction to characteristic speeds in fake EMHD model
  #if (FAKE_EMHD)
    REAL y = 2.*(ADIABATIC_INDEX-1)*elem->primVars[UU]/(bSqr+1.e-16);
    REAL FakeEMHDCoeff = 0.5*(-3.*y-2.+sqrt((3.*y+2)*(3.*y+2)+12*y)); 
    //Speed from Balbusaur, assuming dP = FakeEMHDCoeff*b^2/2
    //with the coefficient chose for the mirror saturation condition
    csSqr  = 3.*(ADIABATIC_INDEX)*(ADIABATIC_INDEX-1)*elem->primVars[UU]
      -FakeEMHDCoeff*bSqr*(ADIABATIC_INDEX-1);
    csSqr = csSqr/
      (3.*elem->primVars[RHO] + 3.*ADIABATIC_INDEX*elem->primVars[UU]-FakeEMHDCoeff*bSqr);
    if(csSqr>1.)
      csSqr=1.;
    
    cAlvenSqr = bSqr*(6.+3.*FakeEMHDCoeff)/
      (6.*elem->primVars[RHO]+6.*ADIABATIC_INDEX*elem->primVars[UU]
    +(6.+FakeEMHDCoeff)*bSqr);
    if(cAlvenSqr>1.)
      cAlvenSqr=1.;
  #endif


  REAL cVisSqr = 0.;
  #if (VISCOSITY)
    REAL beta = elem->tauVis/(2.*elem->eta);
    cVisSqr = 2./3./(elem->primVars[RHO] + ADIABATIC_INDEX*elem->primVars[UU])/beta;
  #endif

  REAL cConSqr = 0.;  
  #if (CONDUCTION)
    cConSqr = (ADIABATIC_INDEX-1.)*elem->kappa/elem->primVars[RHO]/elem->tau;  
  #endif
    
  REAL cmSqr = csSqr + cAlvenSqr - csSqr*cAlvenSqr + cVisSqr + cConSqr;
  if(cmSqr>1.)
    {
      //printf("cSqr = %e; csSqr=%e; cAlvenSqr=%e; cVisSqr=%e; cConSqr=%e\n",
      // 	     cmSqr,csSqr,cAlvenSqr,cVisSqr,cConSqr);
      //exit(1);
      cmSqr=1.;
    }
  
  REAL ACov[NDIM], ACon[NDIM];
  REAL BCov[NDIM], BCon[NDIM];
  for (int mu=0; mu<NDIM; mu++)
  {
    ACov[mu] = 0.;
    BCov[mu] = 0.;
  }
  ACov[dir] = 1.;
  BCov[0] = 1.;
  covToCon(ACov, geom, ACon);
  covToCon(BCov, geom, BCon);

  REAL ASqr = 0., BSqr=0., ADotU=0., BDotU=0., ADotB=0.; 
  for (int mu=0; mu<NDIM; mu++)
  {
    ASqr += ACov[mu]*ACon[mu];
    BSqr += BCov[mu]*BCon[mu];
    ADotU += ACov[mu]*elem->uCon[mu];
    BDotU += BCov[mu]*elem->uCon[mu];
    ADotB += ACov[mu]*BCon[mu];
  }

  REAL A = (BDotU*BDotU) - (BSqr + BDotU*BDotU)*cmSqr;
  REAL B = 2.*(ADotU*BDotU - (ADotB + ADotU*BDotU)*cmSqr);
  REAL C = ADotU*ADotU - (ASqr + ADotU*ADotU)*cmSqr;

  REAL discr = sqrt(B*B - 4.*A*C);

  REAL cPlus = -(-B + discr)/(2.*A);
  REAL cMinus = -(-B - discr)/(2.*A);

  if (cPlus > cMinus)
  {
    cMax[0] = cPlus;
    cMin[0] = cMinus;
  } else
  {
    cMax[0] = cMinus;
    cMin[0] = cPlus;
  }

}
