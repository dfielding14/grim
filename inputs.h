#ifndef GRIM_INPUT_H_
#define GRIM_INPUT_H_

#define COMPUTE_DIM         (2)
#define NDIM                (4)
#define N1                  (64)
#define N2                  (64)
#define NG                  (2)

/* Choose the geometry here:
 * KERRSCHILD or MINKOWSKI
*/
#define KERRSCHILD

/* Double or float precision? */
#define REAL                (double)

/* The entire global domain is divided into tiles that are small enough to fit
 * into the cache of the compute node or an accelerator. This technique
 * optimizes cache usage and makes explicit use of the deep memory hierarchies
 * which are prevalent in all kinds of processors and accelerators. 
 * 
 * Caution: 1) The tile sizes need to divide the global domain size in each
 *             direction exactly!
 *          2) We want a tile size that is as large as the cache. Too large a
 *             size and the code can fail silently because local memory in the
 *             OpenCL specification is a different amount for different 
 *             machines.
 *             NEED TO PUT A CHECK FOR THIS.
 *
*/
#define TILE_SIZE_X1        (8)
#define TILE_SIZE_X2        (8)

#ifdef KERRSCHILD
#define A_SPIN              (0.9375)    // Black hole spin
#define ADIABATIC_INDEX     (4./3.)

// Floor values:
#define SMALL               (1e-15)     // Minimum value of rho and u allowed when
                                        // computing fluxes inside the residual
                                        // evaluation function. Set to the limits of
                                        // REAL precision.

#define RHO_FLOOR_0         (1e-4)      // rho_floor = RHO_FLOOR_0 * radius^RHO_FLOOR_EXPONENT
#define RHO_FLOOR_EXPONENT  (-1.5)           
#define UU_FLOOR_0          (1e-6)      // u_floor = UU_FLOOR_0 * radius^UU_FLOOR_EXPONENT
#define UU_FLOOR_EXPONENT   (-2.5)      
#define GAMMA_MAX           (50.)       // Code cannot handle arbitrary large lorentz factors. 
                                        // Limit the Lorentz factor to this value.
// End of floor values

// Initial conditions:
// Parameters for Fishbone-Moncrief solution: 
// R_DISK_INNER_EDGE, R_PRESSURE_MAX, ENTROPY_CONSTANT
#define R_DISK_INNER_EDGE   (6.)        // Inner edge of the disk in the initial
                                        // conditions

#define R_PRESSURE_MAX      (12.)       // Pressure max of the disk in the initial
                                        // conditions

#define ENTROPY_CONSTANT    (1e-3)      // Constant entropy factor in polytropic 
                                        // EOS for the initial conditions
                                        // P = ENTROPY_CONSTANT rho^ADIABATIC_INDEX

#define BETA                (1e2)       // Plasma beta for the initial conditions
// End of initial conditions

/**
 * Boundaries for the computational domain X^mu = {t, X1, X2, phi}. The
 * physical domain is set by the transformation in XTox(). All
 * computational operations are performed in X^mu coordinates. The
 * discretization is uniform in X^mu coordinates and the solution for
 * all variables are in X^mu coordinates.
 *
 *   +----> X1
 *   |
 *   |
 *   v
 *   X2
 *      (X1_A, X2_A)     (X1_B, X2_B)
 *           +----------------+
 *           |                |
 *           |                |
 *           |                |
 *           |                | 
 *           +----------------+
 *      (X1_D, X2_D)     (X1_C, X2_C)
 *
*/
#define R_A     (.98*(1. + sqrt(1. - A_SPIN*A_SPIN))) 
#define X1_A    (log(R_A)) /* X1_A as a function of R_A must be consistent with
                              functional dependence set in xtoX() */
#define X2_A    (1e-10)

#define R_B     (40.)
#define X1_B    (log(R_B))
#define X2_B    (1e-10)

#define R_C     (40.)
#define X1_C    (log(R_C))
#define X2_C    (M_PI - 1e-10)

#define R_D     (.98*(1. + sqrt(1. - A_SPIN*A_SPIN)))
#define X1_D    (log(R_D))
#define X2_D    (M_PI - 1e-10)
/* End of boundaries for the computational domain.*/

/** Transformation parameters between x^mu and X^mu: 
 *  1) H_SLOPE: 1 <= H_SLOPE <= 0.
 *     More points near midplane as H_SLOPE -> 0
*/
#define H_SLOPE             (0.3)
/* End of transformation parameters. */


#endif /* KERRSCHILD */

/* Numerical differencing parameter for calculating the connection */
#define EPS (1e-5)

/* Variable mnemonics */
#define RHO (0)
#define UU (1)
#define U1 (2)
#define U2 (3)
#define U3 (4)
#define B1 (5)
#define B2 (6)
#define B3 (7)
#define DOF (8)

/* Boundary mnemonics. See boundary.h for description */
#define TILE_BOUNDARY       (99)
#define OUTFLOW             (100)
#define MIRROR              (101)
#define CONSTANT            (102)
#define PERIODIC            (103)
#define NONE                (104)


// iTile, jTile have ranges [-NG, TILE_SIZE+NG)
#define INDEX_LOCAL(iTile,jTile,var) (iTile+NG + \
                                      (TILE_SIZE_X1+2*NG)*(jTile+NG + \
                                      (TILE_SIZE_X2+2*NG)*(var)))

// i, j have ranges [0, N1), [0, N2)
#define INDEX_GLOBAL(i,j,var) (var + DOF*((i)+(N1)*(j)))

#define i_TO_X1_CENTER(i) (X1_START + (i + 0.5)*DX1)
#define j_TO_X2_CENTER(j) (X2_START + (j + 0.5)*DX2)
#define i_TO_X1_FACE(i) (X1_START + (i)*DX1)
#define j_TO_X2_FACE(j) (X2_START + (j)*DX2)


void gammaCalc(REAL* gamma,
               const REAL var[DOF],
               const REAL gcov[NDIM][NDIM])
{
    *gamma = 
        sqrt(1 + gcov[1][1]*var[U1]*var[U1] + 
                 gcov[2][2]*var[U2]*var[U2] + 
                 gcov[3][3]*var[U3]*var[U3] + 
              2*(gcov[1][2]*var[U1]*var[U2] + 
                 gcov[1][3]*var[U1]*var[U3] + 
                 gcov[2][3]*var[U2]*var[U3]));
}

void dgammaCalc_dt(REAL* dgamma_dt,
                   const REAL gamma,
                   const REAL var[DOF],
                   const REAL dvar_dt[DOF],
                   const REAL gcov[NDIM][NDIM])
{
    *dgamma_dt = 
        ((gcov[1][1]*var[U1]*dvar_dt[U1] + 
          gcov[2][2]*var[U2]*dvar_dt[U2] + 
          gcov[3][3]*var[U3]*dvar_dt[U3])+ 
         (gcov[1][2]*dvar_dt[U1]*var[U2] + 
          gcov[1][2]*var[U1]*dvar_dt[U2] + 
          gcov[1][3]*dvar_dt[U1]*var[U3] + 
          gcov[1][3]*var[U1]*dvar_dt[U3] + 
          gcov[2][3]*dvar_dt[U2]*var[U3] + 
          gcov[2][3]*var[U2]*dvar_dt[U3]))/gamma;
}

void uconCalc(REAL ucon[NDIM],
              const REAL gamma,
              const REAL alpha,
              const REAL var[DOF],
              const REAL gcon[NDIM][NDIM])
{
    ucon[0] = gamma/alpha;
    ucon[1] = var[U1] - gamma*gcon[0][1]*alpha;
    ucon[2] = var[U2] - gamma*gcon[0][2]*alpha;
    ucon[3] = var[U3] - gamma*gcon[0][3]*alpha;
}

void duconCalc_dt(REAL ducon_dt[NDIM],
                  const REAL dgamma_dt,
                  const REAL alpha,
                  const REAL dvar_dt[DOF],
                  const REAL gcon[NDIM][NDIM])
{
    ducon_dt[0] = dgamma_dt/alpha;
    ducon_dt[1] = dvar_dt[U1] - dgamma_dt*gcon[0][1]*alpha;
    ducon_dt[2] = dvar_dt[U2] - dgamma_dt*gcon[0][2]*alpha;
    ducon_dt[3] = dvar_dt[U3] - dgamma_dt*gcon[0][3]*alpha;
}

void covFromCon(REAL cov[NDIM],
                const REAL con[NDIM],
                const REAL gcov[NDIM][NDIM])
{
    cov[0] = gcov[0][0]*con[0] + gcov[0][1]*con[1] +
             gcov[0][2]*con[2] + gcov[0][3]*con[3];

    cov[1] = gcov[1][0]*con[0] + gcov[1][1]*con[1] +
             gcov[1][2]*con[2] + gcov[1][3]*con[3];

    cov[2] = gcov[2][0]*con[0] + gcov[2][1]*con[1] +
             gcov[2][2]*con[2] + gcov[2][3]*con[3];

    cov[3] = gcov[3][0]*con[0] + gcov[3][1]*con[1] +
             gcov[3][2]*con[2] + gcov[3][3]*con[3];
}

void conFromCov(REAL con[NDIM],
                const REAL cov[NDIM],
                const REAL gcon[NDIM][NDIM])
{
    con[0] = gcon[0][0]*cov[0] + gcon[0][1]*cov[1] +
             gcon[0][2]*cov[2] + gcon[0][3]*cov[3];

    con[1] = gcon[1][0]*cov[0] + gcon[1][1]*cov[1] +
             gcon[1][2]*cov[2] + gcon[1][3]*cov[3];

    con[2] = gcon[2][0]*cov[0] + gcon[2][1]*cov[1] +
             gcon[2][2]*cov[2] + gcon[2][3]*cov[3];

    con[3] = gcon[3][0]*cov[0] + gcon[3][1]*cov[1] +
             gcon[3][2]*cov[2] + gcon[3][3]*cov[3];

}

void conDotCov(REAL* ans,
               const REAL con[NDIM],
               const REAL cov[NDIM])
{
    *ans = con[0]*cov[0] + con[1]*cov[1] +
           con[2]*cov[2] + con[3]*cov[3];
}

void bconCalc(REAL bcon[NDIM],
              const REAL var[DOF],
              const REAL ucon[NDIM],
              const REAL ucov[NDIM])
{
    bcon[0] = var[B1]*ucov[1]+ var[B2]*ucov[2]+ var[B3]*ucov[3];
    
    bcon[1] = (var[B1] + bcon[0]*ucon[1])/ucon[0];
    bcon[2] = (var[B2] + bcon[0]*ucon[2])/ucon[0];
    bcon[3] = (var[B3] + bcon[0]*ucon[3])/ucon[0];
}

void dbconCalc_dt(REAL dbcon_dt[NDIM],
                  const REAL ucon[NDIM],
                  const REAL ducon_dt[NDIM],
                  const REAL ucov[NDIM],
                  const REAL ducov_dt[NDIM],
                  const REAL bcon[NDIM],
                  const REAL var[DOF],
                  const REAL dvar_dt[DOF])
{

    dbcon_dt[0] = 
        (dvar_dt[B1]*ucov[1] + dvar_dt[B2]*ucov[2] + dvar_dt[B3]*ucov[3] +
         var[B1]*ducov_dt[1] + var[B2]*ducov_dt[2] + var[B3]*ducov_dt[3]);

    dbcon_dt[1] =
        (-(var[B1] + bcon[0]*ucon[1])*ducon_dt[0]/(ucon[0]*ucon[0]) +
         (dvar_dt[B1] + bcon[0]*ducon_dt[1] + dbcon_dt[0]*ucon[1])/ucon[0]);

    dbcon_dt[2] =
        (-(var[B2] + bcon[0]*ucon[2])*ducon_dt[0]/(ucon[0]*ucon[0]) +
         (dvar_dt[B2] + bcon[0]*ducon_dt[2] + dbcon_dt[0]*ucon[2])/ucon[0]);

    dbcon_dt[3] =
        (-(var[B3] + bcon[0]*ucon[3])*ducon_dt[0]/(ucon[0]*ucon[0]) +
         (dvar_dt[B3] + bcon[0]*ducon_dt[3] + dbcon_dt[0]*ucon[3])/ucon[0]);
}

void bSqrCalc(REAL* bsqr,
              const REAL bcon[NDIM],
              const REAL bcov[NDIM])
{
    *bsqr = bcon[0]*bcov[0] + bcon[1]*bcov[1] +
            bcon[2]*bcov[2] + bcon[3]*bcov[3];
}

void mhdCalc(REAL mhd[NDIM][NDIM],
             const REAL var[DOF],
             const REAL ucon[NDIM],
             const REAL ucov[NDIM],
             const REAL bcon[NDIM],
             const REAL bcov[NDIM])
{
    REAL P = (ADIABATIC_INDEX - 1.)*var[UU];
    REAL bsqr;
    bSqrCalc(&bsqr, bcon, bcov);
    
#define DELTA(mu, nu) (mu==nu ? 1 : 0)

    for (int mu=0; mu<NDIM; mu++)
        for (int nu=0; nu<NDIM; nu++) {
            mhd[mu][nu] = (var[RHO] + var[UU] + P + bsqr)*ucon[mu]*ucov[nu] +
                          (P + 0.5*bsqr)*DELTA(mu, nu) - bcon[mu]*bcov[nu];
        }

#undef DELTA
}

void addSources(REAL dU_dt[DOF],
                const REAL ucon[NDIM],
                const REAL ucov[NDIM],
                const REAL bcon[NDIM],
                const REAL bcov[NDIM],
                const REAL gcon[NDIM][NDIM],
                const REAL gcov[NDIM][NDIM],
                const REAL mhd[NDIM][NDIM],
                const REAL var[DOF],
                const REAL g,
                const REAL X1, const REAL X2)
{
    
    REAL gcovh[NDIM][NDIM], gcovl[NDIM][NDIM];
    REAL conntmp[NDIM][NDIM][NDIM], conn[NDIM][NDIM][NDIM];
    REAL Xl[NDIM], Xh[NDIM];

    for (int k = 0; k < NDIM; k++) {
        Xl[0] = 0.; Xl[1] = X1; Xl[2] = X2; Xl[3] = 0.;
        Xh[0] = 0.; Xh[1] = X1; Xh[2] = X2; Xh[3] = 0.;
        Xl[k] = Xl[k] - EPS;
        Xh[k] = Xh[k] + EPS;
        gCovCalc(gcovh, Xh[1], Xh[2]);
        gCovCalc(gcovl, Xl[1], Xl[2]);
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				conn[i][j][k] = (gcovh[i][j] - gcovl[i][j])/(Xh[k] - Xl[k]);
            }
        }
    }

	/* now rearrange to find \Gamma_{ijk} */
	for (int i = 0; i < NDIM; i++)
		for (int j = 0; j < NDIM; j++)
			for (int k = 0; k < NDIM; k++)
			    conntmp[i][j][k] =
				    0.5 * (conn[j][i][k] + conn[k][i][j] -
					   conn[k][j][i]);

	/* finally, raise index */
	for (int i = 0; i < NDIM; i++)
		for (int j = 0; j < NDIM; j++)
			for (int k = 0; k < NDIM; k++) {
				conn[i][j][k] = 0.;
				for (int l = 0; l < NDIM; l++)
					conn[i][j][k] += gcon[i][l]*conntmp[l][j][k];
			}
    
    for (int j=0; j<NDIM; j++)
        for (int k=0; k<NDIM; k++) {
            dU_dt[UU] = dU_dt[UU] - g*(mhd[j][k]*conn[k][0][j]);
            dU_dt[U1] = dU_dt[U1] - g*(mhd[j][k]*conn[k][1][j]);
            dU_dt[U2] = dU_dt[U2] - g*(mhd[j][k]*conn[k][2][j]);
            dU_dt[U3] = dU_dt[U3] - g*(mhd[j][k]*conn[k][3][j]);

        }
}

void ComputeFluxAndU(REAL flux[DOF],
                     REAL U[DOF],
                     const REAL ucon[NDIM],
                     const REAL ucov[NDIM],
                     const REAL bcon[NDIM],
                     const REAL bcov[NDIM],
                     const REAL gcon[NDIM][NDIM],
                     const REAL gcov[NDIM][NDIM],
                     const REAL mhd[NDIM][NDIM],
                     const REAL var[DOF],
                     const REAL g,
                     const int dir)
{
    flux[RHO] = g*var[RHO]*ucon[dir];

    flux[UU] = g*mhd[dir][0];
    flux[U1] = g*mhd[dir][1];
    flux[U2] = g*mhd[dir][2];
    flux[U3] = g*mhd[dir][3];

    flux[B1] = g*(bcon[1]*ucon[dir] - bcon[dir]*ucon[1]);
    flux[B2] = g*(bcon[2]*ucon[dir] - bcon[dir]*ucon[2]);
    flux[B3] = g*(bcon[3]*ucon[dir] - bcon[dir]*ucon[3]);

    U[RHO] = g*var[RHO]*ucon[0];

    U[UU] = g*mhd[0][0];
    U[U1] = g*mhd[0][1];
    U[U2] = g*mhd[0][2];
    U[U3] = g*mhd[0][3];

    U[B1] = g*(bcon[1]*ucon[0] - bcon[0]*ucon[1]);
    U[B2] = g*(bcon[2]*ucon[0] - bcon[0]*ucon[2]);
    U[B3] = g*(bcon[3]*ucon[0] - bcon[0]*ucon[3]);

}

void ComputedU_dt(REAL dU_dt[DOF],
                  const REAL ucon[NDIM],
                  const REAL ducon_dt[NDIM],
                  const REAL ucov[NDIM],
                  const REAL ducov_dt[NDIM],
                  const REAL bcon[NDIM],
                  const REAL dbcon_dt[NDIM],
                  const REAL bcov[NDIM],
                  const REAL dbcov_dt[NDIM],
                  const REAL gcon[NDIM][NDIM],
                  const REAL gcov[NDIM][NDIM],
                  const REAL var[DOF],
                  const REAL dvar_dt[DOF],
                  const REAL gamma,
                  const REAL dgamma_dt,
                  const REAL alpha,
                  const REAL g)
{
    REAL P, dP_dt, bsqr, dbsqr_dt, tmp1, dtmp1_dt, tmp2, dtmp2_dt;

    bSqrCalc(&bsqr, bcon, bcov);

    dbsqr_dt = bcon[0]*dbcov_dt[0] + dbcon_dt[0]*bcov[0] +
               bcon[1]*dbcov_dt[1] + dbcon_dt[1]*bcov[1] +
               bcon[2]*dbcov_dt[2] + dbcon_dt[2]*bcov[2] +
               bcon[3]*dbcov_dt[3] + dbcon_dt[3]*bcov[3];

    P = (ADIABATIC_INDEX-1.)*var[UU];
    dP_dt = (ADIABATIC_INDEX-1.)*dvar_dt[UU];

    tmp1 = P + var[RHO] + var[UU] + bsqr;
    dtmp1_dt = dP_dt + dvar_dt[RHO] + dvar_dt[UU] + dbsqr_dt;

    tmp2 = P + 0.5*bsqr;
    dtmp2_dt = dP_dt + 0.5*dbsqr_dt;

    dU_dt[RHO] = g*(dvar_dt[RHO]*ucon[0] + var[RHO]*ducon_dt[0]);

    dU_dt[UU] = g*(dtmp1_dt*ucon[0]*ucov[0] +
                   tmp1*(ducon_dt[0]*ucov[0] + ucon[0]*ducov_dt[0]) +
                   dtmp2_dt - dbcon_dt[0]*bcov[0] - bcon[0]*dbcov_dt[0]);

    dU_dt[U1] = g*(dtmp1_dt*ucon[0]*ucov[1] +
                   tmp1*(ducon_dt[0]*ucov[1] + ucon[0]*ducov_dt[1]) -
                   dbcon_dt[0]*bcov[1] - bcon[0]*dbcov_dt[1]);

    dU_dt[U2] = g*(dtmp1_dt*ucon[0]*ucov[2] +
                   tmp1*(ducon_dt[0]*ucov[2] + ucon[0]*ducov_dt[2]) -
                   dbcon_dt[0]*bcov[2] - bcon[0]*dbcov_dt[2]);

    dU_dt[U3] = g*(dtmp1_dt*ucon[0]*ucov[3] +
                   tmp1*(ducon_dt[0]*ucov[3] + ucon[0]*ducov_dt[3]) -
                   dbcon_dt[0]*bcov[3] - bcon[0]*dbcov_dt[3]);

    dU_dt[B1] = g*dvar_dt[B1];
    dU_dt[B2] = g*dvar_dt[B2];
    dU_dt[B3] = g*dvar_dt[B3];

}
                  
void VChar(REAL* vmin, REAL* vmax,
           const REAL ucon[DOF], const REAL ucov[DOF],
           const REAL bsqr,
           const REAL gcon[NDIM][NDIM],
           const REAL var[DOF], const int dir)
{
    REAL Acov[NDIM], Acon[NDIM], Bcov[NDIM], Bcon[NDIM];
    REAL Asqr, Bsqr, vasqr, cssqr, cmsqr, Adotu, Bdotu, AdotB;
    REAL A, B, C, discr;
    REAL vp, vm;

    vasqr = bsqr/(bsqr + var[RHO] + ADIABATIC_INDEX*var[UU]);
    cssqr = (ADIABATIC_INDEX)*(ADIABATIC_INDEX-1)*var[UU]/(var[RHO] +
             ADIABATIC_INDEX*var[UU]);

    cmsqr = cssqr + vasqr - cssqr*vasqr;
    
    for (int mu=0; mu<NDIM; mu++)
        Acov[mu] = 0.;
    Acov[dir] = 1.;
    conFromCov(Acon, Acov, gcon);

    for (int mu=0; mu<NDIM; mu++)
        Bcov[mu] = 0.;
    Bcov[0] = 1.;
    conFromCov(Bcon, Bcov, gcon);

    conDotCov(&Asqr, Acon, Acov);
    conDotCov(&Bsqr, Bcon, Bcov);
    conDotCov(&Adotu, ucon, Acov);
    conDotCov(&Bdotu, ucon, Bcov);
    conDotCov(&AdotB, Acon, Bcov);

    A = (Bdotu*Bdotu) - (Bsqr + Bdotu*Bdotu)*cmsqr;
    B = 2.*(Adotu*Bdotu - (AdotB + Adotu*Bdotu)*cmsqr);
    C = Adotu*Adotu - (Asqr + Adotu*Adotu)*cmsqr;

    discr = sqrt(B*B - 4.*A*C);

    vp = -(-B + discr)/(2.*A);
    vm = -(-B - discr)/(2.*A);

    if (vp>vm) {
        *vmax = vp;
        *vmin = vm;
    } else {
        *vmax = vm;
        *vmin = vp;
    }
}

void setFloor(REAL var[DOF], 
              const REAL gamma,
              REAL X1, REAL X2)
{
    REAL r, theta;
    BLCoords(&r, &theta, X1, X2);

    REAL rhoFloor = RHO_MIN*pow(r, -1.5);
    REAL uFloor = U_MIN*pow(r, -2.5);

    if (rhoFloor < RHO_MIN_LIMIT)
        rhoFloor = RHO_MIN_LIMIT;

    if (uFloor < U_MIN_LIMIT)
        uFloor = U_MIN_LIMIT;
    
    if (var[RHO] < rhoFloor)
        var[RHO] = rhoFloor;

    if (var[UU] < uFloor)
        var[UU] = uFloor;

//    if (gamma > GAMMA_MAX) {
//        REAL f = sqrt((GAMMA_MAX*GAMMA_MAX - 1.) /
//				      (gamma * gamma - 1.));
//
//        var[U1] = f*var[U1];
//        var[U2] = f*var[U2];
//        var[U3] = f*var[U3];
//
//    }
}

void setFloorInit(REAL var[DOF], 
                  const REAL gamma,
                  REAL X1, REAL X2)
{
    REAL r, theta;
    BLCoords(&r, &theta, X1, X2);

    REAL rhoFloor = RHO_MIN_INIT*pow(r, -1.5);
    REAL uFloor = U_MIN_INIT*pow(r, -2.5);

    if (rhoFloor < RHO_MIN_LIMIT)
        rhoFloor = RHO_MIN_LIMIT;

    if (uFloor < U_MIN_LIMIT)
        uFloor = U_MIN_LIMIT;
    
    if (var[RHO] < rhoFloor)
        var[RHO] = rhoFloor;

    if (var[UU] < uFloor)
        var[UU] = uFloor;

    if (gamma > GAMMA_MAX) {
        REAL f = sqrt((GAMMA_MAX*GAMMA_MAX - 1.) /
				      (gamma * gamma - 1.));

        var[U1] = f*var[U1];
        var[U2] = f*var[U2];
        var[U3] = f*var[U3];

    }
}

#endif /* GRIM_INPUT_H_ */
