#include "geometry.h"

/* Function that sets the following quantities in X^mu coordinates in a zone:
 * 1) geom.gCov
 * 2) geom.gDet
 * 3) geom.gCon
 * 4) geom.alpha = 1./sqrt(-geom.gCon[0][0])
 * 
 * The functions to set the above quantities need to be called exactly as given
 * in order, since gDet needs gCov, gCon needs gDet and alpha needs gCon.
 *
 * Users call this function when the geometrical quantities are needed and not
 * the individual functions such as gCov and gCon since there is a certain order
 * that needs to be followed when the functions are called.
 *
 * @param Input: X, X^mu = {t, X1, X2, phi}, the coordinates where the
 *                  geometrical quantities are needed, REAL X[NDIM] array.
 * @param Output: geom, The geometry struct at those coordinates.
*/
void setGeometry(const REAL X[NDIM],
                 struct geometry geom)
{
    gCovFunc(X, geom.gCov);
    gDetFunc(geom.gCov, geom.gDet);
    gConFunc(geom.gCov, geom.gDet, geom.gCon);

    geom.alpha = 1./sqrt(-geom.gCon[0][0]);
}

/* The covariant metric in X^mu = {t, X1, X2, phi} coordinates. Inside the
 * function, the metric is first defined in the x^mu Kerr-Schild coordinates and
 * then transformed into the X^mu coordinates using:
 *
 * ds^2 = g_mu_nu dx^mu dx^nu
 *      = (g_mu_nu (dx^mu/dX^alpha) (dx^nu/dX^beta)) dX^alpha dX^beta
 *      = g_alpha_beta dX^alpha dX^beta
 *
 *  Ref: "A Measurement of the Electromagnetic Luminosity of a Kerr Black Hole"
 *       Jonathan C. McKinney and Charles F. Gammie
 *   
 * This function is not directly called by the user. All geometrical quantities
 * are set by calling setGeometry(). 
 *
 * @param Input: X, X^mu coordinates REAL X[NDIM] array
 * @param Output: gCov, The covariant metric, REAL gCov[NDIM][NDIM] array
*/
void gCovFunc(const REAL X[NDIM], REAL gCov[NDIM][NDIM])
{

#ifdef KERRSCHILD
    /* x^mu = {t, r, theta, phi}, X^mu = {t, X1, X2, phi} */

    REAL x[NDIM];

    Xtox(X, x);

    /* Easier to read with (r, theta) than (x[1], x[2]) */
    REAL r = x[1];
    REAL theta = x[2];

    /* r = exp(X1) => dr/dX = exp(X1) = r */
    REAL dr_dX1 = r;

    /* theta = pi*X2 + 0.5*(1 - H_SLOPE)*sin(2*pi*X2) 
       => dtheta/dX2 = pi + pi*(1 - H_SLOPE)*cos(2*pi*X2) */
    REAL dtheta_dX2 = M_PI + M_PI*(1 - H_SLOPE)*cos(2*M_PI*X[2]);

    REAL sigma = r*r + A_SPIN*A_SPIN*cos(theta)*cos(theta);
    REAL delta = r*r - 2*r + A_SPIN*A_SPIN;
    REAL A = pow(r*r + A_SPIN*A_SPIN, 2.0) - 
             A_SPIN*A_SPIN*delta*sin(theta)*sin(theta);

    /* -(1 - 2*r/sigma) dt^2 */
    gCov[0][0] = -(1. - 2.*r/sigma);     

    /* (4*r/sigma * dr/dX1) dt dX1 */
    gCov[0][1] = (2.*r/sigma) * dr_dX1; 

    /* (0) dt dX2 */
    gCov[0][2] = 0.; 

    /* -(4*a*r*sin(theta)^2/sigma) dt dphi */
    gCov[0][3] = -(2.*A_SPIN*r*sin(theta)*sin(theta)/sigma);

    /* (4*r/sigma * dr/dX1) dX1 dt */
    gCov[1][0] = geom.gCov[0][1];

    /* ( (1 + 2*r/sigma)*dr/dX1*dr/dX1) dX1 dX1 */
    gCov[1][1] = (1. + 2*r/sigma) * dr_dX1 * dr_dX1;
    
    /* (0) dX1 dX2 */
    gCov[1][2] = 0.;

    /* -(2*a*(1 + 2.*r/sigma)*sin(theta)^2*dr/dX1) dX1 dphi */
    gCov[1][3] = -A_SPIN*(1. + 2.*r/sigma)*sin(theta)*sin(theta)*dr_dX1;

    /* (0) dX2 dt */
    gCov[2][0] = geom.gCov[0][2];

    /* (0) dX2 dX1 */
    gCov[2][1] = geom.gCov[1][2];

    /* (sigma*dtheta/dX2*dtheta/dX2) dX2 dX2 */
    gCov[2][2] = sigma*dtheta_dX2*dtheta_dX2;

    /* (0) dX2 dphi */
    gCov[2][3] = 0.;

    /* -(4*a*r*sin(theta)^2/sigma) dphi dt */
    gCov[3][0] = geom.gCov[0][3];

    /* -(2*a*(1 + 2.*r/sigma)*sin(theta)^2*dr/dX1) dX1 dphi */
    gCov[3][1] = geom.gCov[1][3];

    /* (0) dphi dX2 */
    gCov[3][2] = geom.gCov[2][3];

    /* (sin(theta)^2*(sigma + a^2*(1. + 2*r/sigma)*sin(theta)^2) dphi dphi */
    gCov[3][3] = (sin(theta)*sin(theta)
                  *(sigma + A_SPIN*A_SPIN*(1. + 2*r/sigma) ) );
#endif /*KERRSCHILD*/

#ifdef MINKOWSKI
    gCov[0][0] = -1.;
    gCov[0][1] = 0.;
    gCov[0][2] = 0.;
    gCov[0][3] = 0.;

    gCov[1][0] = 0.;
    gCov[1][1] = 1.;
    gCov[1][2] = 0.;
    gCov[1][3] = 0.;

    gCov[2][0] = 0.;
    gCov[2][1] = 0.;
    gCov[2][2] = 1.;
    gCov[2][3] = 0.;

    gCov[3][0] = 0.;
    gCov[3][1] = 0.;
    gCov[3][2] = 0.;
    gCov[3][3] = 1.;

#endif /*MINKOWSKI*/
}

/* Determinant of the metric given a gCov. Uses analytic formula for the
 * determinant of a 4x4 matrix.
 * 
 * @param Input: gCov, The covariant metric gCov[NDIM][NDIM] array.
 * @param Output: gDet, The determinant of the metric, a single REAL number.
*/
void gDetFunc(const REAL gCov[NDIM][NDIM], REAL gDet)
{
    gDet = 
            gCov[0][0]*gCov[1][1]*gCov[2][2]*gCov[3][3] 
          + gCov[0][0]*gCov[1][2]*gCov[2][3]*gCov[3][1]
          + gCov[0][0]*gCov[1][3]*gCov[2][1]*gCov[3][2] 
          + gCov[0][1]*gCov[1][0]*gCov[2][3]*gCov[3][2] 
          + gCov[0][1]*gCov[1][2]*gCov[2][0]*gCov[3][3] 
          + gCov[0][1]*gCov[1][3]*gCov[2][2]*gCov[3][0] 
          + gCov[0][2]*gCov[1][0]*gCov[2][1]*gCov[3][3] 
          + gCov[0][2]*gCov[1][1]*gCov[2][3]*gCov[3][0] 
          + gCov[0][2]*gCov[1][3]*gCov[2][0]*gCov[3][1] 
          + gCov[0][3]*gCov[1][0]*gCov[2][2]*gCov[3][1]
          + gCov[0][3]*gCov[1][1]*gCov[2][0]*gCov[3][2] 
          + gCov[0][3]*gCov[1][2]*gCov[2][1]*gCov[3][0]
          - gCov[0][0]*gCov[1][1]*gCov[2][3]*gCov[3][2]
          - gCov[0][0]*gCov[1][2]*gCov[2][1]*gCov[3][3]
          - gCov[0][0]*gCov[1][3]*gCov[2][2]*gCov[3][1]
          - gCov[0][1]*gCov[1][0]*gCov[2][2]*gCov[3][3]
          - gCov[0][1]*gCov[1][2]*gCov[2][3]*gCov[3][0]
          - gCov[0][1]*gCov[1][3]*gCov[2][0]*gCov[3][2]
          - gCov[0][2]*gCov[1][0]*gCov[2][3]*gCov[3][1]
          - gCov[0][2]*gCov[1][1]*gCov[2][0]*gCov[3][3]
          - gCov[0][2]*gCov[1][3]*gCov[2][1]*gCov[3][0]
          - gCov[0][3]*gCov[1][0]*gCov[2][1]*gCov[3][2]
          - gCov[0][3]*gCov[1][1]*gCov[2][2]*gCov[3][0]
          - gCov[0][3]*gCov[1][2]*gCov[2][0]*gCov[3][1];
}

/* Contravariant metric given a covariant metric and it's determinant. Uses
 * analytic formula for the inverse of a 4x4 matrix.
 *
 * @param Input: gCov, The covariant metric: REAL gCov[NDIM][NDIM] array.
 * @param Input: gDet, The determinant of the metric, a single REAL number.
 * @param Output: gCon, The contravariant metric: REAL gCon[NDIM][NDIM] array.
*/
void gConFunc(const REAL gCov[NDIM][NDIM],
              const REAL gDet,
              REAL gCon[NDIM][NDIM])
{

    gCon[0][0] = 
        (  gCov[1][1]*gCov[2][2]*gCov[3][3]
         + gCov[1][2]*gCov[2][3]*gCov[3][1]
         + gCov[1][3]*gCov[2][1]*gCov[3][2]
         - gCov[1][1]*gCov[2][3]*gCov[3][2]
         - gCov[1][2]*gCov[2][1]*gCov[3][3]
         - gCov[1][3]*gCov[2][2]*gCov[3][1])/gDet;

    gCon[0][1] = 
        (  gCov[0][1]*gCov[2][3]*gCov[3][2]
         + gCov[0][2]*gCov[2][1]*gCov[3][3]
         + gCov[0][3]*gCov[2][2]*gCov[3][1]
         - gCov[0][1]*gCov[2][2]*gCov[3][3]
         - gCov[0][2]*gCov[2][3]*gCov[3][1] 
         - gCov[0][3]*gCov[2][1]*gCov[3][2])/gDet;

    gCon[0][2] = 
        (  gCov[0][1]*gCov[1][2]*gCov[3][3]
         + gCov[0][2]*gCov[1][3]*gCov[3][1]
         + gCov[0][3]*gCov[1][1]*gCov[3][2]
         - gCov[0][1]*gCov[1][3]*gCov[3][2]
         - gCov[0][2]*gCov[1][1]*gCov[3][3]
         - gCov[0][3]*gCov[1][2]*gCov[3][1])/gDet;

    gCon[0][3] = 
        (  gCov[0][1]*gCov[1][3]*gCov[2][2]
         + gCov[0][2]*gCov[1][1]*gCov[2][3]
         + gCov[0][3]*gCov[1][2]*gCov[2][1]
         - gCov[0][1]*gCov[1][2]*gCov[2][3]
         - gCov[0][2]*gCov[1][3]*gCov[2][1]
         - gCov[0][3]*gCov[1][1]*gCov[2][2])/gDet;

    gCon[1][0] = gCon[0][1];
    
    gCon[1][1] = 
        (  gCov[0][0]*gCov[2][2]*gCov[3][3]
         + gCov[0][2]*gCov[2][3]*gCov[3][0]
         + gCov[0][3]*gCov[2][0]*gCov[3][2]
         - gCov[0][0]*gCov[2][3]*gCov[3][2]
         - gCov[0][2]*gCov[2][0]*gCov[3][3]
         - gCov[0][3]*gCov[2][2]*gCov[3][0])/gDet;

    gCon[1][2] = 
        (  gCov[0][0]*gCov[1][3]*gCov[3][2]
         + gCov[0][2]*gCov[1][0]*gCov[3][3]
         + gCov[0][3]*gCov[1][2]*gCov[3][0]
         - gCov[0][0]*gCov[1][2]*gCov[3][3]
         - gCov[0][2]*gCov[1][3]*gCov[3][0]
         - gCov[0][3]*gCov[1][0]*gCov[3][2])/gDet;

    gCon[1][3] = 
        (  gCov[0][0]*gCov[1][2]*gCov[2][3]
         + gCov[0][2]*gCov[1][3]*gCov[2][0]
         + gCov[0][3]*gCov[1][0]*gCov[2][2]
         - gCov[0][0]*gCov[1][3]*gCov[2][2]
         - gCov[0][2]*gCov[1][0]*gCov[2][3]
         - gCov[0][3]*gCov[1][2]*gCov[2][0])/gDet;

    gCon[2][0] = gCon[0][2];
    gCon[2][1] = gCon[1][2];

    gCon[2][2] =
        (  gCov[0][0]*gCov[1][1]*gCov[3][3]
         + gCov[0][1]*gCov[1][3]*gCov[3][0]
         + gCov[0][3]*gCov[1][0]*gCov[3][1]
         - gCov[0][0]*gCov[1][3]*gCov[3][1]
         - gCov[0][1]*gCov[1][0]*gCov[3][3]
         - gCov[0][3]*gCov[1][1]*gCov[3][0])/gDet;

    gCon[2][3] =
        (  gCov[0][0]*gCov[1][3]*gCov[2][1]
         + gCov[0][1]*gCov[1][0]*gCov[2][3]
         + gCov[0][3]*gCov[1][1]*gCov[2][0]
         - gCov[0][0]*gCov[1][1]*gCov[2][3]
         - gCov[0][1]*gCov[1][3]*gCov[2][0]
         - gCov[0][3]*gCov[1][0]*gCov[2][1])/gDet;

    gCon[3][0] = gCon[0][3];
    gCon[3][1] = gCon[1][3];
    gCon[3][2] = gCon[2][3];

    gCon[3][3] =
        (  gCov[0][0]*gCov[1][1]*gCov[2][2] + 
         + gCov[0][1]*gCov[1][2]*gCov[2][0] + 
         + gCov[0][2]*gCov[1][0]*gCov[2][1] - 
         - gCov[0][0]*gCov[1][2]*gCov[2][1] - 
         - gCov[0][1]*gCov[1][0]*gCov[2][2] - 
         - gCov[0][2]*gCov[1][1]*gCov[2][0])/gDet;
}

/* Transform covariant vector to a contravariant vector.
 *
 * @param Input: vecCov, The covariant vector, REAL vecCov[NDIM] array.
 * @param Input: geom, The geometry struct.
 * @param Output: vecCon, The contravariant vector, REAL vecCon[NDIM] array.
*/
void covToCon(const REAL vecCov[NDIM],
              const struct geometry geom,
              REAL vecCon[NDIM])
{
    for (int mu=0; mu<NDIM; mu++)
    {
        vecCon[mu] = 0.;
        for (int nu=0; nu<NDIM; nu++)
        {
            vecCon[mu] += geom.gCov[mu][nu]*vecCov[nu];
        }
    }
}

/* Transform contravariant vector to a covariant vector.
 *
 * @param Input: vecCon, The contravariant vector, REAL vecCon[NDIM] array.
 * @param Input: geom, The geometry struct.
 * @param Output: vecCov, The covariant vector, REAL vecCov[NDIM] array.
*/
void conToCov(const REAL vecCon[NDIM],
              const struct geometry geom,
              REAL vecCov[NDIM])
{
    for (int mu=0; mu<NDIM; mu++)
    {
        vecCov[mu] = 0.;
        for (int nu=0; nu<NDIM; nu++)
        {
            vecCov[mu] += geom.gCon[mu][nu]*vecCon[nu];
        }
    }
}

/* Dot product between a covariant and a contravariant vector.
 *
 * @param Input: vecCov, The covariant vector, REAL vecCov[NDIM] array.
 * @param Input: vecCon, The contravariant vector, REAL vecCon[NDIM] array.
 * @return Output: dot, The dot product, a single REAL number.
*/
REAL covDotCon(const REAL vecCov[NDIM],
               const REAL vecCon[NDIM])
{
    REAL ans = 0.;
    for (int mu=0; mu<NDIM; mu++)
    {
        ans += vecCov[mu]*vecCon[mu];
    }

    return ans;
}

/* Christoffel symbols of the first kind: Gamma_eta_mu_nu in a coordinate
 * basis.
 *
 * Gamma_eta_mu_nu = 0.5*(  d(g_eta_mu)/dX^nu 
 *                        + d(g_eta_nu)/dX^mu 
 *                        - d(g_mu_nu)/dX^eta)
 *
 * Implemented numerically as follows:
 *
 * We need to compute Gamma_eta_mu_nu at the given X coordinates X^alpha = {t,
 * X1, X2, phi}. The indices {mu, nu, eta} are each one of {0, 1, 2, 3}. 
 * 
 * dX^nu = (X^alpha[nu] + EPS) - (X^alpha[nu] - EPS) = 2.*EPS
 * Similarly, dX^mu = dX^eta = 2.*EPS
 * 
 * Therefore:
 *
 * Gamma_eta_mu_nu = 0.5
 *                   *(  (  gCov(XNuPlusEpsilon)[eta,mu]
 *                        - gCov(XNuMinusEpsilon)[eta,mu])/(2.*EPS)
 *                     + (  gCov(XMuPlusEpsilon)[eta,nu] 
 *                        - gCov(XMuMinusEpsilon)[eta,nu])/(2.*EPS)
 *                     - (  gCov(XEtaPlusEpsilon)[mu,nu] 
 *                        - gCov(XEtaMinusEpsilon)[mu,nu])/(2.*EPS))
 *
 * where:
 *      XNuPlusEpsilon = X^alpha[nu] + EPS 
 *      XNuMinusEpsilon = X^alpha[nu] - EPS 
 *      (all components other than nu are equal to X^alpha)
 *
 *      XMuPlusEpsilon = X^alpha[mu] + EPS
 *      XMuMinusEpsilon = X^alpha[mu] - EPS
 *      (all components other than mu are equal to X^alpha)
 *
 *      XEtaPlusEpsilon = X^alpha[eta] + EPS
 *      XEtaMinusEpsilon = X^alpha[eta] - EPS
 *      (all components other than eta are equal to X^alpha)
 *
 *      EPS is a small parameter for numerical differentiation and is defined in
 *      inputs.h
 *
 * The input coordinates for this function are typically located at the center
 * of a zone since the connection source terms need to located at the zone
 * center.
 *
 * @param Input: eta, integer belonging to [0, NDIM).
 * @param Input: mu, integer belonging to [0, NDIM).
 * @param Input: nu, integer belonging to [0, NDIM).
 * @param Input: X, X^alpha = {t, X1, X2, phi} the coordinates where
 *                  Gamma_eta_mu_nu is to be computed, REAL X[NDIM] array.
 * @return Output: Gamma_eta_mu_nu, a single REAL number which gives the
 *                 Christoffel symbols of the first kind for the specific
 *                 indices {mu, nu, eta}.
*/
REAL gammaDownDownDown(const int eta, const int mu, const int nu,
                       const REAL X[NDIM])
{
    /* Handle the three differentiations seperately to save storage space */
    REAL XEpsilon[NDIM];
    REAL gCovEpsilon[NDIM][NDIM];

    REAL ans = 0;

    /* First, take care of d(g_eta_mu)/dX^nu. To compute this numerically we
     * first do XEpsilon^alpha[nu] = X^alpha[nu] + EPS and then compute the
     * metric corresponding to XEpsilon^alpha coordinates and add this to *ans
     * Now we do XEpsilon^alpha[nu] = X^alpha[nu] - EPS and then compute the
     * metric corresponding to XEpsilon^alpha coordinates and then subtract it
     * from *ans. And this entire thing is divided by 2.*EPS, corresponding to a
     * centered difference around X^alpha. By doing the computations seperately
     * for +EPS and -EPS, we cut storage requirements by half.*/

    /* Handle +EPS first */
    for (int alpha=0; alpha<NDIM; alpha++)
    {
        XEpsilon[alpha] = X[alpha];
    }
    XEpsilon[nu] += EPS;
    gCovFunc(XEpsilon, gCovEpsilon);
    ans += 0.5*gCovEpsilon[eta][mu]/(2.*EPS);

    /* Now do -EPS */
    for (int alpha=0; alpha<NDIM; alpha++)
    {
        XEpsilon[alpha] = X[alpha];
    }
    XEpsilon[nu] -= EPS;
    gCovFunc(XEpsilon, gCovEpsilon);
    ans -= 0.5*gCovEpsilon[eta][mu]/(2.*EPS);
    /* End of d(g_eta_mu)/dX^nu */



    /* Now, d(g_eta_nu)/dX^mu */
    /* +EPS first */
    for (int alpha=0; alpha<NDIM; alpha++)
    {
        XEpsilon[alpha] = X[alpha];
    }
    XEpsilon[mu] += EPS;
    gCovFunc(XEpsilon, gCovEpsilon);
    ans += 0.5*(gCovEpsilon[eta][nu])/(2.*EPS);

    /* Now do -EPS */
    for (int alpha=0; alpha<NDIM; alpha++)
    {
        XEpsilon[alpha] = X[alpha];
    }
    XEpsilon[mu] -= EPS;
    gCovFunc(XEpsilon, gCovEpsilon);
    ans -= 0.5*gCovEpsilon[eta][nu]/(2.*EPS);
    /* End of d(g_eta_nu)/dX^mu */

    /* Finally, d(g_mu_nu)/dX^eta */
    /* +EPS first */
    for (int alpha=0; alpha<NDIM; alpha++)
    {
        XEpsilon[alpha] = X[alpha];
    }
    XEpsilon[eta] += EPS;
    gCovFunc(XEpsilon, gCovEpsilon);
    ans += 0.5*(gCovEpsilon[mu][nu])/(2.*EPS);

    /* Now do -EPS */
    for (int alpha=0; alpha<NDIM; alpha++)
    {
        XEpsilon[alpha] = X[alpha];
    }
    XEpsilon[eta] -= EPS;
    gCovFunc(XEpsilon, gCovEpsilon);
    ans -= 0.5*gCovEpsilon[mu][nu]/(2.*EPS);
    /* End of d(g_mu_nu)/dX^eta */

    return ans;
}

/*  Returns sqrt(-gDet) * T^kappa_lamda * Gamma^lamda_nu_kappa for a given nu.
 * 
 *  = sqrt(-gDet)*T^kappa_lamda * g^lamda^mu * Gamma_mu_nu_kappa
 *
*/

REAL connectionTerm(const REAL TUpDown[NDIM][NDIM],
                    const struct geometry geom,
                    const int nu, const REAL X[NDIM])
{
    REAL ans = 0., g = geom.sqrt(-gDet);

    for (int kappa=0; kappa<NDIM; kappa++)
    {
        for (int lamda=0; lamda<NDIM; lamda++)
        {
            for (int mu=0; mu<NDIM; mu++)
            {
                ans +=   g * TUpDown[kappa][lamda] * geom.gCon[lamda][mu]
                       * gammaDownDownDown(mu, nu, kappa, X);
            }
        }
    }

    return ans;
}
