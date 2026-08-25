/**
 * @file domaincase_odt_homogeneousStrain.cc
 * @brief Source for class domaincase_odt_homogeneousStrain
 */

#include "domaincase_odt_homogeneousStrain.h"
#include "domain.h"
#include "dv.h"
#include "dv_pos.h"
#include "dv_posf.h"
#include "dv_rho_const.h"
#include "dv_dvisc_const.h"
#include "dv_uvw.h"
#include <random>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/** Initialization: register variables and seed a band-limited, isotropic
 *  fluctuation field from a prescribed 1D spectrum, then whiten it so that
 *  R_ij(0) = (2/3) delta_ij  (k_t = 1), matching Level 0. Replacing the former
 *  white-noise field removes the grid-scale energy that mesh interpolation
 *  damped, and provides a resolved spectrum for the Level 1b distortion study.
 */
void domaincase_odt_homogeneousStrain::init(domain *p_domn){

    domn = p_domn;

    domn->v.push_back(new dv_pos(        domn, "pos",   false, true));
    domn->v.push_back(new dv_posf(       domn, "posf",  false, true));
    domn->v.push_back(new dv_rho_const(  domn, "rho",   false, false));
    domn->v.push_back(new dv_dvisc_const(domn, "dvisc", false, false));
    domn->v.push_back(new dv_uvw(        domn, "uvel",  true,  true));
    domn->v.push_back(new dv_uvw(        domn, "vvel",  true,  true));
    domn->v.push_back(new dv_uvw(        domn, "wvel",  true,  true));

    domn->pos   = domn->v.at(0);
    domn->posf  = domn->v.at(1);
    domn->rho   = domn->v.at(2);
    domn->dvisc = domn->v.at(3);
    domn->uvel  = domn->v.at(4);
    domn->vvel  = domn->v.at(5);
    domn->wvel  = domn->v.at(6);

    //------------------- mesh adaption variables (all three components)

    vector<dv*> phi;
    phi.push_back(domn->uvel);
    phi.push_back(domn->vvel);
    phi.push_back(domn->wvel);
    domn->mesher->init(domn, phi);

    //------------------- seed an isotropic Gaussian field

    int N = domn->ngrd;
    std::mt19937 rng(domn->pram->seed >= 0 ? domn->pram->seed : 22);

    vector<double> &u = domn->uvel->d;
    vector<double> &v = domn->vvel->d;
    vector<double> &w = domn->wvel->d;

    //------------------- band-limited isotropic field from a prescribed
    //  1D spectrum  E(k) = (k/kp)^4 exp(-2 (k/kp)^2)  (Passot-Pouquet form:
    //  compact, peaked at kp, negligible energy near the grid scale, so the
    //  field is well resolved and not damped by mesh interpolation). Each
    //  component is an independent random-phase Fourier sum, giving isotropy.
    //  The grid is uniform at initialization (dv_posf builds it uniform), so
    //  cell centres are known analytically: y_i = xDomainCenter - L/2 + (i+1/2) dx.

    const double L  = domn->pram->domainLength;
    const double x0 = domn->pram->xDomainCenter - 0.5*L;
    const double dx = L / N;
    const double dk = 2.0*M_PI / L;                          // fundamental wavenumber
    const double kp = 2.0*M_PI*domn->pram->specKpWaves / L;  // spectral peak
    const int    Nm = domn->pram->specNmodes;                // number of modes

    vector<double> amp(Nm+1, 0.0);                           // amplitude ~ sqrt(E(k_n))
    for(int n=1;n<=Nm;n++){
        double r = (n*dk)/kp;
        amp[n] = std::sqrt( std::pow(r,4.0)*std::exp(-2.0*r*r) );
    }
    std::uniform_real_distribution<double> uni(0.0, 2.0*M_PI);
    vector<double>* comp[3] = {&u, &v, &w};
    for(int c=0;c<3;c++){
        vector<double> ph(Nm+1);
        for(int n=1;n<=Nm;n++) ph[n] = uni(rng);            // independent phases -> isotropy
        vector<double> &f = *comp[c];
        for(int i=0;i<N;i++){
            double y = x0 + (i+0.5)*dx;
            double s = 0.0;
            for(int n=1;n<=Nm;n++) s += amp[n]*std::cos(n*dk*y + ph[n]);
            f[i] = s;
        }
    }

    //------------------- remove mean, then whiten SYMMETRICALLY to R = (2/3) I
    //  y = C^{-1/2} x treats the three components identically. The previous
    //  ordered Cholesky whitening (u kept pure, v corrected against u, w
    //  against both) imprinted a component-ordered spectral bias that the
    //  eddy dynamics amplified to b_33 ~ +0.015 in nominal HIT (1024-rlz
    //  null-test ensembles + order-reversal attribution, 2026-08-25); the
    //  symmetric root removes the ordering by construction.

    double mu=0, mv=0, mw=0;
    for(int i=0;i<N;i++){ mu+=u[i]; mv+=v[i]; mw+=w[i]; }
    mu/=N; mv/=N; mw/=N;
    for(int i=0;i<N;i++){ u[i]-=mu; v[i]-=mv; w[i]-=mw; }

    double C[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    for(int i=0;i<N;i++){
        double f[3]={u[i],v[i],w[i]};
        for(int a=0;a<3;a++) for(int b=0;b<3;b++) C[a][b]+=f[a]*f[b];
    }
    for(int a=0;a<3;a++) for(int b=0;b<3;b++) C[a][b]/=N;

    // Jacobi eigen-decomposition C = Q diag(lam) Q^T (3x3 symmetric)
    double A[3][3], Q[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    for(int a=0;a<3;a++) for(int b=0;b<3;b++) A[a][b]=C[a][b];
    for(int sweep=0; sweep<60; sweep++){
        int p=0, q=1; double mx=std::abs(A[0][1]);
        if(std::abs(A[0][2])>mx){ mx=std::abs(A[0][2]); p=0; q=2; }
        if(std::abs(A[1][2])>mx){ mx=std::abs(A[1][2]); p=1; q=2; }
        if(mx < 1.0e-15) break;
        double phi = 0.5*std::atan2(2.0*A[p][q], A[q][q]-A[p][p]);
        double c = std::cos(phi), sn = std::sin(phi);
        for(int k=0;k<3;k++){ double akp=A[k][p], akq=A[k][q];
            A[k][p]=c*akp-sn*akq; A[k][q]=sn*akp+c*akq; }
        for(int k=0;k<3;k++){ double apk=A[p][k], aqk=A[q][k];
            A[p][k]=c*apk-sn*aqk; A[q][k]=sn*apk+c*aqk; }
        for(int k=0;k<3;k++){ double qkp=Q[k][p], qkq=Q[k][q];
            Q[k][p]=c*qkp-sn*qkq; Q[k][q]=sn*qkp+c*qkq; }
    }
    // W = Q diag(1/sqrt(lam)) Q^T  (symmetric inverse square root)
    double W[3][3];
    for(int a=0;a<3;a++) for(int b=0;b<3;b++){
        W[a][b]=0.0;
        for(int k=0;k<3;k++) W[a][b] += Q[a][k]*Q[b][k]/std::sqrt(A[k][k]);
    }
    double s=std::sqrt(2.0/3.0);                 // target component std (k_t = 1)
    for(int i=0;i<N;i++){
        double x0=u[i], x1=v[i], x2=w[i];
        u[i]=s*(W[0][0]*x0 + W[0][1]*x1 + W[0][2]*x2);
        v[i]=s*(W[1][0]*x0 + W[1][1]*x1 + W[1][2]*x2);
        w[i]=s*(W[2][0]*x0 + W[2][1]*x1 + W[2][2]*x2);
    }
}

////////////////////////////////////////////////////////////////////////////////
void domaincase_odt_homogeneousStrain::setCaseSpecificVars() {
    domn->rho->setVar();
    domn->dvisc->setVar();
}
