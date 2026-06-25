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

    //------------------- remove mean, then whiten to R = (2/3) I

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

    // lower Cholesky of C
    double L00=std::sqrt(C[0][0]);
    double L10=C[1][0]/L00, L11=std::sqrt(C[1][1]-L10*L10);
    double L20=C[2][0]/L00, L21=(C[2][1]-L20*L10)/L11,
           L22=std::sqrt(C[2][2]-L20*L20-L21*L21);
    double s=std::sqrt(2.0/3.0);                 // target component std (k_t = 1)
    for(int i=0;i<N;i++){
        double x0=u[i], x1=v[i], x2=w[i];
        double y0=x0/L00;
        double y1=(x1-L10*y0)/L11;
        double y2=(x2-L20*y0-L21*y1)/L22;
        u[i]=s*y0; v[i]=s*y1; w[i]=s*y2;
    }
}

////////////////////////////////////////////////////////////////////////////////
void domaincase_odt_homogeneousStrain::setCaseSpecificVars() {
    domn->rho->setVar();
    domn->dvisc->setVar();
}
