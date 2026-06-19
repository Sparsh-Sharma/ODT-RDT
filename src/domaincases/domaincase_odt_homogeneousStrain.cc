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
/** Initialization: register variables and seed an isotropic fluctuation field
 *  whitened so that R_ij(0) = (2/3) delta_ij  (k_t = 1), matching Level 0.
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
    std::normal_distribution<double> gss(0.0, 1.0);

    vector<double> &u = domn->uvel->d;
    vector<double> &v = domn->vvel->d;
    vector<double> &w = domn->wvel->d;
    for(int i=0;i<N;i++){ u[i]=gss(rng); v[i]=gss(rng); w[i]=gss(rng); }

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
