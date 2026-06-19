// strain_harness.cpp
// -----------------------------------------------------------------------------
// Level 1a logic check for the strain-coupled ODT formulation.
//
// Mirrors the C++ ODT code path that we will modify, but on a standalone uniform
// periodic "line" with eddies OFF and viscosity ZERO, so the continuous phase is
//      du_i/dt = A_cal_ij u_j ,   A_cal = -A + B,
// exactly the rapid-distortion limit integrated in the Level 0 Python script.
// Components R_ij are the line average of u_i u_j; B is the rapid pressure-strain
// operator from the Lyapunov solve  B R + R B = Pi^r,  with Pi^r from the chosen
// closure (IP or LRR).  This reproduces the *model* curves of Level 0 (NOT the
// exact-RDT curve, which is a separate benchmark).
//
// Routines named to match the real code: updateStrainOperator(), getRhsSrc(),
// advanceExplicit().  Build:  g++ -O2 -o strain_harness strain_harness.cpp
// -----------------------------------------------------------------------------
#include <vector>
#include <array>
#include <random>
#include <cmath>
#include <cstdio>
#include <string>

using namespace std;
using M3 = array<array<double,3>,3>;

static const double TWO3 = 2.0/3.0;

// ---- 3x3 helpers -----------------------------------------------------------
M3 zero3(){ M3 m{}; for(int i=0;i<3;i++)for(int j=0;j<3;j++)m[i][j]=0; return m; }
double trace3(const M3&a){ return a[0][0]+a[1][1]+a[2][2]; }

// ---- closures: rapid pressure-strain Pi^r given A, R -----------------------
M3 production(const M3&A, const M3&R){            // P_ij = -(A R + R A^T)
    M3 P=zero3();
    for(int i=0;i<3;i++)for(int j=0;j<3;j++){
        double s=0; for(int k=0;k<3;k++) s += A[i][k]*R[k][j] + R[i][k]*A[j][k];
        P[i][j] = -s;
    }
    return P;
}
M3 rapid_IP(const M3&A, const M3&R){              // -C2 (P - 1/3 trP I), C2=3/5
    M3 P=production(A,R); double tr=trace3(P); M3 Pir=zero3();
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)
        Pir[i][j] = -0.6*(P[i][j] - (i==j? tr/3.0:0.0));
    return Pir;
}
M3 rapid_LRR(const M3&A, const M3&R){             // LRR-QI, C2=4/5,C3=7/4,C4=131/100
    double kt = 0.5*trace3(R);
    M3 S=zero3(), W=zero3(), b=zero3();
    for(int i=0;i<3;i++)for(int j=0;j<3;j++){
        S[i][j]=0.5*(A[i][j]+A[j][i]); W[i][j]=0.5*(A[i][j]-A[j][i]);
        b[i][j]=R[i][j]/(2*kt) - (i==j?1.0/3.0:0.0);
    }
    double trbS=0; for(int i=0;i<3;i++)for(int k=0;k<3;k++) trbS += b[i][k]*S[i][k];
    M3 Pir=zero3();
    for(int i=0;i<3;i++)for(int j=0;j<3;j++){
        double bS_Sb=0, Wb_bW=0;
        for(int k=0;k<3;k++){ bS_Sb += b[i][k]*S[j][k]+S[i][k]*b[j][k];   // b_ik S_jk + S_ik b_jk
                              Wb_bW += W[i][k]*b[j][k] - b[i][k]*W[j][k]; }// (Wb - bW)_ij form
        Pir[i][j] = 0.8*kt*S[i][j]
                  + 1.75*kt*(bS_Sb - (i==j? TWO3*trbS:0.0))
                  + 1.31*kt*Wb_bW;
    }
    return Pir;
}

// ---- Lyapunov solve  B R + R B = Pi^r  for symmetric B ---------------------
// General case: diagonalise R = Q diag(l) Q^T, then B = Q M Q^T with
// M_ij = (Q^T Pi^r Q)_ij / (l_i + l_j).  Here R stays diagonal (diagonal A),
// so this reduces to B_ij = Pi^r_ij / (R_ii + R_jj); we use that directly and
// assert the off-diagonal R is negligible.
M3 lyapunovB(const M3&R, const M3&Pir){
    M3 B=zero3();
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)
        B[i][j] = Pir[i][j]/(R[i][i]+R[j][j]);
    return B;
}

// ---- the "domain": a uniform periodic line of velocity vectors -------------
struct Domain {
    int N;
    vector<array<double,3>> u;       // u[p] = (u1,u2,u3) at point p
    M3 Acal;                         // combined strain operator -A + B
    M3 A;                            // imposed mean velocity gradient
    string closure;
    vector<array<double,3>> rhsSrc;  // source term per point (mirrors dv_uvw)

    M3 reynoldsStress(){             // line average of u_i u_j (mode: line average)
        M3 R=zero3();
        for(int p=0;p<N;p++)for(int i=0;i<3;i++)for(int j=0;j<3;j++) R[i][j]+=u[p][i]*u[p][j];
        for(int i=0;i<3;i++)for(int j=0;j<3;j++) R[i][j]/=N;
        return R;
    }
    void updateStrainOperator(){     // recompute B(R) and Acal each substep
        M3 R = reynoldsStress();
        M3 Pir = (closure=="IP") ? rapid_IP(A,R) : rapid_LRR(A,R);
        M3 B = lyapunovB(R, Pir);
        for(int i=0;i<3;i++)for(int j=0;j<3;j++) Acal[i][j] = -A[i][j] + B[i][j];
    }
    void getRhsSrc(){                // rhsSrc_i = sum_j Acal_ij u_j  (acceleration)
        rhsSrc.assign(N,{0,0,0});
        for(int p=0;p<N;p++)for(int i=0;i<3;i++){
            double s=0; for(int j=0;j<3;j++) s += Acal[i][j]*u[p][j];
            rhsSrc[p][i]=s;
        }
    }
    void advanceExplicit(double dt){ // u += dt*rhsSrc  (no rhsMix: nu=0; no eddies)
        for(int p=0;p<N;p++)for(int i=0;i<3;i++) u[p][i] += dt*rhsSrc[p][i];
    }
};

// whiten an isotropic field so that R(0) = (2/3) I exactly (Cholesky rescale)
void initIsotropic(Domain&d, unsigned seed){
    mt19937 rng(seed); normal_distribution<double> g(0,1);
    for(auto&p:d.u) for(int i=0;i<3;i++) p[i]=g(rng);
    // remove mean
    array<double,3> m{0,0,0};
    for(auto&p:d.u)for(int i=0;i<3;i++)m[i]+=p[i];
    for(int i=0;i<3;i++)m[i]/=d.N;
    for(auto&p:d.u)for(int i=0;i<3;i++)p[i]-=m[i];
    // sample covariance -> Cholesky L ; transform x -> sqrt(2/3) L^{-1} x
    M3 C=d.reynoldsStress();
    double L00=sqrt(C[0][0]);
    double L10=C[1][0]/L00, L11=sqrt(C[1][1]-L10*L10);
    double L20=C[2][0]/L00, L21=(C[2][1]-L20*L10)/L11, L22=sqrt(C[2][2]-L20*L20-L21*L21);
    double s=sqrt(TWO3);
    for(auto&p:d.u){
        double x0=p[0], x1=p[1], x2=p[2];
        double y0=x0/L00;
        double y1=(x1-L10*y0)/L11;
        double y2=(x2-L20*y0-L21*y1)/L22;
        p[0]=s*y0; p[1]=s*y1; p[2]=s*y2;
    }
}

int main(){
    const double EMAX=4.0, dt=1e-4; const int N=4000;
    // Level 0 targets at e=4 (from the Python script), [u1^2,u2^2,u3^2]/2kt
    printf("closure   strain         u1^2/2kt  u2^2/2kt  u3^2/2kt   (Level 0 target)\n");
    printf("-------------------------------------------------------------------------\n");
    struct Case{ string name; M3 A; };
    double a=0.5, gax=1.0/sqrt(3.0);
    vector<Case> cases = {
        {"plane",        {{{a,0,0},{0,-a,0},{0,0,0}}}},
        {"axisymmetric", {{{-0.5*gax,0,0},{0,-0.5*gax,0},{0,0,gax}}}}
    };
    // hard-coded Level 0 targets for the printout
    auto target=[&](string c,string s)->string{
        if(c=="IP"  && s=="plane")        return "(0.1183 0.6500 0.2316)";
        if(c=="LRR" && s=="plane")        return "(0.0974 0.6006 0.3020)";
        if(c=="IP"  && s=="axisymmetric") return "(0.4499 0.4499 0.1001)";
        if(c=="LRR" && s=="axisymmetric") return "(0.4683 0.4683 0.0635)";
        return "";
    };
    for(string closure : {"IP","LRR"}){
        for(auto&cs:cases){
            Domain d; d.N=N; d.u.assign(N,{0,0,0}); d.A=cs.A; d.closure=closure;
            initIsotropic(d, 12345);
            int nsteps=(int)llround(EMAX/dt);
            for(int n=0;n<nsteps;n++){      // mirror micromixer explicit substep loop
                d.updateStrainOperator();   // recompute B(R), Acal  (per substep)
                d.getRhsSrc();              // dv_uvw source term
                d.advanceExplicit(dt);      // u += dt*rhsSrc
            }
            M3 R=d.reynoldsStress(); double kt=0.5*trace3(R);
            printf("%-8s  %-13s  %8.4f  %8.4f  %8.4f   %s\n",
                   closure.c_str(), cs.name.c_str(),
                   R[0][0]/(2*kt), R[1][1]/(2*kt), R[2][2]/(2*kt),
                   target(closure,cs.name).c_str());
        }
    }
    return 0;
}
