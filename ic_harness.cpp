// Standalone compile+logic check of the spectral IC from
// domaincase_odt_homogeneousStrain::init (lines 54-121), verbatim math.
#include <vector>
#include <random>
#include <cmath>
#include <cstdio>
using namespace std;
int main(){
    // --- stand-ins for domn->pram and domn->ngrd ---
    int    N           = 4000;
    int    seed        = 22;
    double domainLength= 1.0;
    double xDomainCenter=0.0;
    double specKpWaves = 8.0;
    int    specNmodes  = 64;

    vector<double> u(N), v(N), w(N);
    std::mt19937 rng(seed >= 0 ? seed : 22);

    // ---- verbatim from the domaincase ----
    const double L  = domainLength;
    const double x0 = xDomainCenter - 0.5*L;
    const double dx = L / N;
    const double dk = 2.0*M_PI / L;
    const double kp = 2.0*M_PI*specKpWaves / L;
    const int    Nm = specNmodes;

    vector<double> amp(Nm+1, 0.0);
    for(int n=1;n<=Nm;n++){
        double r = (n*dk)/kp;
        amp[n] = std::sqrt( std::pow(r,4.0)*std::exp(-2.0*r*r) );
    }
    std::uniform_real_distribution<double> uni(0.0, 2.0*M_PI);
    vector<double>* comp[3] = {&u, &v, &w};
    for(int c=0;c<3;c++){
        vector<double> ph(Nm+1);
        for(int n=1;n<=Nm;n++) ph[n] = uni(rng);
        vector<double> &f = *comp[c];
        for(int i=0;i<N;i++){
            double y = x0 + (i+0.5)*dx;
            double s = 0.0;
            for(int n=1;n<=Nm;n++) s += amp[n]*std::cos(n*dk*y + ph[n]);
            f[i] = s;
        }
    }
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

    double L00=std::sqrt(C[0][0]);
    double L10=C[1][0]/L00, L11=std::sqrt(C[1][1]-L10*L10);
    double L20=C[2][0]/L00, L21=(C[2][1]-L20*L10)/L11,
           L22=std::sqrt(C[2][2]-L20*L20-L21*L21);
    double s=std::sqrt(2.0/3.0);
    for(int i=0;i<N;i++){
        double x0v=u[i], x1=v[i], x2=w[i];
        double y0=x0v/L00;
        double y1=(x1-L10*y0)/L11;
        double y2=(x2-L20*y0-L21*y1)/L22;
        u[i]=s*y0; v[i]=s*y1; w[i]=s*y2;
    }
    // ---- check R = (2/3) I ----
    double R[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    for(int i=0;i<N;i++){
        double f[3]={u[i],v[i],w[i]};
        for(int a=0;a<3;a++) for(int b=0;b<3;b++) R[a][b]+=f[a]*f[b];
    }
    for(int a=0;a<3;a++) for(int b=0;b<3;b++) R[a][b]/=N;
    double kt=0.5*(R[0][0]+R[1][1]+R[2][2]);
    printf("R = [%8.5f %8.5f %8.5f]\n", R[0][0],R[0][1],R[0][2]);
    printf("    [%8.5f %8.5f %8.5f]\n", R[1][0],R[1][1],R[1][2]);
    printf("    [%8.5f %8.5f %8.5f]\n", R[2][0],R[2][1],R[2][2]);
    printf("k_t = %.6f  (target 1.0)\n", kt);
    printf("fractions = %.5f %.5f %.5f  (target 0.33333)\n",
           R[0][0]/(2*kt),R[1][1]/(2*kt),R[2][2]/(2*kt));
    return 0;
}
