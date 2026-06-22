/**
 * @file domain.cc
 * @brief Source file for class \ref domain
 */

#include "domain.h"
#include "processor.h"
#include "dv.h"
#include "dv_uvw.h"
#include "dv_pos.h"
#include "dv_posf.h"
#include "dv_rho.h"
#include "dv_dvisc.h"
#include "dv_sca.h"
#include "dv_aDL.h"
#include "domaincase_odt_channel.h"
#include "domaincase_odt_channelScalar.h"
#include "domaincase_odt_isothermalWall.h"
#include "domaincase_odt_jetMixlRxn.h"
#include "domaincase_odt_jetFlame.h"
#include "domaincase_odt_MFjetFlame.h"
#include "domaincase_odt_coldPropaneJet.h"
#include "domaincase_odt_coldJet.h"
#include "domaincase_odt_RT.h"
#include "domaincase_odt_homogeneousStrain.h"
#include <cmath>
#include <iomanip>

extern processor proc;
using Cantera::Solution;

/////////////////////////////////////////////////////////////////////
/** Constructor
 */

domain::domain(domain *p_domn, param *p_pram) {

    domn = p_domn;
    pram = p_pram;
    domc = 0;               // initialize for destruction of eddy domains

 }

/////////////////////////////////////////////////////////////////////
/** Initializer
 */

void domain::init(inputoutput         *p_io,
                  meshManager         *p_mesher,
                  streams             *p_strm,
                  shared_ptr<Solution> csol,
                  micromixer          *p_mimx,
                  eddy                *p_ed,
                  domain              *p_eddl,
                  solver              *p_solv,
                  randomGenerator     *p_rand,
                  bool             LisEddyDomain) {

    //----------------------

    io     = p_io;
    mesher = p_mesher;
    strm   = p_strm;
    mimx   = p_mimx;
    ed     = p_ed;
    eddl   = p_eddl;
    solv   = p_solv;
    rand   = p_rand;

    gas = csol->thermo(); 
    kin = csol->kinetics(); 
    trn = csol->transport(); 

    //----------------------

    ngrd    = pram->ngrd0;
    ngrdf   = ngrd + 1;

    //----------------------

    if(LisEddyDomain) {        // eddy domain needs less data
        initEddyDomain();
        return;
    }

    //----------------------
    io->init(this);
    pram->init(this);
    ed->init(this, eddl);
    solv->init(this);
    // mesher is init below in caseinit for phi
    // strm is init below in caseinit  (domc), (if needed)
    // mimx is init below since it needs v[] set for cvode

    //---------------------- Continue setting up the case using the case_somecase class.
    // Adds to the above variable list, and initializes solution for the run

     if(pram->probType == "CHANNEL")
         domc = new domaincase_odt_channel();    // cold channel flow

     else if(pram->probType == "CHANNEL_SCALAR")
         domc = new domaincase_odt_channelScalar();  // cold channel flow with passive scalar

     else if(pram->probType == "JETMIXL_RXN")
         domc = new domaincase_odt_jetMixlRxn(); // jet, wake, mixing layer with gaseous reaction

     else if(pram->probType == "COLDPROPANEJET")
         domc = new domaincase_odt_coldPropaneJet(); // TNF jet

     else if(pram->probType == "COLDJET")
         domc = new domaincase_odt_coldJet(); // Hussein 1994

     else if(pram->probType == "JETFLAME")
         domc = new domaincase_odt_jetFlame(); // Shaddix jet

     else if(pram->probType == "MF_JETFLAME")
         domc = new domaincase_odt_MFjetFlame(); // jet flame w/ mixt frac density profile

     else if(pram->probType == "ISOTHERMAL_WALL")
         domc = new domaincase_odt_isothermalWall(); // isothermal wall

     else if(pram->probType == "RT")
         domc = new domaincase_odt_RT();      // simple Rayleigh Taylor flow

     else if(pram->probType == "HOMOGENEOUS_STRAIN")
         domc = new domaincase_odt_homogeneousStrain(); // strain-coupled ODT (Level 1a)

     else {
         cout << endl << "ERROR, probType UNKNOWN" << endl;
         exit(0);
     }

    domc->init(this);

    //----------------------

    for(int k=0; k<v.size(); k++)
        varMap[v.at(k)->var_name] = v.at(k);

    nTrans = 0;
    for(int k=0; k<v.size(); k++)
        if(v.at(k)->L_transported)
            nTrans++;

    //----------------------

    mimx->init(this);

    //----------------------

    if(pram->Lrestart) {
        io->loadVarsFromRestartFile();
        io->set_iNextDumpTime(pram->trst);
    }

}

/////////////////////////////////////////////////////////////////////
/** Compute size of domain based on faces.
 */

double domain::Ldomain() {
     return posf->d.at(ngrd) - posf->d.at(0);
}

/////////////////////////////////////////////////////////////////////
/** Initialize data members of the eddy domain.
 *  Note, none of the other members of this domain should be used (like random).
 *  Note, all variables here should have corresponding variables (by var_name) in the
 *     main domn. This is needed for using the eddl object.
 */

void domain::initEddyDomain() {

    v.push_back(new dv_pos(  this, "pos",   false, true));
    v.push_back(new dv_posf( this, "posf",  false, true));
    v.push_back(new dv_uvw(  this, "uvel",  true,  true));   // last are: L_transported, L_output
    v.push_back(new dv_uvw(  this, "vvel",  true,  true));
    v.push_back(new dv_uvw(  this, "wvel",  true,  true));
    v.push_back(new dv_rho(  this, "rho",   false, false));
    v.push_back(new dv_dvisc(this, "dvisc", false, false));
    if(domn->pram->LdoDL)
       v.push_back(new dv_aDL(this, "aDL",   false, false));

    int k = 0;
    pos   = v.at(k++);
    posf  = v.at(k++);
    uvel  = v.at(k++);
    vvel  = v.at(k++);
    wvel  = v.at(k++);
    rho   = v.at(k++);
    dvisc = v.at(k++);
    if(domn->pram->LdoDL)
        aDL   = v.at(k++);

}

/////////////////////////////////////////////////////////////////////
/** Set the domain from a region of the domn.  Normally called by eddy domain.
 *  @param i1 \input index of starting cell of domn to build from
 *  @param i2 \input index of ending cell of domn to build from
 *  If i2 < i1, we have a periodic region (wrap around the domain).
 *     This only happens in planar cases, not cylindrical or sphericial.
 *  nonwrap: |   | * | * | * | * | * | * |   |   |
 *                i1                  i2
 *  new domain consists of *'d cells
 *
 *  Wrap: | 4 | 5 |   |   |   |   | 1 | 2 | 3 |
 *             i2                  i1
 *  New domain consists of #'d cells: 1 2 3 4 5}
 */

void domain::setDomainFromRegion(const int i1, const int i2) {

    ngrd  = i2-i1+1;
    ngrdf = ngrd+1;

    for(int k=0; k<v.size(); k++)
        v.at(k)->setDvFromRegion(i1,i2);
}

///////////////////////////////////////////////////////////////////////////////
/** Find index of cell for given position (residing in cell).
 *  Start search assuming a uniform grid,
 *  then search forward or back till hit the cell index.
 *  If position is on cell face j, then if LowSide true, return j, else j-1.         \n
 * For start of eddy region, set LowSide to true                                     \n
 * For end of eddy region, set LowSide to false                                      \n
 * (This is so triplet maps don't overlap cells)
 *                                                                                   <pre><code>
 * e.g., usual:   | { | | | | } |    5 pts, eddy pos between cell faces
 *       okay:    {   | | | |   }    5 pts, eddy pos on cell faces (1 or both)
 *       bad:     |   { | | }   |    5 pts, eddy pos on internal faces (1 or both)
 *                                                                                   </code></pre>
 * @param position \input position to find corresponding index.
 * @param LowSide  \input flag true, then return j if position is on cell face j, else j-1.
 * @return index of position.
 */

int domain::domainPositionToIndex(double position, const bool LowSide, int dbg) {

    // strain-coupled ODT: with the line dilatation active, the boundary can move
    // inward by up to one timestep's compression between when an eddy edge is
    // sampled/clamped and when it is indexed, leaving the edge a hair outside the
    // current grid. Clamp such sub-step overshoots to the boundary rather than
    // aborting; genuine (large) out-of-domain positions still error below.
    if(pram->Lstrain && pram->Ldilatation) {
        const double tol = 1.0e-3 * Ldomain();
        if(position < posf->d.at(0)    && (posf->d.at(0)    - position) < tol) position = posf->d.at(0);
        if(position > posf->d.at(ngrd) && (position - posf->d.at(ngrd)) < tol) position = posf->d.at(ngrd);
    }

    if(abs(position-posf->d.at(0)) < 1.0E-14)
        return 0;
    if(abs(position-posf->d.at(ngrd)) < 1.0E-14)
        return ngrd-1;

    //if(position < posf->d.at(0))         // for periodic (from eddies only)
    //    position += Ldomain();
    //if(position > posf->d.at(ngrd))
    //    position -= Ldomain();

    if(position < posf->d.at(0) || position > posf->d.at(ngrd)) {
       *io->ostrm << "\ndbg = " << dbg << endl; //doldb
       *io->ostrm << scientific;
       *io->ostrm << setprecision(14);
       *io->ostrm << "\n ERROR odt_grid::domainPositionToIndex position < posf->d.at(0) or > posf->d.at(ngrd) \n"
               " and at processor's id---> " << proc.myid
               <<" Value of position is---> "<<position << " and values of posf->d.at(0) and posf->d.at(ngrd) are "
               <<posf->d.at(0)<< " and "<<posf->d.at(ngrd) <<" respectively "<< endl;
       //io->outputProperties("dbg.dat", 0.0); //doldb
       exit(0);
    }

    int i;
    int ipos = static_cast<int>((position-posf->d.at(0))/Ldomain()*ngrd);

    if(posf->d.at(ipos+1) > position) {      // case 1: grd skewed more pts on right half
        for(i=ipos+1; i>=0; i--)  {
            if(posf->d.at(i) <= position) {
                if(position == posf->d.at(i)) {
                    if(LowSide)
                        return i;
                    else
                        return i-1;
                }
                else
                    return i;
            }
        }
    }

    else  {                           // case 2: grd skewed more pts on left half
        for(i=ipos+1; i<=ngrdf; i++) {
            if(posf->d.at(i) >= position) {
                if(position == posf->d.at(i)) {
                    if(LowSide)
                        return i;
                    else
                        return i-1;
                }
                else
                    return i-1;
            }
        }
    }

    *io->ostrm << "\n\n******** ERROR IN odt_grid::domainPositionToIndex "
         << position << '\t' << posf->d.at(0) << '\t' << posf->d.at(ngrd) << '\t' << endl << endl;

    return -1;
}

/////////////////////////////////////////////////////////////////////
/** Cycle domain for periodic flows.
 *  @param icycle \input move all cells before and including this one
 *   to the end of the domain.
 *  @return the cycle distance (used for backcycling).
 */

double domain::cyclePeriodicDomain(const int icycle) {

    double cycleDistance = posf->d.at(icycle+1)-posf->d.at(0);

    for(int k=0; k<v.size(); k++) {
        if (v.at(k)->var_name=="pos" || v.at(k)->var_name=="posf")
            continue;
        v.at(k)->d.insert(v.at(k)->d.end(),   v.at(k)->d.begin(), v.at(k)->d.begin()+icycle+1);
        v.at(k)->d.erase( v.at(k)->d.begin(), v.at(k)->d.begin()+icycle+1);
    }

    //---------- now do posf, and pos

    double xend = posf->d.at(ngrd);
    for(int i=1; i<=icycle+1; i++)
        posf->d.push_back(xend+(posf->d.at(i)-posf->d.at(0)));
    posf->d.erase(posf->d.begin(), posf->d.begin()+icycle+1);

    pos->setVar();     // does a little extra work (whole domain) but doesn't happen that often
                       //    only when periodic eddies are accepted.

    return cycleDistance;
}

/////////////////////////////////////////////////////////////////////
/** Back cycle domain for periodic flows. Intended to be called some time
 *  after cyclePeriodicDomain is called.
 *  Splits the cell at posf.at(ngrd) - backCycleDistace, then moves end cells
 *     after the split to the beginning of the domain.
 *  @param \input distance from the end to split and move the domain.
 */

void domain::backCyclePeriodicDomain(const double backCycleDistance) {

    double xend = posf->d.at(ngrd) - backCycleDistance;     // end loc.
    double icycle = domainPositionToIndex(xend, true, 1);  // cycle cells greater than this to beginning

    //------------ split the cell where the back cycle happens

    vector<double> interPos(3);
    if(abs(posf->d.at(icycle) - xend) > 1.0e-15) {
        interPos.at(0) = posf->d.at(icycle);
        interPos.at(1) = xend;
        interPos.at(0) = posf->d.at(icycle+1);
        mesher->splitCell(icycle, 1, interPos);
        icycle++;
    }

    //------------ now move the cells

    int nmove = ngrd-icycle+1;

    for(int k=0; k<v.size(); k++) {
        if (v.at(k)->var_name=="pos" || v.at(k)->var_name=="posf")
            continue;
        v.at(k)->d.insert(v.at(k)->d.begin(), v.at(k)->d.begin()+icycle, v.at(k)->d.end());
        v.at(k)->d.erase( v.at(k)->d.begin()+icycle+nmove, v.at(k)->d.end() );
    }

    //---------- now do posf, and pos

    double xstart_orig = posf->d.at(0);

    posf->d.insert(posf->d.begin(), nmove, 0.0);
    icycle += nmove;
    for(int i=0; i<nmove; i++)
        posf->d.at(i) = xstart_orig - (posf->d.at(posf->d.size()-1-i) - xend);

    posf->d.erase(posf->d.begin()+icycle+1, posf->d.end());

    pos->setVar();     // does a little extra work (whole domain) but doesn't happen that often
                       //    only when periodic eddies are accepted.
}

///////////////////////////////////////////////////////////////////////////////
/** Strain-coupled ODT: compute the line Reynolds stress, the rapid
 *  pressure-strain operator B_ij from the Lyapunov solve B R + R B = Pi^r,
 *  and store the combined operator Acal = -A + B in pram->Acal. Called once
 *  per explicit substep from micromixer when pram->Lstrain is true.
 */

namespace {

    // cyclic Jacobi eigensolver for a symmetric 3x3 matrix: A = Q diag(lam) Q^T
    void jacobi3(const double Ain[3][3], double Q[3][3], double lam[3]) {
        double a[3][3];
        for(int i=0;i<3;i++) for(int j=0;j<3;j++){ a[i][j]=Ain[i][j]; Q[i][j]=(i==j)?1.0:0.0; }
        for(int sweep=0; sweep<50; sweep++){
            double off = std::fabs(a[0][1])+std::fabs(a[0][2])+std::fabs(a[1][2]);
            if(off < 1e-30) break;
            for(int p=0;p<2;p++) for(int q=p+1;q<3;q++){
                if(std::fabs(a[p][q]) < 1e-300) continue;
                double th = (a[q][q]-a[p][p])/(2.0*a[p][q]);
                double tt = (th>=0?1.0:-1.0)/(std::fabs(th)+std::sqrt(th*th+1.0));
                double c  = 1.0/std::sqrt(tt*tt+1.0), s=tt*c;
                for(int k=0;k<3;k++){ double akp=a[k][p],akq=a[k][q]; a[k][p]=c*akp-s*akq; a[k][q]=s*akp+c*akq; }
                for(int k=0;k<3;k++){ double apk=a[p][k],aqk=a[q][k]; a[p][k]=c*apk-s*aqk; a[q][k]=s*apk+c*aqk; }
                for(int k=0;k<3;k++){ double qkp=Q[k][p],qkq=Q[k][q]; Q[k][p]=c*qkp-s*qkq; Q[k][q]=s*qkp+c*qkq; }
            }
        }
        for(int i=0;i<3;i++) lam[i]=a[i][i];
    }

    // solve  B R + R B = P  for symmetric B (R symmetric positive-definite)
    void lyapunovSym(const double R[3][3], const double P[3][3], double B[3][3]) {
        double Q[3][3], lam[3]; jacobi3(R, Q, lam);
        double QtP[3][3], Pp[3][3];
        for(int i=0;i<3;i++) for(int j=0;j<3;j++){ double s=0; for(int k=0;k<3;k++) s+=Q[k][i]*P[k][j]; QtP[i][j]=s; }
        for(int i=0;i<3;i++) for(int j=0;j<3;j++){ double s=0; for(int k=0;k<3;k++) s+=QtP[i][k]*Q[k][j]; Pp[i][j]=s; }
        double Bp[3][3];
        for(int i=0;i<3;i++) for(int j=0;j<3;j++) Bp[i][j]=Pp[i][j]/(lam[i]+lam[j]);
        double QB[3][3];
        for(int i=0;i<3;i++) for(int j=0;j<3;j++){ double s=0; for(int k=0;k<3;k++) s+=Q[i][k]*Bp[k][j]; QB[i][j]=s; }
        for(int i=0;i<3;i++) for(int j=0;j<3;j++){ double s=0; for(int k=0;k<3;k++) s+=QB[i][k]*Q[j][k]; B[i][j]=s; }
    }
}

void domain::updateStrainOperator() {

    const vector<double> &u = uvel->d, &v = vvel->d, &w = wvel->d;

    // --- length-weighted Reynolds stress, fluctuation about the line mean ---
    //     (this block is the single averaging-mode seam: swap for an ensemble
    //      or windowed average at Level 2 without touching anything else)
    double Lsum=0, mu=0, mv=0, mw=0;
    for(int i=0; i<ngrd; i++){
        double dx = std::fabs(posf->d.at(i+1)-posf->d.at(i));
        Lsum += dx; mu += u[i]*dx; mv += v[i]*dx; mw += w[i]*dx;
    }
    mu/=Lsum; mv/=Lsum; mw/=Lsum;
    double R[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    for(int i=0; i<ngrd; i++){
        double dx = std::fabs(posf->d.at(i+1)-posf->d.at(i));
        double f[3] = {u[i]-mu, v[i]-mv, w[i]-mw};
        for(int a=0;a<3;a++) for(int b=0;b<3;b++) R[a][b] += f[a]*f[b]*dx;
    }
    for(int a=0;a<3;a++) for(int b=0;b<3;b++) R[a][b] /= Lsum;
    double kt = 0.5*(R[0][0]+R[1][1]+R[2][2]);

    // --- mean strain S and rotation W from pram->Astrain ---
    double A[3][3], S[3][3], W[3][3];
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) A[i][j]=pram->Astrain[i][j];
    for(int i=0;i<3;i++) for(int j=0;j<3;j++){
        S[i][j]=0.5*(A[i][j]+A[j][i]); W[i][j]=0.5*(A[i][j]-A[j][i]);
    }

    // --- production P_ij = -(A R + R A^T) ---
    double P[3][3];
    for(int i=0;i<3;i++) for(int j=0;j<3;j++){
        double s=0; for(int k=0;k<3;k++) s += A[i][k]*R[k][j] + R[i][k]*A[j][k];
        P[i][j] = -s;
    }

    // --- rapid pressure-strain Pi^r per closure ---
    double Pir[3][3];
    if(pram->strainClosure == "IP") {                 // -C2 (P - 1/3 trP I), C2=3/5
        double trP = P[0][0]+P[1][1]+P[2][2];
        for(int i=0;i<3;i++) for(int j=0;j<3;j++)
            Pir[i][j] = -0.6*(P[i][j] - (i==j?trP/3.0:0.0));
    }
    else {                                            // LRR-QI: C2=4/5,C3=7/4,C4=131/100
        double b[3][3];
        for(int i=0;i<3;i++) for(int j=0;j<3;j++) b[i][j]=R[i][j]/(2*kt)-(i==j?1.0/3.0:0.0);
        double trbS=0; for(int i=0;i<3;i++) for(int k=0;k<3;k++) trbS += b[i][k]*S[i][k];
        for(int i=0;i<3;i++) for(int j=0;j<3;j++){
            double bS_Sb=0, Wb_bW=0;
            for(int k=0;k<3;k++){ bS_Sb += b[i][k]*S[j][k]+S[i][k]*b[j][k];
                                  Wb_bW += W[i][k]*b[j][k]-b[i][k]*W[j][k]; }
            Pir[i][j] = 0.8*kt*S[i][j] + 1.75*kt*(bS_Sb-(i==j?(2.0/3.0)*trbS:0.0))
                      + 1.31*kt*Wb_bW;
        }
    }

    // --- Lyapunov solve for B, then Acal = -A + B ---
    double B[3][3]; lyapunovSym(R, Pir, B);
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) pram->Acal[i][j] = -A[i][j] + B[i][j];
}

////////////////////////////////////////////////////////////////////////////////
/** Strain-coupled ODT: apply the mean-strain dilatation of the line.
 *  The ODT line is aligned with the x_2 direction (the compression axis of the
 *  stagnation-blocking strain), so its length evolves as d(ln L)/dt = A_22.
 *  Over a substep this is the affine scaling  x -> xc + (x - xc) exp(A_22 dt)
 *  about the domain centre xc. This geometric compression is what shifts the
 *  resolved wavenumbers (k_2 grows as the line shortens) and is distinct from
 *  the velocity-amplitude strain applied in dv_uvw::getRhsSrc. The mesh length
 *  bounds and eddy-size scales are scaled by the same factor so the relative
 *  resolution and eddy population track the compressing line.
 *  Gated by pram->Ldilatation; a no-op otherwise. With bcType=WALL and zero wall
 *  velocities, meshManager::enforceDomainSize() is a no-op and does not undo it.
 */
void domain::applyStrainDilatation(const double dt) {

    if(!pram->Lstrain || !pram->Ldilatation) return;

    const double A22 = pram->Astrain[1][1];           // d(ln L)/dt
    const double f   = std::exp(A22 * dt);             // affine scale factor
    const double xc  = pram->xDomainCenter;

    for(int i=0; i<ngrdf; i++) posf->d.at(i) = xc + (posf->d.at(i) - xc)*f;
    for(int i=0; i<ngrd;  i++) pos ->d.at(i) = xc + (pos ->d.at(i) - xc)*f;

    // keep mesh-control lengths proportional to the (shrinking) line
    pram->dxmin *= f;
    pram->dxmax *= f;
    pram->Lmin  *= f;
    pram->Lmax  *= f;
    pram->Lp    *= f;
}
