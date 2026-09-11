/**
 * @file eddy.h
 * @brief Header file for class \ref eddy
 */

#pragma once

#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing eddy object
 *
 *  @author David O. Lignell
 */

class eddy {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        domain              *domn;              ///< pointer to domain object
        domain              *eddl;              ///< pointer to eddy line object

        double              eddySize;           ///< size of eddy
        double              leftEdge;           ///< left edge location of eddy
        double              rightEdge;          ///< right edge location of eddy
        double              invTauEddy;         ///< inverse eddy timescale
        double              Pa;                 ///< eddy acceptance probability
        bool                LperiodicEddy;      ///< a wrap-around eddy
        double              curMidFrac;         ///< this candidate's middle-image volume fraction (Option B / mixture)
        double              lastAcceptedSize;   ///< size of the most recent ACCEPTED eddy (0 before the first); rescaled by applyStrainDilatation; feeds LrelaxLastEddySize clock events
        vector<double>      cCoef;              ///< coefficient of K kernel
        vector<double>      bCoef;              ///< coefficient of J kernel
        vector<double>      K;                  ///< eddy kernel K
        vector<double>      dxc;                ///< \delta(x^cCoord) is prop. to cell "volume"
        vector<double>      pos0;               ///< initial eddy cell locations, for kernel

        double              esdp1;              ///< eddy size distribution parameters.
        double              esdp2;
        double              esdp3;
        double              esdp4;

        vector<double> cca,ccb,ccc,ccd;                 ///< polynomial coefficient arrays for cylindricalAnomalyHack

    //////////////////// MEMBER FUNCTIONS /////////////////

        void   sampleEddySize();
        void   sampleEddyPosition();
        void   tripMap(domain *line,    const int iS, int iE, const double C, const bool LsplitAtEddy=false);
        bool   eddyTau(const double Z_value, const double C);
        void   computeEddyAcceptanceProb(const double dtSample);
        void   applyVelocityKernels(domain *line, const int iS, const int iE);
        void   applySubscaleKernels(domain *line, const int iS, const int iE, const int nlevels);
        void   applyKernelOnRange(domain *line, const int i0, const int i1);
        void   applyConcurrentRelax(domain *line, const int N);
        void   applyRelaxEvent(domain *line);
        bool   anisoReductionOK(const double fac);

    private:

        void   fillKernel();
        void   fillKernel_planarAnalytic();
        double eddyFavreAvgVelocity(const vector<double> &dxc);
        void   set_kernel_coefficients();

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        eddy(){}
        void init(domain *p_domn, domain *p_eddl);
        ~eddy(){}

};


////////////////////////////////////////////////////////////////////////////////


