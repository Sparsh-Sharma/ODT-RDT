/**
 * @file domaincase_odt_homogeneousStrain.h
 * @brief Header for class domaincase_odt_homogeneousStrain
 *
 * Strain-coupled ODT verification case (Level 1a): statistically homogeneous
 * turbulence on a periodic line, subjected to a constant imposed mean velocity
 * gradient A_ij (param Astrain) via the strain operator in dv_uvw::getRhsSrc.
 * With LnoEddies=true and kvisc0=0 this reproduces the rapid-distortion limit.
 */

#pragma once

#include "domaincase.h"
#include <vector>
#include <string>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

class domaincase_odt_homogeneousStrain : public domaincase {

    public:

        virtual void init(domain *p_domn);
        virtual void setCaseSpecificVars();

        domaincase_odt_homogeneousStrain(){}
        virtual ~domaincase_odt_homogeneousStrain(){}

};
