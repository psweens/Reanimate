//
// Created by sweene01 on 06/07/2020.
//

#include "Vasculature.hpp"

using namespace reanimate;

Vasculature::Vasculature() {

    varviscosity=true;
    phaseseparation=false;
    memoryeffects=false;

    // Hd and q tolerances for variable Hd
    hdtol = 1.e-3;
    qtol = 1.e-3;

    bifpar = zeros<vec>(3);
    cpar = zeros<vec>(4);
    viscpar = zeros<vec>(6);

}
Vasculature::~Vasculature() = default;