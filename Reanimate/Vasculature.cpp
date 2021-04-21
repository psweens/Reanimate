//
// Created by sweene01 on 06/07/2020.
//

#include "Vasculature.hpp"

using namespace reanimate;

Vasculature::Vasculature() {

    varviscosity = true;
    phaseseparation = false;
    memoryeffects = false;
    updateBoundaryHD = false;
    loadDeadEnds = false;

    // Hd and q tolerances for variable Hd
    hdtol = 1.e-3;
    qtol = 1.e-3;

    bifpar = zeros<vec>(3);
    cpar = zeros<vec>(4);
    viscpar = zeros<vec>(6);

}
Vasculature::~Vasculature() = default;


void Vasculature::printSummary() {



}

void Vasculature::printVisuals()    {

    const char *headers[5] = {"Pressure","Flow","Hd","Velocity","WSS"};

    mat data = zeros<mat>(nseg, 5);
    data.col(0) = segpress;
    data.col(1) = qq;
    data.col(2) = hd;
    data.col(3) = qq;
    data.col(4) = tau;
    printAmira("amiraBloodFlow.am", data, true, headers);

    pictureNetwork("NetworkDiameters.ps", diam);
    pictureNetwork("NetworkPressure.ps", segpress);
    pictureNetwork("NetworkFlow.ps", qq);
    pictureNetwork("NetworkWSS.ps", tau);
    pictureNetwork("NetworkHaematocrit.ps", hd);

}