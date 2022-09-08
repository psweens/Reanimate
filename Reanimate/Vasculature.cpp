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

void Vasculature::printVisuals(bool amira, bool twoDim, bool logDistrib)    {

    qq(find(qq < 1.e-4)).fill(1.e-4);
    tau(find(tau < 1.e-4)).fill(1.e-4);

    if (amira)  {
        const char *headers[5] = {"Pressure","LogFlow","Hd","LogVelocity","LogWSS"};
        mat data = zeros<mat>(nseg, 5);
        data.col(0) = segpress;
        data.col(1) = log(qq);
        data.col(2) = hd;
        data.col(3) = log(vel);
        data.col(4) = log(tau);
        printAmira("amiraBloodFlow.am", data, true, headers);
    }

    if (twoDim) {
        pictureNetwork("Network_Diameters.ps", diam);
        pictureNetwork("Network_BloodPressure.ps", segpress);
        if (logDistrib) {
            pictureNetwork("Network_BloodFlow_log.ps", log(qq));
            pictureNetwork("Network_WSS_log.ps", log(tau));
        }
        else {
            pictureNetwork("Network_BloodFlow.ps", qq);
            pictureNetwork("Network_WSS.ps", tau);
        }
        pictureNetwork("Network_Haematocrit.ps", hd);
    }

}