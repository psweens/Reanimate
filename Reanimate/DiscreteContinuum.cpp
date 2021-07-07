//
// Created by Paul Sweeney on 21/02/2021.
//

#include "DiscreteContinuum.hpp"

using namespace reanimate;

DiscreteContinuum::DiscreteContinuum() = default;
DiscreteContinuum::~DiscreteContinuum() = default;

void DiscreteContinuum::setup_hybridArrays() {sourceTyp = zeros<ivec>(discreteNet.getNnod());}

void DiscreteContinuum::setup_continuumArrays() {

    qout = zeros<vec>(nnodT);
    pout = zeros<vec>(nnodT);
    continuumPress = zeros<vec>(nnodT);
    r0 = zeros<vec>(nnodT);
    Mpress_pred = zeros<vec>(nIterLambda);
    sigmaFlow_pred = zeros<vec>(nIterLambda);
    sigmaPress_pred = zeros<vec>(nIterLambda);
    R2_flow = zeros<vec>(nIterLambda);
    R2_press = zeros<vec>(nIterLambda);

    NINV = zeros<mat>(nnodT,nnodT);
    rnod = zeros<mat>(nnodT,nnodT);
    Mtiss = zeros<mat>(nnodT,nnodT);

}

void DiscreteContinuum::runHybrid()    {

    printText("Discrete-Continuum Module", 3);
    discreteNet.silence = true;

    // Map geometry classification to full network
    graph.mapClassification(discreteNet);
    discreteNet.pictureNetwork("Hybrid_Classification.ps",conv_to<vec>::from(discreteNet.vesstyp));
    Vasculature clone = discreteNet;

    // Calculate no. of arterioles/capillaries/venules
    nart = (int) accu(discreteNet.vesstyp(find(discreteNet.vesstyp == 1)));
    ncap = (int) accu(discreteNet.vesstyp(find(discreteNet.vesstyp == 2))) / 2;
    nven = (int) accu(discreteNet.vesstyp(find(discreteNet.vesstyp == 3))) / 3;

    // Discrete boundary flows
    artIn = accu(discreteNet.BCflow(find(discreteNet.BCgeo == 1)));
    venOut = accu(discreteNet.BCflow(find(discreteNet.BCgeo == 3)));
    capFlow = accu(discreteNet.BCflow(find(discreteNet.BCgeo == 2)));
    qact = -capFlow; // Neg. signs comes from considering inflow/outflow of discrete network compared to sources/sinks in continuum domain

    printNum("Arteriolar inflow (nl/min) =", artIn);
    printNum("Venular outflow (nl/min) =", venOut);
    printNum("Capillary Influx (nl/min) =", capFlow);

    analyseBranches();
    bridgeFlow();
    computeDiscrete();
    computeContinuum();

    plotContour("Hybrid_PressureContourPlot.ps", discreteNet, Pa_Pv.max(), Pa_Pv.min(), false, true, nnodT, nnodT, 60);
    plotContour("Hybrid_SpeedContourPlot.ps", discreteNet, 0., 0., true, false, nnodT, nnodT, 60);

    mapContinuum(clone);

}
