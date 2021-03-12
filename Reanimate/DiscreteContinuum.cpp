//
// Created by Paul Sweeney on 21/02/2021.
//

#include "DiscreteContinuum.hpp"
#include "Vasculature.hpp"

using namespace reanimate;

DiscreteContinuum::DiscreteContinuum() {}
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

    rnod = zeros<mat>(nnodT,nnodT);
    Mtiss = zeros<mat>(nnodT,nnodT);
    snode = zeros<mat>(3,nnodT);
    lambdaRange = zeros<mat>(2,nIterLambda);

/*    xsl0 = zeros<vec>(3);
    xsl1 = zeros<vec>(3);
    xsl2 = zeros<vec>(3);*/

}

void DiscreteContinuum::runHybrid()    {

    printText("Discrete-Continuum Module", 3);
    discreteNet.silence = true;

    // Map geometry classification to full network
    graph.mapClassification(discreteNet);

    // Calculate no. of arterioles/capillaries/venules
    nart = (int) accu(discreteNet.vesstyp(find(discreteNet.vesstyp == 1)));
    ncap = (int) accu(discreteNet.vesstyp(find(discreteNet.vesstyp == 2))) / 2;
    nven = (int) accu(discreteNet.vesstyp(find(discreteNet.vesstyp == 3))) / 3;

    // Discrete boundary flows
    artIn = accu(-discreteNet.BCflow(find(discreteNet.BCgeo == 1)));
    venOut = accu(-discreteNet.BCflow(find(discreteNet.BCgeo == 3)));
    capFlow = abs(artIn) - abs(venOut);
    qact = artIn + venOut;
    printNum("Arteriolar inflow (nl/min) =", artIn);
    printNum("Venular outflow (nl/min) =", venOut);
    printNum("Capillary flow exchange (nl/min) =", capFlow);

    analyseBranches();
    bridgeFlow();
    computeDiscrete();
    computeContinuum();

}