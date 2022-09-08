//
// Created by Paul Sweeney on 21/02/2021.
//

#include "DiscreteContinuum.hpp"

#include <string>
#include <dirent.h>

using namespace reanimate;

DiscreteContinuum::DiscreteContinuum() {}
DiscreteContinuum::~DiscreteContinuum() = default;

void DiscreteContinuum::setBuildPath(const string path, bool deleteFiles) {

    buildPath = path;
    discreteNet.buildPath = path;
    cell.buildPath = path;

    printf("Setting build directory ...\n");
    DIR* dir = opendir(buildPath.c_str());
    if (dir) {
        if (deleteFiles) {
            printf("Directory exists, deleting contents\n");
            string contents = buildPath + "*";
            system(("rm " + contents).c_str());
        }
    }
    else {
        printf("Directory does not exist. Creating folder ...\n");
        system(("mkdir " + buildPath).c_str());
    }
    closedir(dir);

}

void DiscreteContinuum::setLoadPath(const string path) {

    loadPath = path;
    discreteNet.loadPath = path;
    cell.loadPath = path;

    printf("Setting build directory ...\n");
    DIR* dir = opendir(loadPath.c_str());
    if (!dir) {
        printf("Directory does not exist. Creating folder ...\n");
        system(("mkdir " + loadPath).c_str());
    }

    closedir(dir);

}

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

    // Map vesselGeometry classification to full network
    graph.mapClassification(discreteNet);
    discreteNet.pictureNetwork("Hybrid_Classification.ps",conv_to<vec>::from(discreteNet.vesstyp));
    discreteNet.pictureNetwork("Hybrid_InitialBloodPressure.ps", discreteNet.segpress);
    discreteNet.pictureNetwork("Hybrid_InitialBloodFlow.ps", abs(discreteNet.q), true);
    Vasculature clone = discreteNet;

    // Calculate no. of arterioles/capillaries/venules
    nart = (int) accu(discreteNet.vesstyp(find(discreteNet.vesstyp == 1)));
    ncap = (int) accu(discreteNet.vesstyp(find(discreteNet.vesstyp == 2))) / 2;
    nven = (int) accu(discreteNet.vesstyp(find(discreteNet.vesstyp == 3))) / 3;

    // Discrete boundary flows
    artIn = accu(discreteNet.BCflow(find(discreteNet.BCgeo == 1)));
    capFlow = accu(discreteNet.BCflow(find(discreteNet.BCgeo == 2)));
    venOut = accu(discreteNet.BCflow(find(discreteNet.BCgeo == 3)));
    qact = -capFlow; // Neg. signs comes from considering inflow/outflow of discrete network compared to sources/sinks in continuum domain

    printNum("Arteriolar Inflow (nl/min) =", artIn);
    printNum("Venular Outflow (nl/min) =", venOut);
    printNum("Capillary Flux (nl/min) =", capFlow);

    analyseBranches();
    bridgeFlow();
    computeDiscrete();
    computeContinuum();

    plotContour("Hybrid_PressureContourPlot.ps", discreteNet, 55, 20, false, false, nnodT, nnodT, 60);
    plotContour("Hybrid_SpeedContourPlot.ps", discreteNet, 5, -10, true, true, nnodT, nnodT, 60);

    mapContinuum(clone);

}
