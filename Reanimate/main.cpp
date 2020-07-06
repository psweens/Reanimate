//
//  main.cpp
//  Vascular Flows 2.0
//
//  Flow and Pressure Solver in Discrete Vascular Networks making use of the Armadillo
//  C++ library for stability and speed of simulations (requires BLAS, LAPACK and SuperLU).
//
//  Created by Paul Sweeney on 06/05/2015.
//  Copyright (c) 2015 Paul Sweeney - University College London. All rights reserved.
//

#include "Network.hpp"
#include "spatGraph.hpp"
#include "Vasculature.hpp"

using namespace reanimate;
using namespace std;

int main(int argc, char** argv) {

    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) != nullptr)   {
        printf("Current working dir: %s\n", cwd);
    }
    else {
        perror("getcwd() error");
        return 1;
    }

    //gpu_test();

    string root, loadroot;

    root = cwd;
    root = "/home/sweene01/Documents/Reanimate/Build_Data/";
    loadroot = cwd;
    loadroot = "/home/sweene01/Documents/Reanimate/Load_Data/";

    Vasculature cortex;
    cortex.loadNetwork(loadroot + "Network.dat");
    cortex.pictureNetwork(root + "NetworkDiameters.ps",cortex.diam);
    cortex.bloodFlow(true, true, true);
    cortex.pictureNetwork(root + "NetworkFlow.ps",cortex.qq, true);
    cortex.pictureNetwork(root + "NetworkPressure.ps",cortex.segpress);
    cortex.pictureNetwork(root + "NetworkShearStress.ps",cortex.tau, true);
    cortex.pictureNetwork(root + "NetworkHD.ps",cortex.hd);
    cortex.edgeNetwork();

    spatGraph eCortex;
    eCortex.generate(cortex);
    eCortex.pictureNetwork(root + "eNetworkDiameters.ps",eCortex.diam);
    eCortex.loadTrunks(loadroot + "NetworkClass.txt");
    eCortex.pictureNetwork(root + "eNetworkGeometry.ps",conv_to<vec>::from(eCortex.geometry));

}
