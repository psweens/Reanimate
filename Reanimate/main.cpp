#include "Network.hpp"
#include "spatGraph.hpp"
#include "Vasculature.hpp"
#include <sys/resource.h>

using namespace reanimate;
using namespace std;

int main(int argc, char** argv) {

    // Directory setup
    Vasculature vnet;
    vnet.setStackSize(); // Prevent stack overflow when running recursive functions (for large networks)
    vnet.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/";
    vnet.loadPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/";
    vnet.setBuildPath(); // Create build directory / delete contents

    // Load network file
    vnet.loadNetwork("cpNetwork.txt");


    mat extraD = zeros<mat>(vnet.getNseg(), 1);
    vnet.printAmira("amiraNetwork.am", extraD);
    vnet.hd.fill(0.4);
    vnet.bchd.fill(0.4);


    // Solving for blood flow
    //vnet.loadDeadEnds = true;
    vnet.bloodFlow(true, true, false, false);

    extraD = zeros<mat>(vnet.getNseg(), 1);
    extraD.col(0) = vnet.segpress;
    vnet.printAmira("amiraPressure.am", extraD);
    extraD.col(0) = vnet.hd;
    vnet.printAmira("amiraHD.am", extraD);

    vnet.pictureNetwork("NetworkDiameters.ps", vnet.diam);
    vnet.pictureNetwork("NetworkPressure.ps", vnet.segpress);
    vnet.pictureNetwork("NetworkFlow.ps", vnet.qq);
    vnet.pictureNetwork("NetworkHaematocrit.ps", vnet.hd);
    vnet.printNetwork("solved_NetworkBloodFlow.txt");

}
