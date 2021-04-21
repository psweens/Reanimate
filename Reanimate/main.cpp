#include "Network.hpp"
#include "spatGraph.hpp"
#include "Vasculature.hpp"
#include "omp.h"

using namespace reanimate;
using namespace std;

int main(int argc, char** argv) {

    Vasculature vnet;
    vnet.setStackSize(); // Prevent stack overflow when running recursive functions (for large networks)
    //vnet.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data_LS3_" + to_string(file) + "/";
    vnet.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/";
    vnet.loadPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/";
    vnet.setBuildPath(); // Create build directory / delete contents

    vnet.loadNetwork("1pNetwork.dat");
    vnet.bloodFlow(true, true, false, false);
    vnet.printVisuals();
    vnet.printNetwork("solved_NetworkBloodFlow.txt");


    /*mat data;

    #pragma omp parallel for default(none) schedule(dynamic)
    for (int file = 1; file <= 12; file++)    {

        string filename = to_string(file) + "Network.txt";

        // Directory setup

        Vasculature vnet;
        vnet.setStackSize(); // Prevent stack overflow when running recursive functions (for large networks)
        //vnet.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data_LS3_" + to_string(file) + "/";
        vnet.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data_" +to_string(file) + "/";
        vnet.loadPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/LS3/";
        vnet.setBuildPath(); // Create build directory / delete contents

        // Load network file
        vnet.loadNetwork(filename);

        //if (file == 1)  {data = zeros<mat>(vnet.getNseg(),12);}
        //data.col(file-1) = vnet.hd;

        mat extraD = zeros<mat>(vnet.getNseg(), 1);
        vnet.printAmira("amiraNetwork.am", extraD);
        vnet.hd.fill(0.4);
        vnet.bchd.fill(0.4);


        // Solving for blood flow
        //vnet.loadDeadEnds = true;
        vnet.bloodFlow(true, true, false, false);
        vnet.printVisuals();
        vnet.printNetwork("solved_NetworkBloodFlow.txt");

    }*/



}
