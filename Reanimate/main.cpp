#include "Network.hpp"
#include "spatGraph.hpp"
#include "Vasculature.hpp"
#include "omp.h"

using namespace reanimate;
using namespace std;

int main(int argc, char** argv) {

    /*Vasculature vnet;
    vnet.setStackSize(); // Prevent stack overflow when running recursive functions (for large networks)
    //vnet.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data_LS3_" + to_string(file) + "/";
    vnet.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/";
    vnet.loadPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/";
    vnet.setBuildPath(); // Create build directory / delete contents

    vnet.loadNetwork("deadSWNetwork.txt");
    vnet.setTargetPressure(mean(vnet.bcprfl(find(vnet.bctyp == 0))));
    vnet.setTargetTau(15.);
    vnet.loadDeadEnds = true;
    vnet.bloodFlow(true, true, true, false);
    vnet.printVisuals();*/


    mat data;

    //omp_set_num_threads(12);
    //#pragma omp parallel for default(none) schedule(dynamic)
    for (int file = 1; file <= 12; file++)    {

        string filename = "SW1_"+to_string(file) + ".txt";

        // Directory setup

        Vasculature vnet;
        vnet.setStackSize(16*1024*1024); // Prevent stack overflow when running recursive functions (for large networks)
        //vnet.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data_LS3_" + to_string(file) + "/";
        vnet.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data_" +to_string(file) + "/";
        vnet.loadPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/SW1/";
        vnet.setBuildPath(); // Create build directory / delete contents

        // Load network file
        vnet.loadNetwork(filename);

        //if (file == 1)  {data = zeros<mat>(vnet.getNseg(),12);}
        //data.col(file-1) = vnet.hd;

        //mat extraD = zeros<mat>(vnet.getNseg(), 1);
        //vnet.printAmira("amiraNetwork.am", extraD);
        //vnet.hd.fill(0.4);
        //vnet.bchd.fill(0.4);


        // Solving for blood flow
        //vnet.loadDeadEnds = true;
        vnet.bloodFlow(true, false, false, false);
        //vnet.printVisuals();
        vnet.printNetwork("solved_NetworkBloodFlow.txt");

        vnet.printAscii("SW1_"+to_string(file)+"_ascii.txt");

    }



}
