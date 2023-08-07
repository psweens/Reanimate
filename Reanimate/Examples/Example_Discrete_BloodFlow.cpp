#include "Network.hpp"
#include "spatGraph.hpp"
#include "Vasculature.hpp"
#include "DiscreteContinuum.hpp"
#include "MicroCell.hpp"
#include "omp.h"
#include <iostream>
#include <fstream>

using namespace reanimate;
using namespace std;

void example_Discrete_BloodFlow(std::string file_path, std::string file_name, std::string build_path) {

    // Set number of threads for any multithreading
    //omp_set_num_threads(48);

    // Define vasculature
    Vasculature network;
    // Define and set load / build paths
    network.setLoadPath(file_path);
    // "/mnt/ml/anaconda_envs/vessel_growth_38/lib/python3.8/site-packages/Reanimate/Build_Data/"
    network.setBuildPath(build_path,true);

    // Load network file
    printf("> Loading network...\n");
    network.graphOverride = false;
    network.loadNetwork(file_name);

    // Set stack size - used for functions w/ recursive loops (large networks)
    printf("> Setting stack size ...\n");
    network.setStackSize();

    // Solve for blood flow w/ variable viscosity and phase separation
    network.bloodFlow(true, false);
    //network.bloodFlow(false, false);
    
    // Output 2D visuals (false -> true if Amira file needed)
    network.printVisuals(false);

}
