#include "Network.hpp"
#include "spatGraph.hpp"
#include "Vasculature.hpp"
#include "DiscreteContinuum.hpp"
#include "MicroCell.hpp"
#include "omp.h"

using namespace reanimate;
using namespace std;

void example_MicroFluidic_BloodFlow() {

    // Set number of threads for any multithreading
    omp_set_num_threads(48);

    // Define vasculature
    Vasculature network;

    // Set load and build paths
    network.setLoadPath("/Users/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/");
    network.setBuildPath("/Users/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/",true);

    // Load network file - 'cuboidVess=True' indicates network file with width/height columns
    network.loadNetwork("Tube.txt", true);

    // Iterate through various blood viscosity values (cP)
    vec viscRange = zeros<vec>(3);
    viscRange(0) = 2.;
    viscRange(1) = 3.;
    viscRange(2) = 4.;
    for (int i = 0; i < (int) viscRange.n_elem; i++)    {

        // Create output path for each value
        network.setBuildPath("/Users/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/Visc_"+to_string(i)+"/",true);

        // Set blood viscosity value
        network.constvisc = viscRange(i);

        // Solve for blood flow w/ variable viscosity and phase separation
        network.bloodFlow(true, false);

        // Output 2D visuals
        network.printVisuals(false);
    }

    // Reset blood viscosity to default value
    network.constvisc = 3.;

    // Iterate through various inflow conditions (nl/min - can set pressures instead)
    vec inflowRange = zeros<vec>(3);
    inflowRange(0) = 1.;
    inflowRange(1) = 10.;
    inflowRange(2) = 100.;
    network.bctyp(0) = 1.; // Setting inflow boundary node to flow condition (0 value = pressure) - index 0 is inflow node
    for (int i = 0; i < (int) inflowRange.n_elem; i++)  {

        // Create output path for each value
        network.setBuildPath("/Users/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/Inflow_"+to_string(i)+"/",true);

        // Set inflow value
        network.bcprfl(0) = inflowRange(i);

        // Solve for blood flow w/ variable viscosity and phase separation
        network.bloodFlow(true, false);

        // Output 2D visuals
        network.printVisuals(false);

    }

}
