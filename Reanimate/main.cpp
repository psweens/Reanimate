#include "Network.hpp"
#include "spatGraph.hpp"
#include "Vasculature.hpp"
#include "DiscreteContinuum.hpp"
#include "MicroCell.hpp"
#include "Examples/Examples.h"
#//include "Working/Working_Scripts.h"
#include "omp.h"

using namespace reanimate;
using namespace std;

int main(int argc, char** argv) {

    std::string file_path;
    std::string build_path;
    if (argc >= 2) {
        file_path = argv[1];
        build_path = argv[2];
        example_Discrete_BloodFlow(file_path,build_path);
    }
    else {
        file_path = "/mnt/data2/retinasim/cco/graph/";
        build_path = "/mnt/ml/anaconda_envs/vessel_growth_38/lib/python3.8/site-packages/Reanimate/Build_Data/";
        example_Discrete_BloodFlow(file_path,build_path);
    }

    // Finding tumour redundancy
    //tumour_ShortestPath();

    // Extract tumour subnetwork
    //example_ExtractSubnetwork();

    // Tumour haematocrit simulations
    //tumour_PhaseSeparation();

    //example_Discrete_BloodFlow();
    //example_Hybrid_BloodFlow();
    //Medulla();

}
