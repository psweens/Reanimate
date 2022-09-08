#include "Network.hpp"
#include "spatGraph.hpp"
#include "Vasculature.hpp"
#include "DiscreteContinuum.hpp"
#include "MicroCell.hpp"
#include "Examples/Examples.h"
#include "Working/Working_Scripts.h"
#include "omp.h"

using namespace reanimate;
using namespace std;

int main(int argc, char** argv) {

    // Finding tumour redundancy
    //tumour_ShortestPath();

    // Extract tumour subnetwork
    //example_ExtractSubnetwork();

    // Tumour haematocrit simulations
    //tumour_PhaseSeparation();

    //example_Hybrid_BloodFlow();
    Medulla();

}
