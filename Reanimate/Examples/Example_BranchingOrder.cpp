#include "Network.hpp"
#include "spatGraph.hpp"
#include "Vasculature.hpp"
#include "DiscreteContinuum.hpp"
#include "MicroCell.hpp"
#include "omp.h"
#include "math.h"

using namespace std;
using namespace reanimate;

void example_BranchingOrder()    {

    Vasculature vessNet;
    vessNet.buildPath = "/Users/paul/Dropbox/Code/C++/Reanimate/Build_Data/";
    vessNet.loadPath = "/Users/paul/Dropbox/Code/C++/Reanimate/Load_Data/";
    vessNet.setBuildPath(true);
    vessNet.loadNetwork("LS3/LS3_1.txt");
    vessNet.setStackSize();
    vessNet.loadDeadEnds = true;

    // Extract subnetwork
    spatGraph graph;
    graph.generate(vessNet, true);
    int nod{};
    int depth{6};
    for (int inodbc = 0; inodbc < graph.getNnodbc(); inodbc++)   {
        if (graph.bcnodname(inodbc) == 108)    {
            nod = graph.bcnod(inodbc);
            inodbc = graph.getNnodbc();
        }
    }
    ivec order = graph.breadthFirstSearch(nod);
    vec border = zeros<vec>(graph.getNseg());
    for (int iseg = 0; iseg < graph.getNseg(); iseg++) {
        border(iseg) = ceil(0.5*(order(graph.ista(iseg)) + order(graph.iend(iseg))));
    }
    ivec remove = zeros<ivec>(vessNet.getNseg());
    vec neworder = zeros<vec>(vessNet.getNseg());
    for (int iseg = 0; iseg < vessNet.getNseg(); iseg++)   {
        for (int jseg = 0; jseg < graph.getNseg(); jseg++)   {
            if (vessNet.edgeLabels(iseg) == graph.segname(jseg))  {
                if (border(jseg) > depth) {remove(iseg) = 1;}
                neworder(iseg) = border(jseg);
                jseg = graph.getNseg();
            }
        }
    }

    remove(find(neworder > depth)).fill(1);
    neworder = neworder(find(neworder <= depth));
    vessNet.subNetwork(remove);
    vessNet.pictureNetwork("NetworkOrder.ps", neworder);
    vessNet.printNetwork("LS3_Depth6_Tortuous.txt");

    remove = zeros<ivec>(graph.getNseg());
    remove(find(border > depth)).fill(1);
    border = border(find(border <= depth));
    graph.subNetwork(remove, true);
    graph.pictureNetwork("GraphOrder.ps", conv_to<vec>::from(border));
    graph.printNetwork("LS3_Depth6_Graph.txt");

}