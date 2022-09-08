#include "Network.hpp"
#include "spatGraph.hpp"
#include "Vasculature.hpp"
#include "DiscreteContinuum.hpp"
#include "MicroCell.hpp"
#include "omp.h"
#include "math.h"

using namespace reanimate;
using namespace std;

void example_ExtractSubnetwork()    {

    Vasculature vessNet;
    vessNet.buildPath = "/Users/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/";
    vessNet.loadPath = "/Users/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/";
    vessNet.setBuildPath(true);
    vessNet.loadNetwork("LS3Network.txt");
    vessNet.setStackSize();
    vessNet.loadDeadEnds = false;

/*    vessNet.bloodFlow(true);
    vessNet.computeSegpress();
    vessNet.computeBoundaryFlow();
    vessNet.printVisuals(false);*/

    // Generate network graph
    spatGraph graph;
    graph.generate(vessNet, true);

    // Initialise starting point
    int nod{};
    int depth{6};
    for (int inodbc = 0; inodbc < graph.getNnodbc(); inodbc++)   {
        if (graph.bcnodname(inodbc) == 0)    {
            nod = graph.bcnod(inodbc);
            inodbc = graph.getNnodbc();
        }
    }
    ivec order = graph.breadthFirstSearch(nod);

    // Calculate approx. edge branching order
    vec border = zeros<vec>(graph.getNseg());
    for (int iseg = 0; iseg < graph.getNseg(); iseg++) {
        border(iseg) = 0.5*floor(order(graph.ista(iseg)) + order(graph.iend(iseg)));
    }

    // Map branching order to full network
    ivec remove = zeros<ivec>(vessNet.getNseg());
    vec neworder = zeros<vec>(vessNet.getNseg());
    for (int iseg = 0; iseg < vessNet.getNseg(); iseg++)   {
        for (int jseg = 0; jseg < graph.getNseg(); jseg++)   {
            if (vessNet.edgeLabels(iseg) == graph.segname(jseg))  {
                neworder(iseg) = border(jseg);
                jseg = graph.getNseg();
            }
        }
    }

    // Copy across boundary conditions
    for (int inodbc = 0; inodbc < graph.getNnodbc(); inodbc++)  {
        for (int inod = 0; inod < vessNet.getNnod(); inod++)    {
            if (graph.bcnodname(inodbc) == vessNet.nodname(inod))   {
                graph.bctyp(inodbc) = 0;
                graph.bcprfl(inodbc) = vessNet.nodpress(inod);
                inod = vessNet.getNnod();
            }
        }
    }
    for (int inod = 0; inod < graph.getNnod(); inod++)  {
        for (int jnod = 0; jnod < vessNet.getNnod(); jnod++)    {
            if (graph.nodname(inod) == vessNet.nodname(jnod))   {
                graph.nodpress(inod) = vessNet.nodpress(jnod);
                jnod = vessNet.getNnod();
            }
        }
    }

    remove(find(neworder > depth)).fill(1);
    vec tmp = neworder(find(neworder <= depth));
    vessNet.subNetwork(remove);
    vessNet.pictureNetwork("NetworkOrder.ps", tmp);
    vessNet.printNetwork("LS3_Depth"+to_string(depth)+"_Tortuous.txt");

    remove = ones<ivec>(graph.getNseg());
    for (int iseg = 0; iseg < graph.getNseg(); iseg++)  {
        if (order(graph.ista(iseg)) <= depth)   {remove(iseg) = 0;}
        else if (order(graph.iend(iseg)) <= depth)   {remove(iseg) = 0;}
    }

    // remove(find(border > depth)).fill(1);
    tmp = border(find(remove == 0));
    graph.subNetwork(remove, true);
    graph.pictureNetwork("GraphOrder.ps", tmp);
    graph.printNetwork("LS3_Depth"+to_string(depth)+"_Graph.txt");

}