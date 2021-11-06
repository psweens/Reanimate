#include "Network.hpp"
#include "spatGraph.hpp"
#include "Vasculature.hpp"
#include "DiscreteContinuum.hpp"
#include "MicroCell.hpp"
#include "Examples/Examples.h"
#include "omp.h"

using namespace reanimate;
using namespace std;

int main(int argc, char** argv) {

    omp_set_num_threads(48);

    for (int i = 1; i <= 1; i++)    {

        // Find shortest path
        Vasculature test;
        test.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/LS3_" + to_string(i) + "/";
        test.loadPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/LS3/";
        test.setBuildPath(true);
        string netName = "LS3_" + to_string(i) + ".txt";
        test.loadNetwork(netName);
        test.setStackSize();
        //test.loadDeadEnds = true;

/*    for (int inodbc = 0; inodbc < test.getNnodbc(); inodbc++)   {
        test.bcprfl(inodbc) = test.nodpress(test.bcnod(inodbc));
        test.bctyp(inodbc) = 0;
    }*/

        test.bloodFlow(true);
        test.printVisuals(false);

        test.graph.generate(test, true);

        test.rheolParams();
        //test.hd.fill(test.consthd);
        //test.computeConductance();
        test.conductance = zeros<vec>(test.getNseg());
        for (int iseg = 0; iseg < test.getNseg(); iseg++) {
            test.conductance(iseg) = M_PI * pow(test.diam(iseg), 4) / (128. * test.lseg(iseg));
        }

        vec tmp = zeros<vec>(test.graph.getNseg());
        for (int iseg = 0; iseg < test.graph.getNseg(); iseg++) {
            uvec idx = find(test.graph.segname(iseg) == test.edgeLabels);
            tmp(iseg) = mean(test.conductance(idx));
        }
        tmp /= max(tmp);
        tmp = 1./ tmp; // Normalised flow resistance

        test.graph.BCflow = zeros<vec>(test.graph.getNnodbc());
        test.graph.BCpress = zeros<vec>(test.graph.getNnodbc());
        for (int inodbc = 0; inodbc < test.getNnodbc(); inodbc++)   {
            for (int jnodbc = 0; jnodbc < test.graph.getNnodbc(); jnodbc++) {
                if (test.bcnodname(inodbc) == test.graph.bcnodname(jnodbc)) {
                    test.graph.BCflow(jnodbc) = test.BCflow(inodbc);
                    test.graph.BCflow(jnodbc) = test.BCflow(inodbc);
                }
            }
        }

        vec printPath = zeros<vec>(test.getNseg());
        ivec graphPaths = test.graph.findRedundant(tmp);
        graphPaths.save(test.buildPath + "LS3_"+to_string(i)+"_saveGraphRedundancy.txt", raw_ascii);

        test.graph.pictureNetwork("ShortestPath.ps", conv_to<vec>::from(graphPaths));

        const char *headers[1] = {"Redundancy"};
        mat data = zeros<mat>(test.graph.getNseg(), 1);
        data.col(0) = printPath;

        string filename = "LS3_" + to_string(i) + "_graph_redundancy.am";
        test.graph.printAmira(filename, data,false, headers);

        vec fullRedundantPaths = zeros<vec>(test.getNseg());
        for (int iseg = 0; iseg < test.graph.getNseg(); iseg++) {
            if (graphPaths(iseg) == 1) {
                for (int jseg = 0; jseg < test.getNseg(); jseg++) {
                    if (fullRedundantPaths(jseg) == 0) {
                        if (test.graph.segname(iseg) == test.edgeLabels(jseg))    {
                            fullRedundantPaths(jseg) = 1;
                        }
                    }
                }
            }
        }

        data = zeros<mat>(test.getNseg(), 1);
        data.col(0) = fullRedundantPaths;
        test.printAmira("LS3_"+to_string(i)+"_redundancy.am", data, true, headers);

        fullRedundantPaths.save(test.buildPath + "LS3_"+to_string(i)+"_saveRedundancy.txt", raw_ascii);

    }

}
