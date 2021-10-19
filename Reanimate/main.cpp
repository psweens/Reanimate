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

    example_Discrete_BloodFlow();

    // Find shortest path
    Vasculature test;
    test.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/";
    test.loadPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/";
    test.setBuildPath(true);
    test.loadNetwork("1Network.dat");
    test.setStackSize();
    //test.loadDeadEnds = true;
    test.bloodFlow(true);
    test.printVisuals(false);

    test.graph.generate(test, true);

    test.rheolParams();
    //test.hd.fill(test.consthd);
    test.computeConductance();

    vec tmp = zeros<vec>(test.graph.getNseg());
    for (int iseg = 0; iseg < test.graph.getNseg(); iseg++) {
        uvec idx = find(test.graph.segname(iseg) == test.edgeLabels);
        tmp(iseg) = mean(test.conductance(idx));
    }
    tmp /= max(tmp);
    tmp = 1./ tmp; // Normalised flow resistance

    //cout<<test.hd<<endl;

    vec printPath = zeros<vec>(test.graph.getNseg());
    int found{},nod1{},nod2{};
    ivec path;
    // omp parallel for schedule(auto) default(none) shared(tmp,printPath) private(path,test,nod1,nod2,found)
    for (int inodbc = 0; inodbc < test.graph.getNnodbc(); inodbc++)   {
        for (int jnodbc = inodbc+1; jnodbc < test.graph.getNnodbc(); jnodbc++)   {
            path = test.graph.findShortestPath(test.graph.bcnod(inodbc), test.graph.bcnod(jnodbc), tmp);
            for (int iseg = 0; iseg < test.graph.getNseg(); iseg++) {
                found = 0;
                nod1 = test.graph.ista(iseg);
                nod2 = test.graph.iend(iseg);
                for (int inod = 0; inod < (int) path.n_elem; inod++)    {
                    if (nod1 == path(inod)) {found += 1;}
                    else if (nod2 == path(inod))    {found += 1;}
                    if (found == 2) {
                        printPath(iseg) = 1.;
                        inod = (int) path.n_elem;
                    }
                }
            }
        }
    }
    //cout<<"here"<<endl;
    //cout<<test.hd(find(printPath == 1.))<<endl;
    test.graph.pictureNetwork("ShortestPath.ps", printPath);


/*    DiscreteContinuum oNetwork;
    oNetwork.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/";
    oNetwork.loadPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/";
    oNetwork.discreteNet.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/";
    oNetwork.discreteNet.loadPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/";
    oNetwork.discreteNet.setBuildPath(true);
    oNetwork.discreteNet.loadNetwork("SW1/SW1_1.txt");
    oNetwork.discreteNet.setStackSize();
    oNetwork.discreteNet.loadDeadEnds = true;

    // Extract subnetwork
    oNetwork.graph.generate(oNetwork.discreteNet, true);
    int nod{};
    int depth{3};
    for (int inodbc = 0; inodbc < oNetwork.graph.getNnodbc(); inodbc++)   {
        if (oNetwork.graph.bcnodname(inodbc) == 108)    {
            nod = oNetwork.graph.bcnod(inodbc);
            inodbc = oNetwork.graph.getNnodbc();
        }
    }
    ivec order = oNetwork.graph.breadthFirstSearch(nod);
    vec border = zeros<vec>(oNetwork.graph.getNseg());
    for (int iseg = 0; iseg < oNetwork.graph.getNseg(); iseg++) {
        border(iseg) = ceil(0.5*(order(oNetwork.graph.ista(iseg)) + order(oNetwork.graph.iend(iseg))));
    }
    ivec remove = zeros<ivec>(oNetwork.discreteNet.getNseg());
    vec neworder = zeros<vec>(oNetwork.discreteNet.getNseg());
    for (int iseg = 0; iseg < oNetwork.discreteNet.getNseg(); iseg++)   {
        for (int jseg = 0; jseg < oNetwork.graph.getNseg(); jseg++)   {
            if (oNetwork.discreteNet.edgeLabels(iseg) == oNetwork.graph.segname(jseg))  {
                if (border(jseg) > depth) {remove(iseg) = 1;}
                neworder(iseg) = border(jseg);
                jseg = oNetwork.graph.getNseg();
            }
        }
    }

    remove(find(neworder > depth)).fill(1);
    neworder = neworder(find(neworder <= depth));
    oNetwork.discreteNet.subNetwork(remove);
    oNetwork.discreteNet.pictureNetwork("NetworkOrder.ps", neworder);
    oNetwork.discreteNet.printNetwork("SW_Depth4_Tortuous.txt");

    remove = zeros<ivec>(oNetwork.graph.getNseg());
    remove(find(border > depth)).fill(1);
    border = border(find(border <= depth));
    oNetwork.graph.subNetwork(remove, true);
    oNetwork.graph.pictureNetwork("GraphOrder.ps", conv_to<vec>::from(border));
    oNetwork.graph.printNetwork("SW_Depth4_Graph.txt");*/

    // Medulla network
    // Generate network from Amira file
/*    Network medulla;
    medulla.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/";
    medulla.loadPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/";
    medulla.setBuildPath(false);
    medulla.lthresh = 10.;
    medulla.readAmira("GPUdeconvolved_40_iterations_GL1200_substack_div5-SptGraph.am", "MedullaLRG");*/

    // Base data
    FILE *data;
    data = fopen("/home/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/Medulla_Trees.txt","r");
    int n = 32;
    imat array = zeros<imat>(n,2);
    for (int i = 0; i < n; i++) {
        fscanf(data, "%lli %lli", &array(i,0), &array(i,1));
    }
    fclose(data);

    DiscreteContinuum hybrid;
    hybrid.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/";
    hybrid.loadPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/";
    hybrid.discreteNet.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/";
    hybrid.discreteNet.loadPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/";
    hybrid.discreteNet.setBuildPath(true);
    hybrid.discreteNet.loadNetwork("MedullaSolved.txt");
    hybrid.discreteNet.setStackSize();
    hybrid.discreteNet.loadDeadEnds = true;

/*    ivec flag = zeros<ivec>(hybrid.discreteNet.getNseg());
    flag(find(abs(hybrid.discreteNet.q) < 2.)).fill(1);gsuit
    hybrid.discreteNet.subNetwork(flag);*/

/*    const char *headers[1] = {"Vesstyp"};
    mat data3 = zeros<mat>(hybrid.discreteNet.getNseg(), 1);
    data3.col(0) = conv_to<vec>::from(hybrid.discreteNet.vesstyp);
    hybrid.discreteNet.printAmira("test.am", data3, true, headers);

    ivec flag = zeros<ivec>(hybrid.discreteNet.getNseg());
    flag(find(hybrid.discreteNet.vesstyp == 2)).fill(1);
    hybrid.discreteNet.subNetwork(flag);
    hybrid.discreteNet.printAmira("test2.am");*/

//    spatGraph test;
//    test.generate(hybrid.discreteNet);

//    const char *headers[5] = {"Pressure","Flow","Hd","Velocity","WSS"};
//    mat data1 = zeros<mat>(hybrid.discreteNet.getNseg(), 5);
//    hybrid.discreteNet.printAmira("test.am", data1, true, headers);

/*    for (int inod = 0; inod < n; inod++)    {
        uvec jdx = find(hybrid.discreteNet.bcnodname == array(inod, 0));
        hybrid.discreteNet.bctyp(jdx(0)) = 0;
        hybrid.discreteNet.BCgeo(jdx(0)) = array(inod,1);
    }
    hybrid.discreteNet.cortexBoundaryPress();*/


/*    vec diamAvg = zeros<vec>(hybrid.discreteNet.getNnod());
    for (int inod = 0; inod < hybrid.discreteNet.getNnod(); inod++) {
        diamAvg(inod) = hybrid.discreteNet.nodeAverage(inod, hybrid.discreteNet.diam);
    }

    // Initial branches
    ivec track = -ones<ivec>(hybrid.discreteNet.getNnod());
    for (int inodbc = 0; inodbc < hybrid.discreteNet.getNnodbc(); inodbc++) {
        if (diamAvg(hybrid.discreteNet.bcnod(inodbc)) > 10. &&
        (hybrid.discreteNet.cnode(2,hybrid.discreteNet.bcnod(inodbc)) > 0.9*hybrid.discreteNet.alz
        || hybrid.discreteNet.cnode(2,hybrid.discreteNet.bcnod(inodbc)) < 0.1*hybrid.discreteNet.alz))    {
            hybrid.discreteNet.doubleDfs(hybrid.discreteNet.bcnod(inodbc), 1.0, 10., track, diamAvg, "GT");
        }
    }
    vec flag1 = zeros<vec>(hybrid.discreteNet.getNseg());
    for (int iseg = 0; iseg < hybrid.discreteNet.getNseg(); iseg++) {
        if (track(hybrid.discreteNet.ista(iseg)) == 1 || track(hybrid.discreteNet.iend(iseg)) == 1) {
            flag1(iseg) = 1;
        }
    }

    ivec remove = zeros<ivec>(hybrid.discreteNet.getNseg());
    remove(find(flag1 == 0)).fill(1);

    hybrid.discreteNet.subNetwork(remove);


    // Get rid of noise
    track = -ones<ivec>(hybrid.discreteNet.getNnod());

    remove = zeros<ivec>(hybrid.discreteNet.getNseg());
    for (int inodbc = 0; inodbc < hybrid.discreteNet.getNnodbc(); inodbc++) {
        track.ones();
        track *= -1;
        hybrid.discreteNet.dfsBasic(hybrid.discreteNet.bcnod(inodbc), inodbc, track);
        ivec nod = zeros<ivec>(hybrid.discreteNet.getNnod());
        nod(find(track > -1)).fill(1);
        if (accu(nod) <= 10) {
            for (int iseg = 0; iseg < hybrid.discreteNet.getNseg(); iseg++) {
                if (nod(hybrid.discreteNet.ista(iseg)) == 1 && nod(hybrid.discreteNet.iend(iseg)) == 1) {remove(iseg) = 1;}
            }
        }
    }
    hybrid.discreteNet.subNetwork(remove);

    //hybrid.discreteNet.pictureNetwork("Network_Diameter.ps",hybrid.discreteNet.diam);
    flag1 = zeros<vec>(hybrid.discreteNet.getNseg());
    for (int inod = 0; inod < n; inod++)    {
        for (int iseg = 0; iseg < hybrid.discreteNet.getNseg(); iseg++) {
            if (array(inod,0) == hybrid.discreteNet.nodname(hybrid.discreteNet.ista(iseg)) || array(inod,0) == hybrid.discreteNet.nodname(hybrid.discreteNet.iend(iseg)))   {
                flag1(iseg) = 1;
            }
        }
    }
    mat extraD = zeros<mat>(hybrid.discreteNet.getNseg(), 1);
    extraD.col(0) = flag1;*/

/*    hybrid.discreteNet.bloodFlow(true, false);
    hybrid.discreteNet.printNetwork("SolvedNetwork.txt");
    hybrid.discreteNet.printVisuals(true, true);*/

    hybrid.discreteNet.findBoundingBox();

    // Cortex vessel classification
    ivec track = -ones<ivec>(hybrid.discreteNet.getNnod());
    for (int i = 0; i < (int) array.n_rows; i++) {
        for (int inodbc = 0; inodbc < hybrid.discreteNet.getNnodbc(); inodbc++) {
            if (array(i,0) == hybrid.discreteNet.bcnodname(inodbc)) {
                hybrid.discreteNet.branch = 0;
                int tag = array(i,1);
                hybrid.discreteNet.dfsBranch(hybrid.discreteNet.bcnod(inodbc), tag, track, 50);
                inodbc = hybrid.discreteNet.getNnodbc();
            }
        }
    }
    hybrid.discreteNet.vesstyp.fill(2);
    for (int iseg = 0; iseg < hybrid.discreteNet.getNseg(); iseg++)    {
        if (track(hybrid.discreteNet.ista(iseg)) > 0)   {hybrid.discreteNet.vesstyp(iseg) = track(hybrid.discreteNet.ista(iseg));}
        else if (track(hybrid.discreteNet.iend(iseg)) > 0)   {hybrid.discreteNet.vesstyp(iseg) = track(hybrid.discreteNet.iend(iseg));}
    }

/*    const char *headers[1] = {"Vesstyp"};
    mat data3 = zeros<mat>(hybrid.discreteNet.getNseg(), 1);
    data3.col(0) = conv_to<vec>::from(hybrid.discreteNet.vesstyp);
    hybrid.discreteNet.printAmira("amiraVesstype.am", data3, true, headers);*/

    // Run vessel classification
    // Mes. 1
    imat inOutlets = zeros<imat>(2,2);
/*    inOutlets(0,0) = 830;
    inOutlets(0,1) = 1;
    inOutlets(1,0) = 825;
    inOutlets(1,1) = 2;*/

    inOutlets = array;
    for (int i = 0; i < n; i++) {
        if (array(i,1) == 3)    {array(i,1) = 2;}
    }


    // Mes. 3
    /*inOutlets = zeros<imat>(2,2);
    inOutlets(0,0) = 78;
    inOutlets(0,1) = 1;
    inOutlets(1,0) = 76;
    inOutlets(1,1) = 2;*/

    // Mes. 5
    /*imat inOutlets = zeros<imat>(14,2);
    inOutlets(0,0) = 999;
    inOutlets(0,1) = 1;
    inOutlets(1,0) = 2996;
    inOutlets(1,1) = 2;
    inOutlets(2,0) = 3;
    inOutlets(2,1) = 1;
    inOutlets(3,0) = 994;
    inOutlets(3,1) = 2;
    inOutlets(4,0) = 19;
    inOutlets(4,1) = 2;
    inOutlets(5,0) = 992;
    inOutlets(5,1) = 2;
    inOutlets(6,0) = 1109;
    inOutlets(6,1) = 2;
    inOutlets(7,0) = 1102;
    inOutlets(7,1) = 2;
    inOutlets(8,0) = 579;
    inOutlets(8,1) = 1;
    inOutlets(9,0) = 984;
    inOutlets(9,1) = 2;
    inOutlets(10,0) = 992;
    inOutlets(10,1) = 1;
    inOutlets(11,0) = 550;
    inOutlets(11,1) = 2;
    inOutlets(12,0) = 456;
    inOutlets(12,1) = 1;
    inOutlets(13,0) = 297;
    inOutlets(13,1) = 2;*/

    // Cortex
    /*imat inOutlets = zeros<imat>(4,2);
    inOutlets(0,0) = 601;
    inOutlets(0,1) = 1;
    inOutlets(1,0) = 599;
    inOutlets(1,1) = 1;
    inOutlets(2,0) = 595;
    inOutlets(2,1) = 2;
    inOutlets(3,0) = 576;
    inOutlets(3,1) = 2;*/


    // Generate spatial graph
    hybrid.graph.generate(hybrid.discreteNet, true);
    hybrid.graph.findTree(inOutlets);
//    hybrid.graph.analyseTopology(inOutlets, hybrid.discreteNet);
    hybrid.discreteNet.analyseBoundaryType();

    // Compute boundary flow if flow solver has not run
    hybrid.discreteNet.computeBoundaryFlow();

    // Generate micro-cell
    hybrid.cell.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/";
    hybrid.cell.loadPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/";
    hybrid.cell.setBuildPath(false); // Create build directory / delete contents

    hybrid.graph.findLengths();
    hybrid.cell.setEdgeDiamDistrib(hybrid.graph.diam(find(hybrid.graph.vesstyp == 2)));
    hybrid.cell.setEdgeLengthDistrib(hybrid.graph.lseg(find(hybrid.graph.vesstyp == 2)));
    //hybrid.cell.rotationAngle = M_PI / 3.;
    hybrid.cell.computeConductivity("crossCell3D");


    // Run discrete-continuum model
    hybrid.runHybrid();



}
