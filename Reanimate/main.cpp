#include "Network.hpp"
#include "spatGraph.hpp"
#include "Vasculature.hpp"

using namespace reanimate;
using namespace std;

int main(int argc, char** argv) {

    // Directory setup
    Vasculature vessNetwork;
    vessNetwork.buildPath = "/home/sweene01/Documents/Reanimate/Build_Data/";
    vessNetwork.loadPath = "/home/sweene01/Documents/Reanimate/Load_Data/";
    vessNetwork.setBuildPath(); // Create build directory / delete contents

    // Read & process amira binary file
    vessNetwork.lthresh = 10; // Default value in microns (can delete line if using 10um)
    // Below creates network file ('networkname.txt') & new Amira file based on above threshold ('reducedAmira.am')
    vessNetwork.readAmira("LS3.am","LS3_Colorectal_Carcinoma");
    // To load network file straight away, change "LS3Network.txt" to "newNetwork.txt" -> need to assign BCs after loading

    // Load network file
    vessNetwork.loadNetwork("LS3_Colorectal_Carcinoma.txt",true);

    // Solving for blood flow
    // Following is an example of assigning pressure boundary conditions
    for (int inodbc = 0; inodbc < vessNetwork.getNnodbc(); inodbc++)    {
        vessNetwork.bcprfl(inodbc) = 90.;
        vessNetwork.bctyp(inodbc) = 0;
    }
    vessNetwork.bloodFlow(true, true);
    vessNetwork.pictureNetwork("NetworkDiameters.ps", vessNetwork.diam);
    vessNetwork.pictureNetwork("NetworkFlow.ps", vessNetwork.qq, true);
    vessNetwork.pictureNetwork("NetworkPressure.ps", vessNetwork.segpress);
    vessNetwork.pictureNetwork("NetworkShearStress.ps", vessNetwork.tau, true);
    vessNetwork.pictureNetwork("NetworkHD.ps", vessNetwork.hd);
    vessNetwork.printNetwork("solvedNetwork.txt");
    mat extraData = zeros<mat>(vessNetwork.getNseg(), 1);
    extraData.col(0) = vessNetwork.hd;
    vessNetwork.printAmira("amiraHD.am", extraData);

    // Vessel classification analysis
    //spatGraph eCortex;
    //eCortex.generate(vessNetwork, true);
    //eCortex.pictureNetwork(buildPath + "eNetworkDiameters.ps",eCortex.diam);
    //eCortex.loadTrunks(loadPath + "NetworkClass.txt");
    //eCortex.pictureNetwork(buildPath + "eNetworkGeometry.ps",conv_to<vec>::from(eCortex.geometry));*/

}
