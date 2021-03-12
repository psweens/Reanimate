//
//  Network.hpp
//  TEEPEE
//
//  Created by Paul Sweeney on 16/06/2020.
//  Copyright © 2020 Paul Sweeney. All rights reserved.
//

#include <string>
#include <armadillo>
#include <sys/resource.h>

using namespace arma;
using namespace std;

#ifndef Network_hpp
#define Network_hpp

namespace reanimate {

    class Network {

    public:

        string networkName,networkPath,buildPath,loadPath,rLog;
        bool unknownBCs,phaseseparation;
        int mxx{},myy{},mzz{},nodsegm{},nsol{},nnodfl{},track,nitmax{};
        double alx{},aly{},alz{},lb{},maxl{},targPress{},targStress{},tissperfusion{},inflow{},lthresh{10.},tissDensity{},bloodDensity{};
        ivec ista,iend,segname,vesstyp,nodname,bcnodname,bctyp,nodtyp,bcnod,BCgeo,noflow,edgeLabels,nodout,nodrank,nk,deadEnds,articPnt;
        vec diam,rseg,lseg,q,qq,vel,hd,bcprfl,bchd,nodpress,BCflow,BCpress,tau,segpress,elseg,ediam;
        uvec unknownnod_idx,bcpress_idx;
        imat segnodname,nodnod,nodseg;
        mat cnode,bcp;

        void setBuildPath();
        void loadNetwork(const string &filename, const bool directFromAmira=false);
        void analyse_network(bool graph = false, bool print = true);
        void indexSegmentConnectivity();
        void indexNodeConnectivity();
        void subNetwork(ivec &index, bool graph = false, bool print = true);
        void edgeNetwork();
        void pictureNetwork(const string &filename, vec vector, bool logdist = false, int nl=20, bool nodes=false, bool segs=false);
        void fullSolver();
        void estimationSolver();
        void putrank(Network &sGraph);
        int readAmira(const string &filename, const string &networkname, bool stubs=false);
        void processAmira(const bool &stubs);
        void dfsBasic(int v, int tag, ivec &track, bool nodeCondition=false, int type=2);
        void doubleDfs(int v, int tag, double val, ivec &track, vec &param, string cond="GT");
        void dfsArtic(int v, int p=-1);
        void findArticulationPoints();
        ivec findDeadends();
        void removeNewBC(ivec storeBCnodname, bool print=false, bool graph=false);


        // 'Getter' functions
        int getNseg();
        int getNnod();
        int getNnodbc();

        // 'Setter' functions
        void setNseg(int nseg);
        void setNnod(int nnod);
        void setNnodbc(int nnodbc);
        void setStackSize(int stackSize=16*1024*1024);

        // Auxiliary functions
        int detect_col(FILE *ifp);
        void initLog();
        double nodeAverage(const int &inod, const vec &param);

        // Print Functions
        void printText(const string &text, const int type=2, const int newline=1);
        void printNum(const string &text, const double &num, const string unit="");
        void printStat(const string &text, const vec &n, const string &unit);
        void printNetwork(const string& filename, bool resetDomain = false);
        void printHistogram(string filename, mat &data, const field<string> &headers);
        void printRawData(string filename, mat &data, const field<string> &headers);
        void printReducedAmira(const string &filename);
        void printAmira(const string &filename, const mat &extraData=zeros<mat>(0,0), bool smooth=true, const char *headers[]={});
        void printNamira(const string &filename, const string &networkname);

        Network(); // Constructor
        ~Network(); // Destructor

    protected:

        // Dimensional constants - millimetre scaling
        double gamma,alpha,beta,xi;

        // Network
        int nseg{},nnod{},nnodbc{},computeLseg{};

        // Amira variables
        int nvertex{},nedge{},npoint{};
        ivec nedgePoints;
        vec thickness;
        imat edgeConnectivity;
        mat vertexCoordinates,edgePointCoordinates;

        // Flow parameters
        int estimationarraysize{},nIBnod{},nunknown{};
        double constvisc{},consthd{},mcv{},hdtol{},qtol{},ktau{},oldktau{},kp{},targetpress{};
        ivec unknownnod,storeBCtyp;
        vec conductance,c,qold,hdold,flowsign,oldFlowsign,tau0,oldTau,Qo,B,p0,storeBC,storeBChd,storeHD,oldHd,oldNodpress,oldq;
        sp_mat M,L,K,A,H1,H2,W;

        // Topology fn parameters
        double timer;
        ivec visited, tin, low;

        void setup_networkArrays();
        void setup_flowArrays(bool popMatrices=true);
        void setup_estimationArrays();
        double pointAverage(const int &pnt, const ivec &pntIdx, const vec &param);

    private:

        const char* FindAndJump(const char* buffer, const char* SearchString);

    };

}



#endif /* Network_hpp */
