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
        bool unknownBCs,phaseseparation,silence;
        int mxx{},myy{},mzz{},nodsegm{},nsol{},nnodfl{},track{},nitmax{100},nvertex{},nedge{},npoint{},branch{};
        double alx{},aly{},alz{},lb{},maxl{},targPress{},targStress{},tissperfusion{},inflow{},lthresh{10.},tissDensity{},bloodDensity{},consthd{},constvisc{};
        ivec ista,iend,segname,vesstyp,nodname,bcnodname,bctyp,nodtyp,bcnod,BCgeo,noflow,edgeLabels,flagTree,nodout,nodrank,nk,flag,deadends,subGraphs,loops,sgraphTag,ngraphTag,deadEnds,articPnt,edgePnts;
        vec diam,rseg,lseg,q,qq,vel,hd,bcprfl,bchd,nodpress,conductance,BCflow,BCpress,tau,segpress,elseg,ediam;
        uvec unknownnod_idx,bcpress_idx;
        imat segnodname,nodnod,nodseg,segpoints;
        mat cnode,bcp;

        virtual double evalTissPress(vec &x) {return 0.;}
        virtual double evalTissVel(vec &x) {return 0.;}

        void setBuildPath(bool deleteFiles=true);
        void loadNetwork(const string &filename, const bool directFromAmira=false);
        void analyse_network(bool graph = false, bool print = true);
        void indexSegmentConnectivity();
        void indexNodeConnectivity();
        void indexBCconnectivity();
        void findLengths();
        void analyseBoundaryType();
        void scaleNetwork(double sfactor);
        void subNetwork(ivec &index, bool graph = false, bool print = true);
        void edgeNetwork();
        void pictureNetwork(const string &filename, vec vector, bool logdist = false, int nl=20, bool nodes=false, bool segs=false);
        void fullSolver();
        void estimationSolver();
        void putrank(Network &sGraph);
        void computeBoundaryFlow();
        int readAmira(const string &filename, const string &networkname, bool stubs=false);
        void processAmira(const bool &stubs);
        void dfsBasic(int v, int tag, ivec &track, bool nodeCondition=false, int type=2);
        void doubleDfs(int v, int tag, double val, ivec &track, vec &param, string cond="GT");
        void dfsBranch(int v, int tag, ivec &track, int maxBranch=10);
        void dfsArtic(int v, int p=-1);
        ivec breadthFirstSearch(int nod);
        ivec findShortestPath(int startnode, int endNode, vec &edgeWeight, bool printPaths=false);
        void printShortestPaths(int startNode, int endNode, vec distance, ivec pred);

        void findArticulationPoints();
        ivec findDeadends();
        void removeNewBC(ivec storeBCnodname, bool print=false, bool graph=false);
        void setup_networkArrays();
        void setup_estimationArrays();
        void setup_flowArrays(bool popMatrices=true);


        // 'Getter' functions
        int getNseg();
        int getNnod();
        int getNnodbc();

        // 'Setter' functions
        void setNseg(int nseg);
        void setNnod(int nnod);
        void setNnodbc(int nnodbc);
        void setStackSize(int stackSize=16*1024*1024);

        // Math functions
        double eucDistance(vec &x, vec &y);
        double MSE(vec x, vec y);
        vec MaxSE(vec x, vec y);
        double MedianSE(vec x, vec y);
        double SPHI(int n, double &x);
        double SPHK(int n, double &x);
        vec lognormal(double &mean, double &SD, int &nseg);
        mat eulerRotation(const vec &c, const vec &d);

        // Auxiliary functions
        int detect_col(FILE *ifp);
        void initLog();
        double nodeAverage(const int &inod, const vec &param);
        void findBoundingBox();

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
        void plotContour(const string filename, Network &graph, double maxval, double minval, bool vector=false, bool overlay=true, const int xgrid=1e2, const int ygrid=1e2, const int NL=10);
        void shadeContour(FILE *ofp, const int &m, const int &n, double &scalefac, int &nl, const double pint, const double &xmin, const double &xmax, const double &ymin, const double &ymax, const vec &cl, const mat &zv);

        Network(); // Constructor
        ~Network(); // Destructor

    protected:

        // Dimensional constants - millimetre scaling
        double gamma,alpha,beta,xi;

        // Network
        int nseg{},nnod{},nnodbc{},computeLseg{};

        // Amira variables
        ivec nedgePoints;
        vec thickness;
        imat edgeConnectivity;
        mat vertexCoordinates,edgePointCoordinates;

        // Flow parameters
        int estimationarraysize{},nIBnod{},nunknown{};
        double mcv{},hdtol{},qtol{},ktau{},oldktau{},kp{},targetpress{};
        ivec unknownnod,storeBCtyp;
        vec c,qold,hdold,flowsign,oldFlowsign,tau0,oldTau,Qo,B,p0,storeBC,storeBChd,storeHD,oldHd,oldNodpress,oldq;
        sp_mat M,L,K,A,H1,H2,W;

        // Topology fn parameters
        double timer{};
        ivec visited, tin, low;

        double pointAverage(const int &pnt, const ivec &pntIdx, const vec &param);

    private:

        void findBoundaryNodes();
        const char* FindAndJump(const char* buffer, const char* SearchString);

    };

}



#endif /* Network_hpp */
