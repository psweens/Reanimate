//
//  Network.hpp
//  TEEPEE
//
//  Created by Paul Sweeney on 16/06/2020.
//  Copyright © 2020 Paul Sweeney. All rights reserved.
//

#include <string>
#include <armadillo>

using namespace arma;
using namespace std;

#ifndef Network_hpp
#define Network_hpp

namespace reanimate {

    class Network {

    public:

        string networkName,networkPath,buildPath,loadPath,rLog;
        bool unknownBCs,phaseseparation;
        int mxx{},myy{},mzz{},nodsegm{},nsol{},nnodfl{};
        float alx{},aly{},alz{},lb{},maxl{};
        double targPress{},targStress{},tissperfusion{},inflow{},lthresh{10.};
        ivec ista,iend,segname,vesstyp,nodname,bcnodname,bctyp,nodtyp,bcnod,BCgeo,noflow,edgeLabels,nodout,nodrank,nk,flag,deadends,subGraphs;
        vec diam,rseg,lseg,q,qq,hd,bcprfl,bchd,nodpress,BCflow,BCpress,tau,segpress,elseg,ediam;
        imat segnodname,nodnod,nodseg;
        mat cnode,bcp;

        void setBuildPath();
        void loadNetwork(const string &filename, const bool directFromAmira=false);
        void analyse_network(bool graph = false);
        void subNetwork(ivec &index, bool graph = false);
        void edgeNetwork();
        void pictureNetwork(const string &filename, vec vector, bool logdist = false, int nl=20, bool nodes=false, bool segs=false);
        void fullSolver();
        void estimationSolver();
        void putrank(Network &graph);
        void printNetwork(const string& filename, bool resetDomain = false);
        void printHistogram(const string &filename, mat &data, const field<string> &headers);
        void printReducedAmira(const string &filename);
        void printAmira(const string &filename, const mat &extraData=zeros<mat>(0,0));
        void printNamira(const string &filename, const string &networkname);
        int readAmira(const string &filename, const string &networkname, bool stubs=false);
        void processAmira(const bool &stubs);


        // 'Getter' functions
        int getNseg();
        int getNnod();
        int getNnodbc();

        // 'Setter' functions
        void setNseg(int nseg);
        void setNnod(int nnod);
        void setNnodbc(int nnodbc);

        // Auxiliary functions
        int detect_col(FILE *ifp);
        void initLog();
        void printText(const string &text, const int type=2, const int newline=1);
        void printNum(const string &text, const double &num, const string unit="");
        void printStat(const string &text, const vec &n, const string &unit);

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

        void setup_networkArrays();
        void setup_flowArrays();
        void setup_estimationArrays();

    private:
        const char* FindAndJump(const char* buffer, const char* SearchString);

    };

}



#endif /* Network_hpp */
