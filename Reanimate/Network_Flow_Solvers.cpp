#include "Network.hpp"

using namespace reanimate;

void Network::fullSolver()    {

    // Construct conductance matrix, M
    if (!phaseseparation)   {printText("Computing conductance matrix, M",2, 0);}
    double tcond{};
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (deadends(iseg) == 1)    {tcond = 1.0;}
        else {tcond = conductance(iseg);}
        //tcond = conductance(iseg);
        M(iseg,ista(iseg)) = tcond;
        M(iseg,iend(iseg)) = -tcond;
    }

    // Construct K matrix
    if (!phaseseparation)   {printText("Computing K matrix",2, 0);}
    K = L * M;

    // Define qo - pressure and flow boundary conditions
    if (!phaseseparation)   {printText("Assigning boundary conditions",2, 0);}
    int idx{};
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        idx = bcnod(inodbc);
        if (bctyp(inodbc) == 0)  {
            Qo(idx) = bcprfl(inodbc)*alpha;
            K.row(idx).zeros();
            K(idx,idx) = 1.;
        }
        else    {
            Qo(idx) = -bcprfl(inodbc)*gamma;
        }
    }

    superlu_opts settings;
    if (phaseseparation)    {
        //settings.allow_ugly = true;
        settings.equilibrate = true;
        settings.refine = settings.REF_EXTRA;
    }

    if (!phaseseparation)   {printText("Solving linear system",2, 0);}
    nodpress = spsolve(K,Qo)/alpha;
    q = (M * nodpress)*(alpha/gamma);

}


void Network::estimationSolver() {

    // Construct conductance, M
    double tcond{};
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (deadends(iseg) == 1)    {tcond = 1.0;}
        else {tcond = conductance(iseg);}
        //tcond = conductance(iseg);
        M(iseg,ista(iseg)) = tcond;
        M(iseg,iend(iseg)) = -tcond;
    }

    // Construct K matrix but inputting into matrix A to reduce run-time
    A(span(0,nIBnod-1),span(0,nnod-1)) = L * M;
    int cntr{};
    for (int inod = 0; inod < nnod; inod++) {
        if (unknownnod(inod) == 0)  {
            for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
                if (bcnod(inodbc) == inod && bctyp(inodbc) == 0)    {
                    A.row(inod - cntr).zeros();
                    A(inod - cntr,inod) = 1.;
                }
            }
        }
        else {
            cntr += 1;
        }
    }

    // Constructing H matrix - replacing with matrix A for reduced run-time
    if (nseg < 1e3) {
        A(span(nIBnod,estimationarraysize-1),span(0,nnod-1)) = ((repmat((square(c).t()),nnod,1) % (M.t()))*(repmat(lseg,1,nnod) % M))*ktau + W;
    }
    else    {

        for (int iseg = 0; iseg < nseg; iseg++) {
            double tlseg = lseg(iseg);
            double tc = pow(c(iseg),2);
            long tsta = ista(iseg);
            long tend = iend(iseg);
            double tMsta = M(iseg,tsta);
            double tMend = M(iseg,tend);
            H2(iseg,tsta) = tMsta * tlseg;
            H1(tsta,iseg) = tMsta * tc;
            H2(iseg,tend) = tMend * tlseg;
            H1(tend,iseg) = tMend * tc;
        }

        A(span(nIBnod,estimationarraysize-1),span(0,nnod-1)) = ktau*(H1 * H2) + W;
    }

    // Final A matrix entry - transpose of K
    A(span(nIBnod,estimationarraysize-1),span(nnod,estimationarraysize-1)) = A(span(0,nIBnod-1),span(0,nnod-1)).t();

    // Construct B vector
    B(span(nIBnod,estimationarraysize-1)) = (W*p0)+(ktau*(M.t()*(tau0 % (c % lseg))));


    vec x;
    // Equilibrates system - scales rows and columns to unit norm in order to prevent singular system
    superlu_opts settings;
    //settings.allow_ugly = true;
    settings.equilibrate = true;
    settings.refine = settings.REF_EXTRA	;
    x = spsolve(A,B,"superlu",settings);


    nodpress = x(span(0,nnod-1));
    q = (M * nodpress);
    tau = c % q;

}