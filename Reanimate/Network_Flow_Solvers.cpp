#include "Network.hpp"

using namespace reanimate;

void Network::fullSolver()    {

    // Construct conductance matrix, M
    if (!phaseseparation)   {printText("Computing conductance matrix, M",2, 0);}
    double tcond{};
    for (int iseg = 0; iseg < nseg; iseg++) {
        tcond = conductance(iseg);
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
        else    {Qo(idx) = -bcprfl(inodbc)*gamma;}
    }


    superlu_opts settings;
    if (phaseseparation)    {
        //settings.allow_ugly = true;
        settings.equilibrate = true;
        settings.refine = settings.REF_EXTRA;
    }

    if (!phaseseparation)   {printText("Solving linear system",2, 0);}
    nodpress = spsolve(K,Qo, "superlu", settings)/alpha;
    q = (M * nodpress)*(alpha/gamma);
    tau = c % abs(q) * (gamma/beta);

}


void Network::estimationSolver() {

    // Construct conductance, M
    double tcond{};
    for (int iseg = 0; iseg < nseg; iseg++) {
        tcond = conductance(iseg);
        M(iseg,ista(iseg)) = tcond;
        M(iseg,iend(iseg)) = -tcond;
    }

    // Construct K matrix but inputting into matrix A to reduce run-time
    A(span(0,nIBnod-1),span(0,nnod-1)) = L * M;
    int nod1{},nod2{};
    uvec idx{};
    rowvec aRow = zeros<rowvec>(A.n_cols);
    for (int inodbc = 0; inodbc < (int) bcpress_idx.n_elem; inodbc++) {
        nod1 = bcnod(bcpress_idx(inodbc)); // Node idx w/ pressure BC
        idx = find(nod1 == unknownnod_idx);
        nod2 = idx(0);  // Node idx of known nodes
        aRow = A.row(nod2);
        aRow(find(aRow > 0.)).fill(0.);
        aRow(nod1) = 1.;
        A.row(nod2) = aRow;
    }

    // Constructing H matrix - replacing with matrix A for reduced run-time
    if (nseg < 1e3) {
        A(span(nIBnod,estimationarraysize-1),span(0,nnod-1)) = ((repmat((square(c).t()),nnod,1) % (M.t()))*(repmat(lseg,1,nnod) % M))*ktau + W;
    }
    else    {
        double tlseg{}, tc{}, tMsta{}, tMend{};
        long tsta{}, tend{};
        for (int iseg = 0; iseg < nseg; iseg++) {
            tlseg = lseg(iseg);
            tc = pow(c(iseg),2);
            tsta = ista(iseg);
            tend = iend(iseg);
            tMsta = M(iseg,tsta);
            tMend = M(iseg,tend);
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
    B(span(nIBnod,estimationarraysize-1)) = (W.diag() % p0) + (ktau * (M.t() * (tau0 % (c % lseg))));


    // Equilibrates system - scales rows and columns to unit norm in order to prevent singular system
    superlu_opts settings;
    //settings.allow_ugly = true;
    settings.equilibrate = true;
    settings.refine = settings.REF_EXTRA;
    vec x = spsolve(A,B,"superlu",settings);


    nodpress = x(span(0,nnod-1));
    q = (M * nodpress);
    tau = c % abs(q);
}