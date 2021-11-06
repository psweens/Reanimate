#include "Tissue.hpp"

using namespace reanimate;

/*double Tissue::bicgstab(vec &b, vec &x, int &n, double &eps, int &itmax)     {

    // Ax = b

    int isp = 0;

    int i{},j{},kk{},subn{n-1};
    double lu,lunew,beta,delta,gamma1,t1,t2,err=0.;
    vec r = zeros<vec>(n);
    vec rs = ones<vec>(n);
    vec v = zeros<vec>(n);
    vec s = zeros<vec>(n);
    vec t = zeros<vec>(n);
    vec p = zeros<vec>(n);
    vec er = zeros<vec>(n);


    double rstore=0.,vstore=0.,tstore=0.,rstore2=0.,vstore2=0.,tstore2=0.,tempStore=0.;
    double dist = 0.;
    double rr0;
    double fac = 0.;
    vec coord;

    double subnvar = x(subn);

#pragma omp parallel for default(none) shared(subn,isp) private(coord,rr0,fac,dist,i,j) firstprivate(subnvar,vseg,snode,x) schedule(dynamic) reduction(+:tempStore,rstore,rstore2)
    for(i = 0; i < subn; i++)  {
        rstore = 0.;
        rstore2 = 0.;
        coord = snode.col(i);
        rr0 = r0(i);
        fac = gfac(i,isp);
        //r(i) = greensMat(i, j, subn, rstore, rstore2, rr0, fac, dist, coord, snode, x) + (1. - volDt * sumbgvt(i)) * subnvar; // Inc G0 col
        tempStore += (1. - volDt * sumbgtv(i)) * x(i); // Bottom G0 row
    }
    r(subn) = tempStore - volDt * (nt - volDt * accubgtt) * x(subn); // G0 diagonal (final entry)

    r -= b;
    p = r;
    lu = accu(r);

    kk = 1;
    do  {

        v.zeros();
        tempStore = 0.;
        subnvar = p(subn);
#pragma omp parallel for default(none) shared(subn,isp) private(coord,rr0,fac,dist,i,j) firstprivate(subnvar,vseg,snode,p) schedule(dynamic) reduction(+:tempStore,vstore,vstore2)
        for(i = 0; i < subn; i++)  {
            vstore = 0.;
            vstore2 = 0.;
            coord = snode.col(i);
            rr0 = r0(i);
            fac = gfac(i,isp);
            //v(i) = greensMat(i, j, subn, vstore, vstore2, rr0, fac, dist, coord, snode, p) + (1. - volDt * sumbgvt(i)) * subnvar;
            tempStore += (1. - volDt * sumbgtv(i)) * p(i);
        }
        v(subn) = tempStore - volDt * (nt - volDt * accubgtt) * p(subn);


        t1 = accu(v % rs);
        delta = -lu/t1;

        s = r + delta * v;

        t.zeros();
        tempStore = 0.;
        subnvar = s(subn);
#pragma omp parallel for default(none) shared(subn,isp) private(coord,rr0,fac,dist,i,j) firstprivate(subnvar,vseg,snode,s) schedule(dynamic) reduction(+:tempStore,tstore,tstore2)
        for(i = 0; i < subn; i++)  {
            tstore = 0.;
            tstore2 = 0.;
            coord = snode.col(i);
            rr0 = r0(i);
            fac = gfac(i,isp);
            //t(i) = greensMat(i, j, subn, tstore, tstore2, rr0, fac, dist, coord, snode, s) + (1. - volDt * sumbgvt(i)) * subnvar;
            tempStore += (1. - volDt * sumbgtv(i)) * s(i);
        }
        t(subn) = tempStore - volDt * (nt - volDt * accubgtt) * s(subn);

        t1 = accu(s % t);
        t2 = accu(pow(t,2));
        gamma1 = -t1/t2;

        r = s + gamma1 * t;
        er = delta * p + gamma1 * s;
        err = max(abs(er));
        x += er;

        lunew = accu(r % rs);
        beta = lunew*delta/(lu*gamma1);
        lu = lunew;

        p = r + beta * (p + gamma1 * v);

        printf("Loop %i: Error: %e\n",kk,err);
        kk += 1;
    }
    while(kk <= itmax && err > eps);

    //if(err > eps) outputf(ltg, slog, "*** Warning: linear solution using BICGSTB not converged, err = ", err, "");
    return err;
}*/