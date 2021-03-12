//
// Created by Paul Sweeney on 14/02/2021.
//

#include "Network.hpp"

using namespace reanimate;

// Euclidean distance between two points
double Network::eucDistance(vec &x, vec &y) {

    double val{};
    for (int i = 0; i < (int) x.n_elem; i++)    {val += pow(x(i)-y(i),2);}

    return sqrt(val);

}

// Lognormal random generator
vec Network::lognormal(double &mean, double &SD, int &nseg)  {

    double v = pow(SD,2);
    double mu = log(mean/sqrt(1 + v/pow(mean,2)));
    double sig = sqrt(log(1 + v/pow(mean,2)));

    default_random_engine engine(std::random_device{}());
    lognormal_distribution<double> distribution(mu,sig);

    vec output = zeros<vec>(nseg);
    for (int i = 0; i < nseg; i++)    {
        output(i) = distribution(engine);
    }

    return output;
}

mat Network::eulerRotation(const vec &c, const vec &d)    {

    vec a = c / norm(c);
    vec b = d / norm(d);

    vec v = cross(a,b);
    mat ssc = zeros<mat>(3,3);
    ssc(0,1) = -v(2);
    ssc(0,2) = v(1);
    ssc(1,0) = v(2);
    ssc(1,2) = -v(0);
    ssc(2,0) = -v(1);
    ssc(2,1) = v(0);

    mat R = eye<mat>(3,3) + ssc + (ssc * ssc)*(1-dot(a,b))/pow(norm(v),2);
    return R;
}

// Modified spherical bessel function of a first kind - computes either zeroth or first order
double Network::SPHI(int &n, double &x) {

    double SI{};
    if (n != 0 && n != 1)   {printText("Modified spherical Bessel function (1st kind) order out of bounds", 4);}

    if (fabs(x) < 1e-100) {
        if (n == 0) {SI = 1.0;}
        else if (n == 1)    {SI = 0.0;}
        return SI;
    }
    if (n == 0) {SI = sinh(x) / x;}
    else if (n == 1)    {SI = (x * cosh(x) - sinh(x)) / pow(x, 2);}

    return SI;

}


// Modified spherical bessel function of a second kind - computers either zeroth or first order
double Network::SPHK(int &n, double &x) {

    double SK{};
    if (n != 0 && n != 1)   {printText("Modified spherical Bessel function (2nd kind) order out of bounds", 4);}

    if (x < 1e-60) {
        SK = 1.0e+300;
        return  SK;
    }

    if (n == 0) {SK = 0.5 * M_PI / x * exp(-x);}
    else if (n == 1)    {
        double SK0 = 0.5 * M_PI / x * exp(-x);
        SK = SK0*(1.0+ 1.0 / x);
    }

    return SK;
}
