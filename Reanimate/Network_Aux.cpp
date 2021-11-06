#include "Network.hpp"

using namespace std;
using namespace reanimate;

/*double pearson_coeff(vec vec1, vec vec2, int num)    {
    
    double r = (num*accu(vec1%vec2) - accu(vec1)*accu(vec2))/(sqrt(num*accu(pow(vec1,2)) - pow(accu(vec1),2))*sqrt(num*accu(pow(vec2,2)) - pow(accu(vec2),2)));
    
    return r;
}*/

int Network::detect_col(FILE *ifp)    {
    
    int max=200;
    char bb[200];
    string s;
    fpos_t position;
    fgetpos(ifp,&position);
    fgets(bb,max,ifp);
    s = bb;
    istringstream is( s );
    double n;
    int cnt = 0;
    while( is >> n ) {
        cnt += 1;
    }
    fsetpos(ifp,&position);
    
    return cnt;
    
}

// Initiate log file
void Network::initLog()    {

    FILE *ofp;
    ofp = fopen((buildPath + rLog).c_str(),"w");
    fclose(ofp);
    
}

// Output and log text
void Network::printText(const string &text, const int type, const int newline)   {

    // Type: (1) normal text, (2) loading description, (3) module, (4) error, (5) warning
    if (!silence)   {

        FILE *ofp;
        ofp = fopen((buildPath + rLog).c_str(),"a");
        if (type == 2)  {
            if (newline == 1)   {
                printf("\n%s ...\n", text.c_str());
                fprintf(ofp, "\n%s ... \n", text.c_str());
            }
            else {
                printf("%s ...\n", text.c_str());
                fprintf(ofp, "%s ... \n", text.c_str());
            }
        }
        else if (type == 3) {
            printf("\n\n");
            fprintf(ofp, "\n\n");
            for (int i = 0; i < (int) text.length(); i++)  {
                printf("-");
                fprintf(ofp, "-");
            }
            printf("\n%s\n", text.c_str());
            fprintf(ofp, "\n%s\n", text.c_str());
            for (int i = 0; i < (int) text.length(); i++) {
                printf("-");
                fprintf(ofp, "-");
            }
            printf("\n");
            fprintf(ofp, "\n");
        }
        else if (type == 4) {
            printf("*** ERROR: %s ***\n",text.c_str());
            fprintf(ofp, "*** ERROR: %s ***\n",text.c_str());
        }
        else if (type == 5) {
            printf("*** WARNING: %s ***\n",text.c_str());
            fprintf(ofp, "*** WARNING: %s ***\n",text.c_str());
        }
        else if (type == 6) {
            printf("\n\n");
            fprintf(ofp, "\n\n");
            for (int i = 0; i < (int) text.length(); i++)  {
                printf("-");
                fprintf(ofp, "-");
            }
            printf(" %s ", text.c_str());
            fprintf(ofp, " %s ", text.c_str());
            for (int i = 0; i < (int) text.length(); i++) {
                printf("-");
                fprintf(ofp, "-");
            }

            printf("\n");
            fprintf(ofp, "\n");

        }
        else {
            if (newline == 1)   {
                printf("\n%s\n", text.c_str());
                fprintf(ofp, "\n%s\n", text.c_str());
            }
            else if (newline == -1) {
                printf("%s", text.c_str());
                fprintf(ofp, "%s", text.c_str());
            }
            else    {
                printf("%s\n", text.c_str());
                fprintf(ofp, "%s\n", text.c_str());
            }
        }

        fclose(ofp);

    }
    
}

// Output and log "string : float dimensions"
void Network::printNum(const string &text, const double &num, const string unit)   {

    if (!silence)   {

        FILE *ofp;
        ofp = fopen((buildPath + rLog).c_str(),"a");

        if (round(num) == num)  {
            printf("%s %i %s\n",text.c_str(),(int) num,unit.c_str());
            fprintf(ofp,"%s %i %s\n",text.c_str(),(int) num,unit.c_str());
        }
        else if (num < 1e-3)    {
            printf("%s %.2e %s\n",text.c_str(),num,unit.c_str());
            fprintf(ofp,"%s %.4e %s\n",text.c_str(),num,unit.c_str());
        }
        else {
            printf("%s %.3f %s\n",text.c_str(),num,unit.c_str());
            fprintf(ofp,"%s %.5f %s\n",text.c_str(),num,unit.c_str());
        }

        fclose(ofp);

    }

}

// Output and log "string : float dimensions"
void Network::printStat(const string &text, const vec &n, const string &unit)   {

    if (!silence)   {

        FILE *ofp;
        ofp = fopen((buildPath + rLog).c_str(),"a");

        if (mean(n) < 1e-3)    {
            printf("%s %.2e ± %.4e %s\n",text.c_str(),mean(n),stddev(n),unit.c_str());
            fprintf(ofp,"%s %.4e ± %.4e %s\n",text.c_str(),mean(n),stddev(n),unit.c_str());
        }
        else {
            printf("%s %.3f ± %.3f %s\n",text.c_str(),mean(n),stddev(n),unit.c_str());
            fprintf(ofp,"%s %.5f ± %.5f %s\n",text.c_str(),mean(n),stddev(n),unit.c_str());
        }

        fclose(ofp);

    }
    
}

void Network::findBoundingBox() {

    double xmin = cnode.row(0).min();
    double ymin = cnode.row(1).min();
    double zmin = cnode.row(2).min();

    if (xmin >= 0.) {cnode.row(0) -= xmin;}
    else {cnode.row(0) += abs(xmin);}

    if (ymin >= 0.) {cnode.row(1) -= ymin;}
    else {cnode.row(1) += abs(ymin);}

    if (zmin >= 0.) {cnode.row(2) -= zmin;}
    else {cnode.row(2) += abs(zmin);}

    alx = cnode.row(0).max();
    aly = cnode.row(1).max();
    alz = cnode.row(2).max();
    if (alz == 0.)  {alz = diam.max();}

}


/*
int max_double(const double &a, const double &b) {
    
    if (a > b)  {
        return 1;
    }
    else if (a < b) {
        return 2;
    }
    else {
        return 1;
    }
    
}

double max_double_val(const double &a, const double &b) {
    
    if (a > b)  {
        return a;
    }
    else if (a < b) {
        return b;
    }
    else {
        return 0.;
    }
    
}

double min_double_val(const double &a, const double &b) {
    
    if (a > b)  {
        return b;
    }
    else if (a < b) {
        return a;
    }
    else {
        return 0.;
    }
    
}

int sgn(const double &a)   {
    
    int b = 0;
    if (a < 0)  {
        b = -1;
    }
    else if (a > 0) {
        b = 1;
    }
    
    return b;
}


ivec sphere_packing(const double grow_decay, vec &r0, const mat &snode)   {
    
    int circ_count = 0;
    int ns = (int) snode.n_cols;
    ivec flag_sphere = zeros<ivec>(ns);
    vec old_r0 = zeros<vec>(ns);
    flag_sphere = zeros<ivec>(ns);
    
    int itmax = 1e3;
    int iter = 0.;
    
    int old_circ_count = -1;
    int cnt = 0;
    while (circ_count < ns && iter < itmax)  {
        
        double dist = 0.0;
        for (int is = 0; is < ns; is++) {
            if (flag_sphere(is) == 0)   {
                int cntr = 0;
                for (int js = 0; js < ns; js++) {
                    if (is != js && flag_sphere(js) == 0 && grow_decay < 1.)   {
                        dist = sqrt(pow(snode(0,is)-snode(0,js),2) + pow(snode(1,is)-snode(1,js),2) + pow(snode(2,is)-snode(2,js),2)) - r0(is) - r0(js);
                        if (dist < 0.0)  {
                            cntr += 1;
                            js = ns;
                        }
                    }
                    else if (is != js && grow_decay > 1.)   {
                        dist = sqrt(pow(snode(0,is)-snode(0,js),2) + pow(snode(1,is)-snode(1,js),2) + pow(snode(2,is)-snode(2,js),2)) - r0(is) - r0(js);
                        if (dist < 0.0)  {
                            cntr += 1;
                            js = ns;
                        }
                    }
                }
                if (cntr == 0 && grow_decay < 1.)  {
                    flag_sphere(is) = 1;
                }
                else if (cntr == 0 && grow_decay > 1.)  {
                    flag_sphere(is) = 0;
                }
                else if (cntr > 0 && grow_decay > 1.)   {
                    flag_sphere(is) = 1;
                    r0(is) *= (2 - grow_decay);
                }
                if (cnt > 10 && grow_decay < 1.)   {
                    //flag_sphere(is) = 2;
                }
                else if (cnt > 1e2 && grow_decay > 1.)  {
                    //flag_sphere(is) = 2;
                }
            }
        }

        
        old_r0 = r0;
        r0(find(flag_sphere == 0)) *= grow_decay;
        circ_count = (int) accu(flag_sphere);
        
        if (old_circ_count == circ_count)   {
            cnt += 1;
        }
        else {
            cnt = 0;
        }
        old_circ_count = circ_count;
        
        iter += 1;
    }
    
    
    return flag_sphere;
}


// Produces radially dependent vector by normalising vascular network
vec var_param(const int &ns, const double &inner, const double &outer, const mat &snode, const double &alx, const double &aly, const double &alz)   {
    
    vec y = zeros<vec>(ns);
    mat norm = snode;
    
    // Normalise coordinates
    norm.row(0) /=  max(norm.row(0));
    norm.row(1) /=  max(norm.row(1));
    norm.row(2) /=  max(norm.row(2));
    norm.row(0) -= 0.5;
    norm.row(1) -= 0.5;
    norm.row(2) -= 0.5;
    
    double diff = abs(outer - inner);
    for (int i = 0; i < ns; i++)    {
        double r = sqrt(pow(norm(0,i),2) + pow(norm(1,i),2) + pow(norm(2,i),2));
        r *= 2;
        if (r > 1.) {
            r = 1;
        }
        if (outer > inner)  {
            y(i) = inner + r * diff;
        }
        else {
            y(i) = inner - r * diff;
        }
        
    }
    
    return y;
    
}


void time_check(FILE *ift, const string &local, double &run_start, const string &text)   {
    
    double time_span = (omp_get_wtime()-run_start);
    if (time_span > 3600)    {
        outputf(ift,local,text, time_span/3600,"hours");
    }
    else if (time_span > 60)    {
        outputf(ift,local,text, time_span/60,"minutes");
    }
    else    {
        outputf(ift,local,text, time_span,"seconds");
    }
    
}*/


