//
//  misc_func.h
//  Vascular-Flow
//
//  Created by Paul Sweeney on 04/04/2017.
//  Copyright © 2017 Paul Sweeney - University College London. All rights reserved.
//

#ifndef misc_func_h
#define misc_func_h

#include <stdio.h>
#include <armadillo>

using namespace std;
using namespace arma;

extern int detect_col(FILE *ifp);
extern int max_double(const double &a, const double &b);
extern int sgn(const double &a);

extern double pearson_coeff(vec vec1, vec vec2, int num);
extern double max_double_val(const double &a, const double &b);
extern double min_double_val(const double &a, const double &b);

extern ivec sphere_packing(const double grow_decay, vec &r0, const mat &snode);

extern vec var_param(const int &ns, const double &inner, const double &outer, const mat &snode, const double &alx, const double &aly, const double &alz);

extern void init_log(FILE *ofp, const string &local);
extern void outputt(FILE *ofp, const string &local, const string &state);
extern void outputf(FILE *ofp, const string &local, const string &text, const double &num, const string &unit);
extern void output2f(FILE *ofp, const string &local, const string &text, const double &num1, const double &num2, const string &unit);
extern void time_check(FILE *ift, const string &local, double &run_start, const string &text);
extern void ascii_output(const string &filename);
extern void vessel_order(ivec master,vec v_lseg);


#endif /* misc_func_h */
