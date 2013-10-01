//
//  Functions for evaluating the pareto density, distribution, and 
//  quantile function
//
//  Filename   : paretoC.c
//  Assignment : hw10
//
//  This version created by Philip Erickson on 4/12/13.
//

#include <stdio.h>
#include "math.h"
#include "stdlib.h"
#include <stdbool.h>
#include <ctype.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Random.h"

static double errorCheck(double x, double a, double b, double val) {
    // Error handling
    if (a <= 0 || b <= 0) {
        return NA_REAL;
    }
    return val;
}

static double *vectorize(double *v, int nv, int n) {
    // Vectorizor
    int ind, i;
    double *vnew = NULL;
    vnew = malloc(sizeof(double)*n);
    if (nv == n)
        return v;
    else {
        for (i = 0; i < n; i++) {
            ind = i % nv;
            vnew[i] = v[ind];
        }
        return vnew;
    }
}

/***********************/
/* Evaluation routines */
/***********************/
static double dparetoEval(double x, double a, double b, int l) {
    // Density
    double paretoi;
    if (x <= a)
        paretoi = log(0);
    else
        paretoi = log(b) + b * log(a) - (b + 1) * log(x);
    if (l == false)
        paretoi = exp(paretoi);
    return paretoi;
}

static double pparetoEvalL(double q, double a, double b, int l) {
    // Distribution (lower tail)
    double paretoi;
    if (q <= a)
        paretoi = 0;
    else
        paretoi = 1 - pow((a / q), b);
    if (l == true)
        paretoi = log(paretoi);
    return paretoi;
}

static double pparetoEvalU(double q, double a, double b, int l) {
    // Distribution (upper tail)
    double paretoi;
    if (q <= a)
        paretoi = 1;
    else
        paretoi = pow((a / q), b);
    if (l = true)
        paretoi = log(paretoi);
    return paretoi;
}

static double FReturn(double p, int logp) {
    // Control for logged inputs in quantile function
    double f;
    if (logp == true)
        f = exp(p);
    else
        f = p;
    return f;
}

static double qparetoEvalL(double p, double a, double b, int logp) {
    // Quantile (lower tail)
    double paretoi, f;
    f = FReturn(p, logp);
    if (logp == false && (p < 0 || p > 1)) {
        return NA_REAL;
    }
    else {
        paretoi = log(a) - log(1 - f) / b;
        paretoi = exp(paretoi);
    }
    return paretoi;
}

static double qparetoEvalU(double p, double a, double b, int logp) {
    // Quantile (upper tail)
    double paretoi, f;
    f = FReturn(p, logp);
    if (logp == false && (p < 0 || p > 1)) {
        return NA_REAL;
    }
    else {
        paretoi = log(a) - log(f) / b;
        paretoi = exp(paretoi);
    }
    return paretoi;
}

/*****************************/
/* Primary Functions (no MP) */
/*****************************/
static void paretodens(double *xi, int *nx, double *alphai, int *na, 
                       double *betai, int *nb, int *logp, int *nMax,
                       double *val) {
    // Density
    int n, i;
    int naflag = FALSE;
    n  = nMax[0];
    xi = vectorize(xi, nx[0], n);
    alphai = vectorize(alphai, na[0], n);
    betai = vectorize(betai, nb[0], n);
    for (i = 0; i < n; i++) {
        val[i] = dparetoEval(xi[i], alphai[i], betai[i], logp[0]);
        val[i] = errorCheck(xi[i], alphai[i], betai[i], val[i]);
        if (ISNAN(val[i]))
            naflag = TRUE;           
    }
    if (naflag) warning("NAs generated");
}

static void paretodist(double *xi, int *nx, double *alphai, int *na, 
                       double *betai, int *nb, int *ltail, int *logp, 
                       int *nMax, double *val) {
    // CDF
    int n, i;
    int naflag = FALSE;
    n  = nMax[0];
    xi = vectorize(xi, nx[0], n);
    alphai = vectorize(alphai, na[0], n);
    betai = vectorize(betai, nb[0], n);
    for (i = 0; i < n; i++) {
        if (ltail[0] == true){
            val[i] = pparetoEvalL(xi[i], alphai[i], betai[i], logp[0]);
        }
        else {
            val[i] = pparetoEvalU(xi[i], alphai[i], betai[i], logp[0]);
        }
        val[i] = errorCheck(xi[i], alphai[i], betai[i], val[i]);
        if (ISNAN(val[i]))
            naflag = TRUE; 
    }
    if (naflag) warning("NAs generated");
}

static void paretoquant(double *xi, int *nx, double *alphai, int *na, 
                        double *betai, int *nb, int *ltail, int *logp, 
                        int *nMax, double *val) {
    // Quantile function
    int n, i;
    int naflag = FALSE;
    n  = nMax[0];
    xi = vectorize(xi, nx[0], n);
    alphai = vectorize(alphai, na[0], n);
    betai = vectorize(betai, nb[0], n);
    for (i = 0; i < n; i++) {  
        if (ltail[0] == true){
            val[i] = qparetoEvalL(xi[i], alphai[i], betai[i], logp[0]);
        }
        else {
            val[i] = qparetoEvalU(xi[i], alphai[i], betai[i], logp[0]);
        }
        val[i] = errorCheck(xi[i], alphai[i], betai[i], val[i]);
        if (ISNAN(val[i]))
            naflag = TRUE;
    }
    if (naflag) warning("NAs generated");
}

/**************************/
/* Primary Functions (MP) */
/**************************/
static void paretodensP(double *xi, int *nx, double *alphai, int *na, 
                       double *betai, int *nb, int *logp, int *nMax,
                       double *val, int *P) {
    // Density
    int n, i;
    int naflag = FALSE;
    n  = nMax[0];
    xi = vectorize(xi, nx[0], n);
    alphai = vectorize(alphai, na[0], n);
    betai = vectorize(betai, nb[0], n);
#pragma omp parallel for num_threads(P[0]) default(none) \
    firstprivate(n, xi, alphai, betai, logp, val) private(i) \
    reduction(||:naflag)
    for (i = 0; i < n; i++) {
        val[i] = dparetoEval(xi[i], alphai[i], betai[i], logp[0]);
        val[i] = errorCheck(xi[i], alphai[i], betai[i], val[i]);
        if (ISNAN(val[i]))
            naflag = TRUE;           
    }
    if (naflag) warning("NAs generated");
}

static void paretodistP(double *xi, int *nx, double *alphai, int *na, 
                       double *betai, int *nb, int *ltail, int *logp, 
                       int *nMax, double *val, int *P) {
    // CDF
    int n, i;
    int naflag = FALSE;
    n  = nMax[0];
    xi = vectorize(xi, nx[0], n);
    alphai = vectorize(alphai, na[0], n);
    betai = vectorize(betai, nb[0], n);
#pragma omp parallel for num_threads(P[0]) default(none) \
    firstprivate(n, xi, alphai, betai, logp, val) private(i) \
    reduction(||:naflag)
    for (i = 0; i < n; i++) {
        if (ltail[0] == true){
            val[i] = pparetoEvalL(xi[i], alphai[i], betai[i], logp[0]);
        }
        else {
            val[i] = pparetoEvalU(xi[i], alphai[i], betai[i], logp[0]);
        }
        val[i] = errorCheck(xi[i], alphai[i], betai[i], val[i]);
        if (ISNAN(val[i]))
            naflag = TRUE; 
    }
    if (naflag) warning("NAs generated");
}

static void paretoquantP(double *xi, int *nx, double *alphai, int *na, 
                        double *betai, int *nb, int *ltail, int *logp, 
                        int *nMax, double *val, int *P) {
    // Quantile function
    int n, i;
    int naflag = FALSE;
    n  = nMax[0];
    xi = vectorize(xi, nx[0], n);
    alphai = vectorize(alphai, na[0], n);
    betai = vectorize(betai, nb[0], n);
#pragma omp parallel for num_threads(P[0]) default(none) \
    firstprivate(n, xi, alphai, betai, logp, val) private(i) \
    reduction(||:naflag)
    for (i = 0; i < n; i++) {  
        if (ltail[0] == true){
            val[i] = qparetoEvalL(xi[i], alphai[i], betai[i], logp[0]);
        }
        else {
            val[i] = qparetoEvalU(xi[i], alphai[i], betai[i], logp[0]);
        }
        val[i] = errorCheck(xi[i], alphai[i], betai[i], val[i]);
        if (ISNAN(val[i]))
            naflag = TRUE;
    }
    if (naflag) warning("NAs generated");
}

static void paretorand(int *ni, double *alphai, double *betai, double *val) {
    // Pseudorandom number generator
    int n = ni[0];
    int naflag = FALSE;
    int i;
    GetRNGstate(); 
    for (i = 0; i < n; i++) {
        val[i] = unif_rand();
        val[i] = qparetoEvalL(val[i], alphai[0], betai[0], FALSE);
        val[i] = errorCheck(1, alphai[0], betai[0], val[i]);
        if (ISNAN(val[i]))
            naflag = TRUE;
    }
    if (naflag) warning("NAs generated");
    PutRNGstate(); 
}

/* Data structure for registering routines. */
static R_CMethodDef DotCEntries[] = {
    {"paretodens", (DL_FUNC) paretodens, 9},
    {"paretodist", (DL_FUNC) paretodist, 10},
    {"paretoquant", (DL_FUNC) paretoquant, 10},
    {"paretodensP", (DL_FUNC) paretodensP, 10},
    {"paretodistP", (DL_FUNC) paretodistP, 11},
    {"paretoquantP", (DL_FUNC) paretoquantP, 11},
    {"paretorand", (DL_FUNC) paretorand, 4},
    {NULL}
};

/* This is called by the dynamic loader to register the routine. */
void R_init_pareto(DllInfo *info)
{
    R_registerRoutines(info, DotCEntries, NULL, NULL, NULL);
}





