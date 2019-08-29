#ifndef FUNCTIONS_H
#define FUNCTIONS_H

/*
 *
 * functions.h: Special functions
 * BRIAN Software Package Version 3.0
 *
 * $Id: functions.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.30 (21/04/13): released version 2.4
 * 0.40 (19/12/13): documented
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Defines interfaces to special math functions.
*/

extern double qGamma(const double a, const double x);
extern double pGamma(const double a, const double x);
extern double iGamma(const double a, const double x);
extern double MarcumQ(const double nu, const double a, const double b, const double tol = 1e-16);
extern double digamma(const double x);
extern double iBeta(const double a, const double b, const double x);
extern double lBeta(const double x, const double y);
extern double Beta(const double x, const double y);
extern double rBeta(const double a, const double b, const double x);
extern double besseli_sc(const double nu, const double x);
extern double inverff(const double y);
extern double lmvGamma(const unsigned int p, const unsigned int n);
extern double arch(const double alpha);
extern double Hyper1F1(const double x, const double a, const double b, const double eps = 1e-10);
extern double logHyper1F1(const double x, const double a, const double b, const double eps = 1e-10);
extern double stirlerr(const double n);
extern double bd0(const double x, const double np);
extern double sBesselI0(const double x);
extern double BesselRatio(double v);
extern double sBesselI1(const double x);
extern double BesselI1(const double x);
extern double LaguerreHalf(const double x);
extern double HyperRatio(const double k0, const double k1, const double ki);
extern double LegendrePlm(const int l, const int m, const double x);
extern std::complex<double> YlmC(const int l, const int m, const double theta, const double phi);
extern double YlmR(const int l, const int m, const double theta, const double phi);
extern double ShapiroTest(const dvecD& x);
extern double splineK(const double a);
#endif

