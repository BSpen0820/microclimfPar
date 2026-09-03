// microclimfCppfunctions.h
// Shared includes, constants, and forward declarations for helper functions
// defined in microclimfCpp.cpp. Included by both microclimfCpp.cpp and
// microclimfParallel.cpp so the parallel workers can resolve all symbols at link time.
#pragma once
#include <Rcpp.h>
#include "microclimfheaders.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <tuple>
#include <cfloat>
#include <numeric>
#include <queue>
using namespace Rcpp;

// Physical / mathematical constants
constexpr double pi    = 3.14159265358979323846;
constexpr double torad = 3.14159265358979323846 / 180.0;
constexpr double sb    = 5.67e-8;
constexpr double thetam = 0.365;
constexpr double ka    = 0.4;
constexpr double omdy  = (2.0 * 3.14159265358979323846) / (24.0 * 3600.0);

// ---- Forward declarations of all non-exported (and some exported) helpers ----
// Solar / geometry
double radem(double tc);
int    juldayCpp(int year, int month, int day);
double soltimeCpp(int jd, double lt, double lon);
solmodel solpositionCpp(double lat, double lon, int year, int month, int day, double lt);
double solarindexCpp(double slope, double aspect, double zend, double azid, bool shadowmask = false);
kstruct cankCpp(double zenr, double x, double si);

// Two-stream radiation
tsdifstruct twostreamdifCpp(double pait, double x, double lref, double ltra, double gref);
tsdirstruct twostreamdirCpp(double pait, double om, double a, double gma, double J,
                            double del, double h, double gref, double kd,
                            double u1, double S1, double D1, double D2);
radmodel RadswabsCpp(double pai, double x, double lref, double ltra, double clump,
                     double gref, double slope, double aspect, double lat, double lon,
                     std::vector<int> year, std::vector<int> month, std::vector<int> day,
                     std::vector<double> lt, std::vector<double> Rsw, std::vector<double> Rdif);

// Air properties
double phairCpp(double tc, double pk);
double cpairCpp(double tc);

// Exported scalar helpers (defined in microclimfCpp.cpp)
double zeroplanedisCpp(double h, double pai);
double roughlengthCpp(double h, double pai, double d, double psi_h);
double satvapCpp(double tc);
double dewpointCpp(double ea);
double mincondCpp(double leafabs, double gs, double tc, double leafd);
NumericMatrix soildCppm(NumericMatrix twi, double Smin, double Smax, double tfact);
std::vector<double> hourtodayCpp(std::vector<double> hourly, std::string stat, bool rephour = true);
std::vector<double> manCpp(std::vector<double> x, int n);

// Diabatic / stability
double dpsimCpp(double ze);
double dpsihCpp(double ze);
double dphihCpp(double ze);
double gfreeCpp(double leafd, double H);
double gturbCpp(double uf, double d, double zm, double zref, double ph, double psi_h, double gmin);

// Stomatal / conductance
double psiwfromthetaCpp(double theta, double Smax, double psi_e, double b);
stompstruct stomparamsCpp(double hgt, double lat, double x);
double stomcondCpp(double Rswabs, double theta, double gsmax, double Smax,
                   double psi_e, double b, stompstruct stomp);
double canopycondCpp(double Rsw, double Rdif, double k, double om, double theta,
                     double gsmax, double PAI, double Smax, double psi_e,
                     double b, stompstruct stomp);

// Penman-Monteith
double PenmanMonteithCpp(double Rabs, double gHa, double gV, double tc, double te,
                         double pk, double ea, double em, double G, double erh);
penmonstruct PenmanMonteith2Cpp(double Rabs, double gHa, double gV, double tc,
                                double mxtc, double pk, double ea, double es,
                                double G, double surfwet, double tdew);

// Rolling means
std::vector<double> maCpp(std::vector<double> x, int n);
std::vector<double> mayCpp(std::vector<double> x);

// Soil
soilstruct soilpfun(double Vm, double Vq, double Mc, double rho);
double soildCpp(double soilm, double Smin, double Smax, double tadd);
soilkstruct soilcondCpp(double rho, double soilm, soilstruct sp);
soilmodelG0 soiltempG0(double tc, double es, double ea, double pk,
                       double radGsw, double radGlw, double tdew, double gHa,
                       double soilm, double mxtc, soilpstruct soilparams);
soilmodel soiltemp_hrCpp(double tc, double es, double ea, double pk,
                         double radabs, double surfwet, double tdew, double gHa,
                         double soilm, double mxtc, double Gp, double dtr,
                         double dtrp, double muGp, double kp, double Rdmx,
                         soilstruct sp, soilpstruct soilparams);

// Radiation (combined model)
tirstruct twostreamdif(double pai, double paia, double x, double lref, double ltra,
                       double clump, double gref);
radmodel2 twostreamCpp(double pai, double clump, double gref, double svfa, double si,
                       double tc, double Rsw, double Rdif, double Rlw,
                       solmodel solp, kstruct kp, tsdirstruct tspdir, tirstruct tir);

// Wind
tiwstruct windtiCpp(double h, double pai);
windmodel windCpp(double reqhgt, double zref, double h, double pai,
                  double uref, double umu, double ws, tiwstruct tiw);

// Above-ground temperature
abovecanstruct TVabove(double reqhgt, double zref, double h, double d, double zm,
                       double T0, double tc, double ea, double surfwet = 1.0);
leaftempstruct leaftemp(double Tcan, double Tg, double tc, double mxtc, double pk,
                        double ea, double es, double uz, double tdew, double surfwet,
                        double radLsw, double Rddown, double Rbdown, double Rlw,
                        double pai, double paia, double leafd, double gsmax,
                        double PARabs, double theta, double Smax, double psi_e,
                        double soilb, stompstruct stomp);
double rhcanopy(double uf, double h, double d, double z, double H, double Tk);
double TVbelow(double zref, double z, double d, double h, double pai, double uf,
               double leafden, double Flux, double Fluxz, double SH, double SG,
               double mxnear, double Hflux, double Tk);
abovemodel TVaboveground(double reqhgt, double zref, double tc, double pk, double ea,
                         double es, double tdew, double Rsw, double Rdif, double Rlw,
                         double soilm, double hgt, double pai, double paia, double vegx,
                         double leafd, double leafden, double Smin, double Smax,
                         double psi_e, double soilb, double gsmax, double mxtc,
                         stompstruct stomp, tirstruct tir, radmodel2 rvars,
                         tiwstruct tiw, windmodel wvars, soilmodel Gvars);
std::vector<double> Tbelowgroundv(double reqhgt, std::vector<double> Tg,
                                  std::vector<double> Tgp, std::vector<double> Tbp,
                                  double meanD, double mat, int hiy, bool complete);

// Snow helpers
double canopysnowintCpp(double hgt, double pai, double uf, double prec,
                        double tc, double Li, double Sh = 6.2);
std::vector<double> snowdenp(std::string snowenv);
NumericVector snowalbCpp(NumericVector prec);
snowrad radoneB(obspoint obstime, climpoint clim, vegpoint vegp,
                snowpoint snow, otherpoint other);
snowmodpoint snowoneB(obspoint obstime, climpoint clim, vegpoint vegp,
                      snowpoint snow, otherpoint other,
                      std::vector<double> sdp, double umu = 1.0);
NumericVector GFluxCppsnow(NumericVector snowt, NumericVector snowden);
NumericMatrix bioclimfill(int rows, int cols);
// Snow microclimate helpers (used by gridmicrosnow workers)
NumericVector snowdayan(NumericVector stempg);
NumericMatrix meanDsnow(NumericVector snowden);
NumericVector snowdaymov(NumericVector Tgref, NumericVector snowdays, NumericMatrix meanD, double reqhgt);
double belowpointsnow(double reqhgt, double meanD, double snowtempg, double Tzd, double Tza, double hiy);
snowmicro snowabovepoint(double reqhgt, double zref, double tc, double relhum, double pk, double u2,
                         double Rsw, double Rdif, double Rlw,
                         double hgt, double pai, double paia, double leafd, double clump, double ltra, double leafden,
                         solmodel solp, double si, double svfa, int shadowmask, double ws, double umu, double mxtc, snowpoint2 snowp);

// Ground heat flux model (used by BigLeafCpp only, declared for completeness)
Gmodel GFluxCpp(std::vector<double> Tg, std::vector<double> soilm,
                double rho, double Vm, double Vq, double Mc,
                std::vector<double> Gmax, std::vector<double> Gmin,
                int iter, bool yearG = true);
