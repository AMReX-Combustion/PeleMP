#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

/* Elements
0  C
1  H
2  O
3  N
*/

// Species
#define NC5H12_ID 0
#define NC6H14_ID 1
#define NC8H18_ID 2
#define IC8H18_ID 3
#define NC10H22_ID 4
#define NC7H16_ID 5
#define NC12H26_ID 6
#define O2_ID 7
#define N2_ID 8
#define NC16H34_ID 9
#define H2O_ID 10
#define CO2_ID 11

#define NUM_ELEMENTS 4
#define NUM_SPECIES 12
#define NUM_IONS 0
#define NUM_REACTIONS 0

#define NUM_FIT 4

//  ALWAYS on CPU stuff -- can have different def depending on if we are CPU or GPU based. Defined in mechanism.cpp 
void atomicWeight(amrex::Real *  awt);
//  MISC 
void CKAWT(amrex::Real *  awt);
void CKNCF(int * ncf);
void CKSYME_STR(amrex::Vector<std::string>& ename);
void CKSYMS_STR(amrex::Vector<std::string>& kname);
void GET_RMAP(int * _rmap);
void CKINU(const int i, int &nspec, int * ki, int * nu);
void CKKFKR(const amrex::Real P, const amrex::Real T,const amrex::Real * x, amrex::Real *  q_f, amrex::Real *  q_r);
void progressRateFR(amrex::Real *  q_f, amrex::Real *  q_r,amrex::Real *  sc, amrex::Real T);
//  SPARSE INFORMATION 
void SPARSITY_INFO(int * nJdata, const int * consP, int NCELLS);
void SPARSITY_INFO_SYST(int * nJdata, const int * consP, int NCELLS);
void SPARSITY_INFO_SYST_SIMPLIFIED(int * nJdata, const int * consP);
void SPARSITY_PREPROC_CSC(int * rowVals, int * colPtrs, const int * consP, int NCELLS);
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, const int * consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtrs, const int * consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, const int * consP);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, const int * consP, int base);


// A few mechanism parameters
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKINDX(int& mm, int& kk, int& ii, int& nfit)
{
mm = 4;
kk = 12;
ii = 0;
nfit = -1; // Why do you need this anyway ? 
}

//  inverse molecular weights 
AMREX_GPU_CONSTANT const amrex::Real global_imw[12]={
0.0138598217626921,// NC5H12
0.0116038896238019,// NC6H14
0.0087541144337839,// NC8H18
0.0087541144337839,// IC8H18
0.0070280983371519,// NC10H22
0.0099795419390250,// NC7H16
0.0058706117177410,// NC12H26
0.0312519532470779,// O2
0.0356964374955379,// N2
0.0044160248710521,// NC16H34
0.0555092978073827,// H2O
0.0227226249176305,// CO2
};
const amrex::Real h_global_imw[12]={
0.0138598217626921,// NC5H12
0.0116038896238019,// NC6H14
0.0087541144337839,// NC8H18
0.0087541144337839,// IC8H18
0.0070280983371519,// NC10H22
0.0099795419390250,// NC7H16
0.0058706117177410,// NC12H26
0.0312519532470779,// O2
0.0356964374955379,// N2
0.0044160248710521,// NC16H34
0.0555092978073827,// H2O
0.0227226249176305,// CO2
};

//  molecular weights 
AMREX_GPU_CONSTANT const amrex::Real global_mw[12]={
72.151000,// NC5H12
86.178000,// NC6H14
114.232000,// NC8H18
114.232000,// IC8H18
142.286000,// NC10H22
100.205000,// NC7H16
170.340000,// NC12H26
31.998000,// O2
28.014000,// N2
226.448000,// NC16H34
18.015000,// H2O
44.009000,// CO2
};
const amrex::Real h_global_mw[12]={
72.151000,// NC5H12
86.178000,// NC6H14
114.232000,// NC8H18
114.232000,// IC8H18
142.286000,// NC10H22
100.205000,// NC7H16
170.340000,// NC12H26
31.998000,// O2
28.014000,// N2
226.448000,// NC16H34
18.015000,// H2O
44.009000,// CO2
};

//  inverse molecular weights 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_imw(amrex::Real *imw_new){
imw_new[0] = 0.0138598217626921;// NC5H12
imw_new[1] = 0.0116038896238019;// NC6H14
imw_new[2] = 0.0087541144337839;// NC8H18
imw_new[3] = 0.0087541144337839;// IC8H18
imw_new[4] = 0.0070280983371519;// NC10H22
imw_new[5] = 0.0099795419390250;// NC7H16
imw_new[6] = 0.0058706117177410;// NC12H26
imw_new[7] = 0.0312519532470779;// O2
imw_new[8] = 0.0356964374955379;// N2
imw_new[9] = 0.0044160248710521;// NC16H34
imw_new[10] = 0.0555092978073827;// H2O
imw_new[11] = 0.0227226249176305;// CO2
}

//  inverse molecular weight 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real imw(const int n){
#if AMREX_DEVICE_COMPILE
return global_imw[n];
#else
return h_global_imw[n];
#endif
}
//  molecular weights 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_mw(amrex::Real *mw_new){
mw_new[0] = 72.151000;// NC5H12
mw_new[1] = 86.178000;// NC6H14
mw_new[2] = 114.232000;// NC8H18
mw_new[3] = 114.232000;// IC8H18
mw_new[4] = 142.286000;// NC10H22
mw_new[5] = 100.205000;// NC7H16
mw_new[6] = 170.340000;// NC12H26
mw_new[7] = 31.998000;// O2
mw_new[8] = 28.014000;// N2
mw_new[9] = 226.448000;// NC16H34
mw_new[10] = 18.015000;// H2O
mw_new[11] = 44.009000;// CO2
}

//  molecular weight 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real mw(const int n){
#if AMREX_DEVICE_COMPILE
return global_mw[n];
#else
return h_global_mw[n];
#endif
}

//  Returns R, Rc, Patm
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRP(amrex::Real& ru, amrex::Real& ruc, amrex::Real& pa)
{
 ru  = 8.31446261815324e+07; 
 ruc = 1.98721558317399615845; 
 pa  = 1.01325e+06; 
}

// compute Cv/R at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cv_R(amrex::Real * species, const amrex::Real *  tc)
{

// temperature
const amrex::Real T = tc[1];

// species with midpoint at T=1390 kelvin
if (T < 1390) {
// species 0: NC5H12
species[0] =
-1.73676612e+00
+6.07200973e-02 * tc[1]
-3.57592761e-05 * tc[2]
+1.04907042e-08 * tc[3]
-1.21487315e-12 * tc[4]
;
// species 1: NC6H14
species[1] =
-1.96960618e+00
+7.29085608e-02 * tc[1]
-4.38853919e-05 * tc[2]
+1.32312807e-08 * tc[3]
-1.58437423e-12 * tc[4]
;} else {
// species 0: NC5H12
species[0] =
+1.47257603e+01
+2.61086045e-02 * tc[1]
-8.90970996e-06 * tc[2]
+1.38102248e-09 * tc[3]
-8.00296536e-14 * tc[4]
;
// species 1: NC6H14
species[1] =
+1.79634117e+01
+3.04480204e-02 * tc[1]
-1.03794829e-05 * tc[2]
+1.60775457e-09 * tc[3]
-9.31269728e-14 * tc[4]
;}

// species with midpoint at T=1391 kelvin
if (T < 1391) {
// species 2: NC8H18
species[2] =
-2.54218406e+00
+9.78112063e-02 * tc[1]
-6.09318358e-05 * tc[2]
+1.92005591e-08 * tc[3]
-2.42996250e-12 * tc[4]
;
// species 4: NC10H22
species[4] =
-3.08416969e+00
+1.22535012e-01 * tc[1]
-7.76815739e-05 * tc[2]
+2.49834877e-08 * tc[3]
-3.23548038e-12 * tc[4]
;
// species 5: NC7H16
species[5] =
-2.26836187e+00
+8.54355820e-02 * tc[1]
-5.25346786e-05 * tc[2]
+1.62945721e-08 * tc[3]
-2.02394925e-12 * tc[4]
;
// species 6: NC12H26
species[6] =
-3.62181594e+00
+1.47237711e-01 * tc[1]
-9.43970271e-05 * tc[2]
+3.07441268e-08 * tc[3]
-4.03602230e-12 * tc[4]
;} else {
// species 2: NC8H18
species[2] =
+2.44710194e+01
+3.90887037e-02 * tc[1]
-1.33038777e-05 * tc[2]
+2.05867527e-09 * tc[3]
-1.19167174e-13 * tc[4]
;
// species 4: NC10H22
species[4] =
+3.09882239e+01
+4.77244922e-02 * tc[1]
-1.62276391e-05 * tc[2]
+2.50963259e-09 * tc[3]
-1.45215772e-13 * tc[4]
;
// species 5: NC7H16
species[5] =
+2.12148969e+01
+3.47675750e-02 * tc[1]
-1.18407129e-05 * tc[2]
+1.83298478e-09 * tc[3]
-1.06130266e-13 * tc[4]
;
// species 6: NC12H26
species[6] =
+3.75095037e+01
+5.63550048e-02 * tc[1]
-1.91493200e-05 * tc[2]
+2.96024862e-09 * tc[3]
-1.71244150e-13 * tc[4]
;}

// species with midpoint at T=1396 kelvin
if (T < 1396) {
// species 3: IC8H18
species[3] =
-5.20868893e+00
+1.11440581e-01 * tc[1]
-7.91346582e-05 * tc[2]
+2.92406242e-08 * tc[3]
-4.43743191e-12 * tc[4]
;} else {
// species 3: IC8H18
species[3] =
+2.61373590e+01
+3.79004890e-02 * tc[1]
-1.29437358e-05 * tc[2]
+2.00760372e-09 * tc[3]
-1.16400580e-13 * tc[4]
;}

// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 7: O2
species[7] =
+2.78245636e+00
-2.99673416e-03 * tc[1]
+9.84730201e-06 * tc[2]
-9.68129509e-09 * tc[3]
+3.24372837e-12 * tc[4]
;
// species 8: N2
species[8] =
+2.29867700e+00
+1.40824000e-03 * tc[1]
-3.96322200e-06 * tc[2]
+5.64151500e-09 * tc[3]
-2.44485500e-12 * tc[4]
;
// species 10: H2O
species[10] =
+2.38684249e+00
+3.47498246e-03 * tc[1]
-6.35469633e-06 * tc[2]
+6.96858127e-09 * tc[3]
-2.50658847e-12 * tc[4]
;
// species 11: CO2
species[11] =
+1.35677352e+00
+8.98459677e-03 * tc[1]
-7.12356269e-06 * tc[2]
+2.45919022e-09 * tc[3]
-1.43699548e-13 * tc[4]
;} else {
// species 7: O2
species[7] =
+2.28253784e+00
+1.48308754e-03 * tc[1]
-7.57966669e-07 * tc[2]
+2.09470555e-10 * tc[3]
-2.16717794e-14 * tc[4]
;
// species 8: N2
species[8] =
+1.92664000e+00
+1.48797700e-03 * tc[1]
-5.68476100e-07 * tc[2]
+1.00970400e-10 * tc[3]
-6.75335100e-15 * tc[4]
;
// species 10: H2O
species[10] =
+1.67214561e+00
+3.05629289e-03 * tc[1]
-8.73026011e-07 * tc[2]
+1.20099639e-10 * tc[3]
-6.39161787e-15 * tc[4]
;
// species 11: CO2
species[11] =
+2.85746029e+00
+4.41437026e-03 * tc[1]
-2.21481404e-06 * tc[2]
+5.23490188e-10 * tc[3]
-4.72084164e-14 * tc[4]
;}

// species with midpoint at T=1800 kelvin
if (T < 1800) {
// species 9: NC16H34
species[9] =
-4.64057977e+00
+1.96191732e-01 * tc[1]
-1.26758759e-04 * tc[2]
+4.13645475e-08 * tc[3]
-5.40570458e-12 * tc[4]
;} else {
// species 9: NC16H34
species[9] =
+4.88210237e+01
+7.73881688e-02 * tc[1]
-2.77557896e-05 * tc[2]
+4.69678103e-09 * tc[3]
-3.12959243e-13 * tc[4]
;}
}

// compute Cp/R at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cp_R(amrex::Real * species, const amrex::Real *  tc)
{

// temperature
const amrex::Real T = tc[1];

// species with midpoint at T=1390 kelvin
if (T < 1390) {
// species 0: NC5H12
species[0] =
-7.36766122e-01
+6.07200973e-02 * tc[1]
-3.57592761e-05 * tc[2]
+1.04907042e-08 * tc[3]
-1.21487315e-12 * tc[4]
;
// species 1: NC6H14
species[1] =
-9.69606184e-01
+7.29085608e-02 * tc[1]
-4.38853919e-05 * tc[2]
+1.32312807e-08 * tc[3]
-1.58437423e-12 * tc[4]
;} else {
// species 0: NC5H12
species[0] =
+1.57257603e+01
+2.61086045e-02 * tc[1]
-8.90970996e-06 * tc[2]
+1.38102248e-09 * tc[3]
-8.00296536e-14 * tc[4]
;
// species 1: NC6H14
species[1] =
+1.89634117e+01
+3.04480204e-02 * tc[1]
-1.03794829e-05 * tc[2]
+1.60775457e-09 * tc[3]
-9.31269728e-14 * tc[4]
;}

// species with midpoint at T=1391 kelvin
if (T < 1391) {
// species 2: NC8H18
species[2] =
-1.54218406e+00
+9.78112063e-02 * tc[1]
-6.09318358e-05 * tc[2]
+1.92005591e-08 * tc[3]
-2.42996250e-12 * tc[4]
;
// species 4: NC10H22
species[4] =
-2.08416969e+00
+1.22535012e-01 * tc[1]
-7.76815739e-05 * tc[2]
+2.49834877e-08 * tc[3]
-3.23548038e-12 * tc[4]
;
// species 5: NC7H16
species[5] =
-1.26836187e+00
+8.54355820e-02 * tc[1]
-5.25346786e-05 * tc[2]
+1.62945721e-08 * tc[3]
-2.02394925e-12 * tc[4]
;
// species 6: NC12H26
species[6] =
-2.62181594e+00
+1.47237711e-01 * tc[1]
-9.43970271e-05 * tc[2]
+3.07441268e-08 * tc[3]
-4.03602230e-12 * tc[4]
;} else {
// species 2: NC8H18
species[2] =
+2.54710194e+01
+3.90887037e-02 * tc[1]
-1.33038777e-05 * tc[2]
+2.05867527e-09 * tc[3]
-1.19167174e-13 * tc[4]
;
// species 4: NC10H22
species[4] =
+3.19882239e+01
+4.77244922e-02 * tc[1]
-1.62276391e-05 * tc[2]
+2.50963259e-09 * tc[3]
-1.45215772e-13 * tc[4]
;
// species 5: NC7H16
species[5] =
+2.22148969e+01
+3.47675750e-02 * tc[1]
-1.18407129e-05 * tc[2]
+1.83298478e-09 * tc[3]
-1.06130266e-13 * tc[4]
;
// species 6: NC12H26
species[6] =
+3.85095037e+01
+5.63550048e-02 * tc[1]
-1.91493200e-05 * tc[2]
+2.96024862e-09 * tc[3]
-1.71244150e-13 * tc[4]
;}

// species with midpoint at T=1396 kelvin
if (T < 1396) {
// species 3: IC8H18
species[3] =
-4.20868893e+00
+1.11440581e-01 * tc[1]
-7.91346582e-05 * tc[2]
+2.92406242e-08 * tc[3]
-4.43743191e-12 * tc[4]
;} else {
// species 3: IC8H18
species[3] =
+2.71373590e+01
+3.79004890e-02 * tc[1]
-1.29437358e-05 * tc[2]
+2.00760372e-09 * tc[3]
-1.16400580e-13 * tc[4]
;}

// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 7: O2
species[7] =
+3.78245636e+00
-2.99673416e-03 * tc[1]
+9.84730201e-06 * tc[2]
-9.68129509e-09 * tc[3]
+3.24372837e-12 * tc[4]
;
// species 8: N2
species[8] =
+3.29867700e+00
+1.40824000e-03 * tc[1]
-3.96322200e-06 * tc[2]
+5.64151500e-09 * tc[3]
-2.44485500e-12 * tc[4]
;
// species 10: H2O
species[10] =
+3.38684249e+00
+3.47498246e-03 * tc[1]
-6.35469633e-06 * tc[2]
+6.96858127e-09 * tc[3]
-2.50658847e-12 * tc[4]
;
// species 11: CO2
species[11] =
+2.35677352e+00
+8.98459677e-03 * tc[1]
-7.12356269e-06 * tc[2]
+2.45919022e-09 * tc[3]
-1.43699548e-13 * tc[4]
;} else {
// species 7: O2
species[7] =
+3.28253784e+00
+1.48308754e-03 * tc[1]
-7.57966669e-07 * tc[2]
+2.09470555e-10 * tc[3]
-2.16717794e-14 * tc[4]
;
// species 8: N2
species[8] =
+2.92664000e+00
+1.48797700e-03 * tc[1]
-5.68476100e-07 * tc[2]
+1.00970400e-10 * tc[3]
-6.75335100e-15 * tc[4]
;
// species 10: H2O
species[10] =
+2.67214561e+00
+3.05629289e-03 * tc[1]
-8.73026011e-07 * tc[2]
+1.20099639e-10 * tc[3]
-6.39161787e-15 * tc[4]
;
// species 11: CO2
species[11] =
+3.85746029e+00
+4.41437026e-03 * tc[1]
-2.21481404e-06 * tc[2]
+5.23490188e-10 * tc[3]
-4.72084164e-14 * tc[4]
;}

// species with midpoint at T=1800 kelvin
if (T < 1800) {
// species 9: NC16H34
species[9] =
-3.64057977e+00
+1.96191732e-01 * tc[1]
-1.26758759e-04 * tc[2]
+4.13645475e-08 * tc[3]
-5.40570458e-12 * tc[4]
;} else {
// species 9: NC16H34
species[9] =
+4.98210237e+01
+7.73881688e-02 * tc[1]
-2.77557896e-05 * tc[2]
+4.69678103e-09 * tc[3]
-3.12959243e-13 * tc[4]
;}
}

// compute the g/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void gibbs(amrex::Real * species, const amrex::Real *  tc)
{

// temperature
const amrex::Real T = tc[1];
const amrex::Real invT = 1.0 / T;

// species with midpoint at T=1390 kelvin
if (T < 1390) {
// species 0: NC5H12
species[0] =
-1.989349340000000e+04 * invT
-3.027259482200000e+01
+7.367661220000000e-01 * tc[0]
-3.036004865000000e-02 * tc[1]
+5.959879350000000e-06 * tc[2]
-8.742253500000000e-10 * tc[3]
+6.074365750000000e-14 * tc[4]
;
// species 1: NC6H14
species[1] =
-2.278038620000000e+04 * invT
-3.327658598400000e+01
+9.696061840000000e-01 * tc[0]
-3.645428040000000e-02 * tc[1]
+7.314231983333333e-06 * tc[2]
-1.102606725000000e-09 * tc[3]
+7.921871150000000e-14 * tc[4]
;} else {
// species 0: NC5H12
species[0] =
-2.605195430000000e+04 * invT
+7.606230600000001e+01
-1.572576030000000e+01 * tc[0]
-1.305430225000000e-02 * tc[1]
+1.484951660000000e-06 * tc[2]
-1.150852066666667e-10 * tc[3]
+4.001482680000000e-15 * tc[4]
;
// species 1: NC6H14
species[1] =
-3.016287390000000e+04 * invT
+9.524734720000001e+01
-1.896341170000000e+01 * tc[0]
-1.522401020000000e-02 * tc[1]
+1.729913816666667e-06 * tc[2]
-1.339795475000000e-10 * tc[3]
+4.656348640000000e-15 * tc[4]
;}

// species with midpoint at T=1391 kelvin
if (T < 1391) {
// species 2: NC8H18
species[2] =
-2.853956410000000e+04 * invT
-3.987498186000000e+01
+1.542184060000000e+00 * tc[0]
-4.890560315000000e-02 * tc[1]
+1.015530596666667e-05 * tc[2]
-1.600046591666667e-09 * tc[3]
+1.214981250000000e-13 * tc[4]
;
// species 4: NC10H22
species[4] =
-3.430218630000000e+04 * invT
-4.631018369000000e+01
+2.084169690000000e+00 * tc[0]
-6.126750600000000e-02 * tc[1]
+1.294692898333333e-05 * tc[2]
-2.081957308333334e-09 * tc[3]
+1.617740190000000e-13 * tc[4]
;
// species 5: NC7H16
species[5] =
-2.565865650000000e+04 * invT
-3.664165307000000e+01
+1.268361870000000e+00 * tc[0]
-4.271779100000000e-02 * tc[1]
+8.755779766666667e-06 * tc[2]
-1.357881008333333e-09 * tc[3]
+1.011974625000000e-13 * tc[4]
;
// species 6: NC12H26
species[6] =
-4.006542530000000e+04 * invT
-5.272127854000000e+01
+2.621815940000000e+00 * tc[0]
-7.361885550000000e-02 * tc[1]
+1.573283785000000e-05 * tc[2]
-2.562010566666667e-09 * tc[3]
+2.018011150000000e-13 * tc[4]
;} else {
// species 2: NC8H18
species[2] =
-3.839627550000000e+04 * invT
+1.338321134000000e+02
-2.547101940000000e+01 * tc[0]
-1.954435185000000e-02 * tc[1]
+2.217312950000000e-06 * tc[2]
-1.715562725000000e-10 * tc[3]
+5.958358700000000e-15 * tc[4]
;
// species 4: NC10H22
species[4] =
-4.663928400000000e+04 * invT
+1.724923449000000e+02
-3.198822390000000e+01 * tc[0]
-2.386224610000000e-02 * tc[1]
+2.704606516666666e-06 * tc[2]
-2.091360491666667e-10 * tc[3]
+7.260788599999999e-15 * tc[4]
;
// species 5: NC7H16
species[5] =
-3.427600810000000e+04 * invT
+1.145189165000000e+02
-2.221489690000000e+01 * tc[0]
-1.738378750000000e-02 * tc[1]
+1.973452150000000e-06 * tc[2]
-1.527487316666667e-10 * tc[3]
+5.306513300000000e-15 * tc[4]
;
// species 6: NC12H26
species[6] =
-5.488434650000000e+04 * invT
+2.111804257000000e+02
-3.850950370000000e+01 * tc[0]
-2.817750240000000e-02 * tc[1]
+3.191553333333333e-06 * tc[2]
-2.466873850000000e-10 * tc[3]
+8.562207500000000e-15 * tc[4]
;}

// species with midpoint at T=1396 kelvin
if (T < 1396) {
// species 3: IC8H18
species[3] =
-2.994468750000000e+04 * invT
-4.916085903000000e+01
+4.208688930000000e+00 * tc[0]
-5.572029050000000e-02 * tc[1]
+1.318910970000000e-05 * tc[2]
-2.436718683333333e-09 * tc[3]
+2.218715955000000e-13 * tc[4]
;} else {
// species 3: IC8H18
species[3] =
-4.079581770000000e+04 * invT
+1.504148540000000e+02
-2.713735900000000e+01 * tc[0]
-1.895024450000000e-02 * tc[1]
+2.157289300000000e-06 * tc[2]
-1.673003100000000e-10 * tc[3]
+5.820029000000000e-15 * tc[4]
;}

// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 7: O2
species[7] =
-1.063943560000000e+03 * invT
+1.247806300000001e-01
-3.782456360000000e+00 * tc[0]
+1.498367080000000e-03 * tc[1]
-1.641217001666667e-06 * tc[2]
+8.067745908333334e-10 * tc[3]
-1.621864185000000e-13 * tc[4]
;
// species 8: N2
species[8] =
-1.020900000000000e+03 * invT
-6.516950000000001e-01
-3.298677000000000e+00 * tc[0]
-7.041200000000000e-04 * tc[1]
+6.605369999999999e-07 * tc[2]
-4.701262500000001e-10 * tc[3]
+1.222427500000000e-13 * tc[4]
;
// species 10: H2O
species[10] =
-3.020811330000000e+04 * invT
+7.966096399999998e-01
-3.386842490000000e+00 * tc[0]
-1.737491230000000e-03 * tc[1]
+1.059116055000000e-06 * tc[2]
-5.807151058333333e-10 * tc[3]
+1.253294235000000e-13 * tc[4]
;
// species 11: CO2
species[11] =
-4.837196970000000e+04 * invT
-7.544278700000000e+00
-2.356773520000000e+00 * tc[0]
-4.492298385000000e-03 * tc[1]
+1.187260448333333e-06 * tc[2]
-2.049325183333333e-10 * tc[3]
+7.184977399999999e-15 * tc[4]
;} else {
// species 7: O2
species[7] =
-1.088457720000000e+03 * invT
-2.170693450000000e+00
-3.282537840000000e+00 * tc[0]
-7.415437700000000e-04 * tc[1]
+1.263277781666667e-07 * tc[2]
-1.745587958333333e-11 * tc[3]
+1.083588970000000e-15 * tc[4]
;
// species 8: N2
species[8] =
-9.227977000000000e+02 * invT
-3.053888000000000e+00
-2.926640000000000e+00 * tc[0]
-7.439885000000000e-04 * tc[1]
+9.474601666666666e-08 * tc[2]
-8.414199999999999e-12 * tc[3]
+3.376675500000000e-16 * tc[4]
;
// species 10: H2O
species[10] =
-2.989920900000000e+04 * invT
-4.190671200000001e+00
-2.672145610000000e+00 * tc[0]
-1.528146445000000e-03 * tc[1]
+1.455043351666667e-07 * tc[2]
-1.000830325000000e-11 * tc[3]
+3.195808935000000e-16 * tc[4]
;
// species 11: CO2
species[11] =
-4.875916600000000e+04 * invT
+1.585822230000000e+00
-3.857460290000000e+00 * tc[0]
-2.207185130000000e-03 * tc[1]
+3.691356733333334e-07 * tc[2]
-4.362418233333334e-11 * tc[3]
+2.360420820000000e-15 * tc[4]
;}

// species with midpoint at T=1800 kelvin
if (T < 1800) {
// species 9: NC16H34
species[9] =
-5.160528910000000e+04 * invT
-6.814307517000000e+01
+3.640579770000000e+00 * tc[0]
-9.809586600000000e-02 * tc[1]
+2.112645983333333e-05 * tc[2]
-3.447045625000000e-09 * tc[3]
+2.702852290000000e-13 * tc[4]
;} else {
// species 9: NC16H34
species[9] =
-7.085146630000000e+04 * invT
+2.746638787000000e+02
-4.982102370000000e+01 * tc[0]
-3.869408440000000e-02 * tc[1]
+4.625964933333333e-06 * tc[2]
-3.913984191666667e-10 * tc[3]
+1.564796215000000e-14 * tc[4]
;}
}

// compute the a/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void helmholtz(amrex::Real * species, const amrex::Real *  tc)
{

// temperature
const amrex::Real T = tc[1];
const amrex::Real invT = 1.0 / T;

// species with midpoint at T=1390 kelvin
if (T < 1390) {
// species 0: NC5H12
species[0] =
-1.98934934e+04 * invT
-3.12725948e+01
+7.36766122e-01 * tc[0]
-3.03600486e-02 * tc[1]
+5.95987935e-06 * tc[2]
-8.74225350e-10 * tc[3]
+6.07436575e-14 * tc[4]
;
// species 1: NC6H14
species[1] =
-2.27803862e+04 * invT
-3.42765860e+01
+9.69606184e-01 * tc[0]
-3.64542804e-02 * tc[1]
+7.31423198e-06 * tc[2]
-1.10260672e-09 * tc[3]
+7.92187115e-14 * tc[4]
;} else {
// species 0: NC5H12
species[0] =
-2.60519543e+04 * invT
+7.50623060e+01
-1.57257603e+01 * tc[0]
-1.30543023e-02 * tc[1]
+1.48495166e-06 * tc[2]
-1.15085207e-10 * tc[3]
+4.00148268e-15 * tc[4]
;
// species 1: NC6H14
species[1] =
-3.01628739e+04 * invT
+9.42473472e+01
-1.89634117e+01 * tc[0]
-1.52240102e-02 * tc[1]
+1.72991382e-06 * tc[2]
-1.33979548e-10 * tc[3]
+4.65634864e-15 * tc[4]
;}

// species with midpoint at T=1391 kelvin
if (T < 1391) {
// species 2: NC8H18
species[2] =
-2.85395641e+04 * invT
-4.08749819e+01
+1.54218406e+00 * tc[0]
-4.89056031e-02 * tc[1]
+1.01553060e-05 * tc[2]
-1.60004659e-09 * tc[3]
+1.21498125e-13 * tc[4]
;
// species 4: NC10H22
species[4] =
-3.43021863e+04 * invT
-4.73101837e+01
+2.08416969e+00 * tc[0]
-6.12675060e-02 * tc[1]
+1.29469290e-05 * tc[2]
-2.08195731e-09 * tc[3]
+1.61774019e-13 * tc[4]
;
// species 5: NC7H16
species[5] =
-2.56586565e+04 * invT
-3.76416531e+01
+1.26836187e+00 * tc[0]
-4.27177910e-02 * tc[1]
+8.75577977e-06 * tc[2]
-1.35788101e-09 * tc[3]
+1.01197462e-13 * tc[4]
;
// species 6: NC12H26
species[6] =
-4.00654253e+04 * invT
-5.37212785e+01
+2.62181594e+00 * tc[0]
-7.36188555e-02 * tc[1]
+1.57328378e-05 * tc[2]
-2.56201057e-09 * tc[3]
+2.01801115e-13 * tc[4]
;} else {
// species 2: NC8H18
species[2] =
-3.83962755e+04 * invT
+1.32832113e+02
-2.54710194e+01 * tc[0]
-1.95443518e-02 * tc[1]
+2.21731295e-06 * tc[2]
-1.71556272e-10 * tc[3]
+5.95835870e-15 * tc[4]
;
// species 4: NC10H22
species[4] =
-4.66392840e+04 * invT
+1.71492345e+02
-3.19882239e+01 * tc[0]
-2.38622461e-02 * tc[1]
+2.70460652e-06 * tc[2]
-2.09136049e-10 * tc[3]
+7.26078860e-15 * tc[4]
;
// species 5: NC7H16
species[5] =
-3.42760081e+04 * invT
+1.13518917e+02
-2.22148969e+01 * tc[0]
-1.73837875e-02 * tc[1]
+1.97345215e-06 * tc[2]
-1.52748732e-10 * tc[3]
+5.30651330e-15 * tc[4]
;
// species 6: NC12H26
species[6] =
-5.48843465e+04 * invT
+2.10180426e+02
-3.85095037e+01 * tc[0]
-2.81775024e-02 * tc[1]
+3.19155333e-06 * tc[2]
-2.46687385e-10 * tc[3]
+8.56220750e-15 * tc[4]
;}

// species with midpoint at T=1396 kelvin
if (T < 1396) {
// species 3: IC8H18
species[3] =
-2.99446875e+04 * invT
-5.01608590e+01
+4.20868893e+00 * tc[0]
-5.57202905e-02 * tc[1]
+1.31891097e-05 * tc[2]
-2.43671868e-09 * tc[3]
+2.21871595e-13 * tc[4]
;} else {
// species 3: IC8H18
species[3] =
-4.07958177e+04 * invT
+1.49414854e+02
-2.71373590e+01 * tc[0]
-1.89502445e-02 * tc[1]
+2.15728930e-06 * tc[2]
-1.67300310e-10 * tc[3]
+5.82002900e-15 * tc[4]
;}

// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 7: O2
species[7] =
-1.06394356e+03 * invT
-8.75219370e-01
-3.78245636e+00 * tc[0]
+1.49836708e-03 * tc[1]
-1.64121700e-06 * tc[2]
+8.06774591e-10 * tc[3]
-1.62186418e-13 * tc[4]
;
// species 8: N2
species[8] =
-1.02090000e+03 * invT
-1.65169500e+00
-3.29867700e+00 * tc[0]
-7.04120000e-04 * tc[1]
+6.60537000e-07 * tc[2]
-4.70126250e-10 * tc[3]
+1.22242750e-13 * tc[4]
;
// species 10: H2O
species[10] =
-3.02081133e+04 * invT
-2.03390360e-01
-3.38684249e+00 * tc[0]
-1.73749123e-03 * tc[1]
+1.05911606e-06 * tc[2]
-5.80715106e-10 * tc[3]
+1.25329424e-13 * tc[4]
;
// species 11: CO2
species[11] =
-4.83719697e+04 * invT
-8.54427870e+00
-2.35677352e+00 * tc[0]
-4.49229839e-03 * tc[1]
+1.18726045e-06 * tc[2]
-2.04932518e-10 * tc[3]
+7.18497740e-15 * tc[4]
;} else {
// species 7: O2
species[7] =
-1.08845772e+03 * invT
-3.17069345e+00
-3.28253784e+00 * tc[0]
-7.41543770e-04 * tc[1]
+1.26327778e-07 * tc[2]
-1.74558796e-11 * tc[3]
+1.08358897e-15 * tc[4]
;
// species 8: N2
species[8] =
-9.22797700e+02 * invT
-4.05388800e+00
-2.92664000e+00 * tc[0]
-7.43988500e-04 * tc[1]
+9.47460167e-08 * tc[2]
-8.41420000e-12 * tc[3]
+3.37667550e-16 * tc[4]
;
// species 10: H2O
species[10] =
-2.98992090e+04 * invT
-5.19067120e+00
-2.67214561e+00 * tc[0]
-1.52814644e-03 * tc[1]
+1.45504335e-07 * tc[2]
-1.00083033e-11 * tc[3]
+3.19580894e-16 * tc[4]
;
// species 11: CO2
species[11] =
-4.87591660e+04 * invT
+5.85822230e-01
-3.85746029e+00 * tc[0]
-2.20718513e-03 * tc[1]
+3.69135673e-07 * tc[2]
-4.36241823e-11 * tc[3]
+2.36042082e-15 * tc[4]
;}

// species with midpoint at T=1800 kelvin
if (T < 1800) {
// species 9: NC16H34
species[9] =
-5.16052891e+04 * invT
-6.91430752e+01
+3.64057977e+00 * tc[0]
-9.80958660e-02 * tc[1]
+2.11264598e-05 * tc[2]
-3.44704562e-09 * tc[3]
+2.70285229e-13 * tc[4]
;} else {
// species 9: NC16H34
species[9] =
-7.08514663e+04 * invT
+2.73663879e+02
-4.98210237e+01 * tc[0]
-3.86940844e-02 * tc[1]
+4.62596493e-06 * tc[2]
-3.91398419e-10 * tc[3]
+1.56479622e-14 * tc[4]
;}
}

// compute the e/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesInternalEnergy(amrex::Real * species, const amrex::Real *  tc)
{

// temperature
const amrex::Real T = tc[1];
const amrex::Real invT = 1.0 / T;

// species with midpoint at T=1390 kelvin
if (T < 1390) {
// species 0: NC5H12
species[0] =
-1.73676612e+00
+3.03600486e-02 * tc[1]
-1.19197587e-05 * tc[2]
+2.62267605e-09 * tc[3]
-2.42974630e-13 * tc[4]
-1.98934934e+04 * invT
;
// species 1: NC6H14
species[1] =
-1.96960618e+00
+3.64542804e-02 * tc[1]
-1.46284640e-05 * tc[2]
+3.30782017e-09 * tc[3]
-3.16874846e-13 * tc[4]
-2.27803862e+04 * invT
;} else {
// species 0: NC5H12
species[0] =
+1.47257603e+01
+1.30543023e-02 * tc[1]
-2.96990332e-06 * tc[2]
+3.45255620e-10 * tc[3]
-1.60059307e-14 * tc[4]
-2.60519543e+04 * invT
;
// species 1: NC6H14
species[1] =
+1.79634117e+01
+1.52240102e-02 * tc[1]
-3.45982763e-06 * tc[2]
+4.01938643e-10 * tc[3]
-1.86253946e-14 * tc[4]
-3.01628739e+04 * invT
;}

// species with midpoint at T=1391 kelvin
if (T < 1391) {
// species 2: NC8H18
species[2] =
-2.54218406e+00
+4.89056031e-02 * tc[1]
-2.03106119e-05 * tc[2]
+4.80013978e-09 * tc[3]
-4.85992500e-13 * tc[4]
-2.85395641e+04 * invT
;
// species 4: NC10H22
species[4] =
-3.08416969e+00
+6.12675060e-02 * tc[1]
-2.58938580e-05 * tc[2]
+6.24587193e-09 * tc[3]
-6.47096076e-13 * tc[4]
-3.43021863e+04 * invT
;
// species 5: NC7H16
species[5] =
-2.26836187e+00
+4.27177910e-02 * tc[1]
-1.75115595e-05 * tc[2]
+4.07364302e-09 * tc[3]
-4.04789850e-13 * tc[4]
-2.56586565e+04 * invT
;
// species 6: NC12H26
species[6] =
-3.62181594e+00
+7.36188555e-02 * tc[1]
-3.14656757e-05 * tc[2]
+7.68603170e-09 * tc[3]
-8.07204460e-13 * tc[4]
-4.00654253e+04 * invT
;} else {
// species 2: NC8H18
species[2] =
+2.44710194e+01
+1.95443518e-02 * tc[1]
-4.43462590e-06 * tc[2]
+5.14668817e-10 * tc[3]
-2.38334348e-14 * tc[4]
-3.83962755e+04 * invT
;
// species 4: NC10H22
species[4] =
+3.09882239e+01
+2.38622461e-02 * tc[1]
-5.40921303e-06 * tc[2]
+6.27408147e-10 * tc[3]
-2.90431544e-14 * tc[4]
-4.66392840e+04 * invT
;
// species 5: NC7H16
species[5] =
+2.12148969e+01
+1.73837875e-02 * tc[1]
-3.94690430e-06 * tc[2]
+4.58246195e-10 * tc[3]
-2.12260532e-14 * tc[4]
-3.42760081e+04 * invT
;
// species 6: NC12H26
species[6] =
+3.75095037e+01
+2.81775024e-02 * tc[1]
-6.38310667e-06 * tc[2]
+7.40062155e-10 * tc[3]
-3.42488300e-14 * tc[4]
-5.48843465e+04 * invT
;}

// species with midpoint at T=1396 kelvin
if (T < 1396) {
// species 3: IC8H18
species[3] =
-5.20868893e+00
+5.57202905e-02 * tc[1]
-2.63782194e-05 * tc[2]
+7.31015605e-09 * tc[3]
-8.87486382e-13 * tc[4]
-2.99446875e+04 * invT
;} else {
// species 3: IC8H18
species[3] =
+2.61373590e+01
+1.89502445e-02 * tc[1]
-4.31457860e-06 * tc[2]
+5.01900930e-10 * tc[3]
-2.32801160e-14 * tc[4]
-4.07958177e+04 * invT
;}

// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 7: O2
species[7] =
+2.78245636e+00
-1.49836708e-03 * tc[1]
+3.28243400e-06 * tc[2]
-2.42032377e-09 * tc[3]
+6.48745674e-13 * tc[4]
-1.06394356e+03 * invT
;
// species 8: N2
species[8] =
+2.29867700e+00
+7.04120000e-04 * tc[1]
-1.32107400e-06 * tc[2]
+1.41037875e-09 * tc[3]
-4.88971000e-13 * tc[4]
-1.02090000e+03 * invT
;
// species 10: H2O
species[10] =
+2.38684249e+00
+1.73749123e-03 * tc[1]
-2.11823211e-06 * tc[2]
+1.74214532e-09 * tc[3]
-5.01317694e-13 * tc[4]
-3.02081133e+04 * invT
;
// species 11: CO2
species[11] =
+1.35677352e+00
+4.49229839e-03 * tc[1]
-2.37452090e-06 * tc[2]
+6.14797555e-10 * tc[3]
-2.87399096e-14 * tc[4]
-4.83719697e+04 * invT
;} else {
// species 7: O2
species[7] =
+2.28253784e+00
+7.41543770e-04 * tc[1]
-2.52655556e-07 * tc[2]
+5.23676387e-11 * tc[3]
-4.33435588e-15 * tc[4]
-1.08845772e+03 * invT
;
// species 8: N2
species[8] =
+1.92664000e+00
+7.43988500e-04 * tc[1]
-1.89492033e-07 * tc[2]
+2.52426000e-11 * tc[3]
-1.35067020e-15 * tc[4]
-9.22797700e+02 * invT
;
// species 10: H2O
species[10] =
+1.67214561e+00
+1.52814644e-03 * tc[1]
-2.91008670e-07 * tc[2]
+3.00249098e-11 * tc[3]
-1.27832357e-15 * tc[4]
-2.98992090e+04 * invT
;
// species 11: CO2
species[11] =
+2.85746029e+00
+2.20718513e-03 * tc[1]
-7.38271347e-07 * tc[2]
+1.30872547e-10 * tc[3]
-9.44168328e-15 * tc[4]
-4.87591660e+04 * invT
;}

// species with midpoint at T=1800 kelvin
if (T < 1800) {
// species 9: NC16H34
species[9] =
-4.64057977e+00
+9.80958660e-02 * tc[1]
-4.22529197e-05 * tc[2]
+1.03411369e-08 * tc[3]
-1.08114092e-12 * tc[4]
-5.16052891e+04 * invT
;} else {
// species 9: NC16H34
species[9] =
+4.88210237e+01
+3.86940844e-02 * tc[1]
-9.25192987e-06 * tc[2]
+1.17419526e-09 * tc[3]
-6.25918486e-14 * tc[4]
-7.08514663e+04 * invT
;}
}

// compute the h/(RT) at the given temperature (Eq 20)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesEnthalpy(amrex::Real * species, const amrex::Real *  tc)
{

// temperature
const amrex::Real T = tc[1];
const amrex::Real invT = 1.0 / T;

// species with midpoint at T=1390 kelvin
if (T < 1390) {
// species 0: NC5H12
species[0] =
-7.36766122e-01
+3.03600486e-02 * tc[1]
-1.19197587e-05 * tc[2]
+2.62267605e-09 * tc[3]
-2.42974630e-13 * tc[4]
-1.98934934e+04 * invT
;
// species 1: NC6H14
species[1] =
-9.69606184e-01
+3.64542804e-02 * tc[1]
-1.46284640e-05 * tc[2]
+3.30782017e-09 * tc[3]
-3.16874846e-13 * tc[4]
-2.27803862e+04 * invT
;} else {
// species 0: NC5H12
species[0] =
+1.57257603e+01
+1.30543023e-02 * tc[1]
-2.96990332e-06 * tc[2]
+3.45255620e-10 * tc[3]
-1.60059307e-14 * tc[4]
-2.60519543e+04 * invT
;
// species 1: NC6H14
species[1] =
+1.89634117e+01
+1.52240102e-02 * tc[1]
-3.45982763e-06 * tc[2]
+4.01938643e-10 * tc[3]
-1.86253946e-14 * tc[4]
-3.01628739e+04 * invT
;}

// species with midpoint at T=1391 kelvin
if (T < 1391) {
// species 2: NC8H18
species[2] =
-1.54218406e+00
+4.89056031e-02 * tc[1]
-2.03106119e-05 * tc[2]
+4.80013978e-09 * tc[3]
-4.85992500e-13 * tc[4]
-2.85395641e+04 * invT
;
// species 4: NC10H22
species[4] =
-2.08416969e+00
+6.12675060e-02 * tc[1]
-2.58938580e-05 * tc[2]
+6.24587193e-09 * tc[3]
-6.47096076e-13 * tc[4]
-3.43021863e+04 * invT
;
// species 5: NC7H16
species[5] =
-1.26836187e+00
+4.27177910e-02 * tc[1]
-1.75115595e-05 * tc[2]
+4.07364302e-09 * tc[3]
-4.04789850e-13 * tc[4]
-2.56586565e+04 * invT
;
// species 6: NC12H26
species[6] =
-2.62181594e+00
+7.36188555e-02 * tc[1]
-3.14656757e-05 * tc[2]
+7.68603170e-09 * tc[3]
-8.07204460e-13 * tc[4]
-4.00654253e+04 * invT
;} else {
// species 2: NC8H18
species[2] =
+2.54710194e+01
+1.95443518e-02 * tc[1]
-4.43462590e-06 * tc[2]
+5.14668817e-10 * tc[3]
-2.38334348e-14 * tc[4]
-3.83962755e+04 * invT
;
// species 4: NC10H22
species[4] =
+3.19882239e+01
+2.38622461e-02 * tc[1]
-5.40921303e-06 * tc[2]
+6.27408147e-10 * tc[3]
-2.90431544e-14 * tc[4]
-4.66392840e+04 * invT
;
// species 5: NC7H16
species[5] =
+2.22148969e+01
+1.73837875e-02 * tc[1]
-3.94690430e-06 * tc[2]
+4.58246195e-10 * tc[3]
-2.12260532e-14 * tc[4]
-3.42760081e+04 * invT
;
// species 6: NC12H26
species[6] =
+3.85095037e+01
+2.81775024e-02 * tc[1]
-6.38310667e-06 * tc[2]
+7.40062155e-10 * tc[3]
-3.42488300e-14 * tc[4]
-5.48843465e+04 * invT
;}

// species with midpoint at T=1396 kelvin
if (T < 1396) {
// species 3: IC8H18
species[3] =
-4.20868893e+00
+5.57202905e-02 * tc[1]
-2.63782194e-05 * tc[2]
+7.31015605e-09 * tc[3]
-8.87486382e-13 * tc[4]
-2.99446875e+04 * invT
;} else {
// species 3: IC8H18
species[3] =
+2.71373590e+01
+1.89502445e-02 * tc[1]
-4.31457860e-06 * tc[2]
+5.01900930e-10 * tc[3]
-2.32801160e-14 * tc[4]
-4.07958177e+04 * invT
;}

// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 7: O2
species[7] =
+3.78245636e+00
-1.49836708e-03 * tc[1]
+3.28243400e-06 * tc[2]
-2.42032377e-09 * tc[3]
+6.48745674e-13 * tc[4]
-1.06394356e+03 * invT
;
// species 8: N2
species[8] =
+3.29867700e+00
+7.04120000e-04 * tc[1]
-1.32107400e-06 * tc[2]
+1.41037875e-09 * tc[3]
-4.88971000e-13 * tc[4]
-1.02090000e+03 * invT
;
// species 10: H2O
species[10] =
+3.38684249e+00
+1.73749123e-03 * tc[1]
-2.11823211e-06 * tc[2]
+1.74214532e-09 * tc[3]
-5.01317694e-13 * tc[4]
-3.02081133e+04 * invT
;
// species 11: CO2
species[11] =
+2.35677352e+00
+4.49229839e-03 * tc[1]
-2.37452090e-06 * tc[2]
+6.14797555e-10 * tc[3]
-2.87399096e-14 * tc[4]
-4.83719697e+04 * invT
;} else {
// species 7: O2
species[7] =
+3.28253784e+00
+7.41543770e-04 * tc[1]
-2.52655556e-07 * tc[2]
+5.23676387e-11 * tc[3]
-4.33435588e-15 * tc[4]
-1.08845772e+03 * invT
;
// species 8: N2
species[8] =
+2.92664000e+00
+7.43988500e-04 * tc[1]
-1.89492033e-07 * tc[2]
+2.52426000e-11 * tc[3]
-1.35067020e-15 * tc[4]
-9.22797700e+02 * invT
;
// species 10: H2O
species[10] =
+2.67214561e+00
+1.52814644e-03 * tc[1]
-2.91008670e-07 * tc[2]
+3.00249098e-11 * tc[3]
-1.27832357e-15 * tc[4]
-2.98992090e+04 * invT
;
// species 11: CO2
species[11] =
+3.85746029e+00
+2.20718513e-03 * tc[1]
-7.38271347e-07 * tc[2]
+1.30872547e-10 * tc[3]
-9.44168328e-15 * tc[4]
-4.87591660e+04 * invT
;}

// species with midpoint at T=1800 kelvin
if (T < 1800) {
// species 9: NC16H34
species[9] =
-3.64057977e+00
+9.80958660e-02 * tc[1]
-4.22529197e-05 * tc[2]
+1.03411369e-08 * tc[3]
-1.08114092e-12 * tc[4]
-5.16052891e+04 * invT
;} else {
// species 9: NC16H34
species[9] =
+4.98210237e+01
+3.86940844e-02 * tc[1]
-9.25192987e-06 * tc[2]
+1.17419526e-09 * tc[3]
-6.25918486e-14 * tc[4]
-7.08514663e+04 * invT
;}
}

// compute the S/R at the given temperature (Eq 21)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesEntropy(amrex::Real * species, const amrex::Real *  tc)
{

// temperature
const amrex::Real T = tc[1];

// species with midpoint at T=1390 kelvin
if (T < 1390) {
// species 0: NC5H12
species[0] =
-7.36766122e-01 * tc[0]
+6.07200973e-02 * tc[1]
-1.78796381e-05 * tc[2]
+3.49690140e-09 * tc[3]
-3.03718288e-13 * tc[4]
+2.95358287e+01
;
// species 1: NC6H14
species[1] =
-9.69606184e-01 * tc[0]
+7.29085608e-02 * tc[1]
-2.19426959e-05 * tc[2]
+4.41042690e-09 * tc[3]
-3.96093557e-13 * tc[4]
+3.23069798e+01
;} else {
// species 0: NC5H12
species[0] =
+1.57257603e+01 * tc[0]
+2.61086045e-02 * tc[1]
-4.45485498e-06 * tc[2]
+4.60340827e-10 * tc[3]
-2.00074134e-14 * tc[4]
-6.03365457e+01
;
// species 1: NC6H14
species[1] =
+1.89634117e+01 * tc[0]
+3.04480204e-02 * tc[1]
-5.18974145e-06 * tc[2]
+5.35918190e-10 * tc[3]
-2.32817432e-14 * tc[4]
-7.62839355e+01
;}

// species with midpoint at T=1391 kelvin
if (T < 1391) {
// species 2: NC8H18
species[2] =
-1.54218406e+00 * tc[0]
+9.78112063e-02 * tc[1]
-3.04659179e-05 * tc[2]
+6.40018637e-09 * tc[3]
-6.07490625e-13 * tc[4]
+3.83327978e+01
;
// species 4: NC10H22
species[4] =
-2.08416969e+00 * tc[0]
+1.22535012e-01 * tc[1]
-3.88407869e-05 * tc[2]
+8.32782923e-09 * tc[3]
-8.08870095e-13 * tc[4]
+4.42260140e+01
;
// species 5: NC7H16
species[5] =
-1.26836187e+00 * tc[0]
+8.54355820e-02 * tc[1]
-2.62673393e-05 * tc[2]
+5.43152403e-09 * tc[3]
-5.05987313e-13 * tc[4]
+3.53732912e+01
;
// species 6: NC12H26
species[6] =
-2.62181594e+00 * tc[0]
+1.47237711e-01 * tc[1]
-4.71985136e-05 * tc[2]
+1.02480423e-08 * tc[3]
-1.00900558e-12 * tc[4]
+5.00994626e+01
;} else {
// species 2: NC8H18
species[2] =
+2.54710194e+01 * tc[0]
+3.90887037e-02 * tc[1]
-6.65193885e-06 * tc[2]
+6.86225090e-10 * tc[3]
-2.97917935e-14 * tc[4]
-1.08361094e+02
;
// species 4: NC10H22
species[4] =
+3.19882239e+01 * tc[0]
+4.77244922e-02 * tc[1]
-8.11381955e-06 * tc[2]
+8.36544197e-10 * tc[3]
-3.63039430e-14 * tc[4]
-1.40504121e+02
;
// species 5: NC7H16
species[5] =
+2.22148969e+01 * tc[0]
+3.47675750e-02 * tc[1]
-5.92035645e-06 * tc[2]
+6.10994927e-10 * tc[3]
-2.65325665e-14 * tc[4]
-9.23040196e+01
;
// species 6: NC12H26
species[6] =
+3.85095037e+01 * tc[0]
+5.63550048e-02 * tc[1]
-9.57466000e-06 * tc[2]
+9.86749540e-10 * tc[3]
-4.28110375e-14 * tc[4]
-1.72670922e+02
;}

// species with midpoint at T=1396 kelvin
if (T < 1396) {
// species 3: IC8H18
species[3] =
-4.20868893e+00 * tc[0]
+1.11440581e-01 * tc[1]
-3.95673291e-05 * tc[2]
+9.74687473e-09 * tc[3]
-1.10935798e-12 * tc[4]
+4.49521701e+01
;} else {
// species 3: IC8H18
species[3] =
+2.71373590e+01 * tc[0]
+3.79004890e-02 * tc[1]
-6.47186790e-06 * tc[2]
+6.69201240e-10 * tc[3]
-2.91001450e-14 * tc[4]
-1.23277495e+02
;}

// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 7: O2
species[7] =
+3.78245636e+00 * tc[0]
-2.99673416e-03 * tc[1]
+4.92365101e-06 * tc[2]
-3.22709836e-09 * tc[3]
+8.10932092e-13 * tc[4]
+3.65767573e+00
;
// species 8: N2
species[8] =
+3.29867700e+00 * tc[0]
+1.40824000e-03 * tc[1]
-1.98161100e-06 * tc[2]
+1.88050500e-09 * tc[3]
-6.11213750e-13 * tc[4]
+3.95037200e+00
;
// species 10: H2O
species[10] =
+3.38684249e+00 * tc[0]
+3.47498246e-03 * tc[1]
-3.17734817e-06 * tc[2]
+2.32286042e-09 * tc[3]
-6.26647117e-13 * tc[4]
+2.59023285e+00
;
// species 11: CO2
species[11] =
+2.35677352e+00 * tc[0]
+8.98459677e-03 * tc[1]
-3.56178134e-06 * tc[2]
+8.19730073e-10 * tc[3]
-3.59248870e-14 * tc[4]
+9.90105222e+00
;} else {
// species 7: O2
species[7] =
+3.28253784e+00 * tc[0]
+1.48308754e-03 * tc[1]
-3.78983334e-07 * tc[2]
+6.98235183e-11 * tc[3]
-5.41794485e-15 * tc[4]
+5.45323129e+00
;
// species 8: N2
species[8] =
+2.92664000e+00 * tc[0]
+1.48797700e-03 * tc[1]
-2.84238050e-07 * tc[2]
+3.36568000e-11 * tc[3]
-1.68833775e-15 * tc[4]
+5.98052800e+00
;
// species 10: H2O
species[10] =
+2.67214561e+00 * tc[0]
+3.05629289e-03 * tc[1]
-4.36513005e-07 * tc[2]
+4.00332130e-11 * tc[3]
-1.59790447e-15 * tc[4]
+6.86281681e+00
;
// species 11: CO2
species[11] =
+3.85746029e+00 * tc[0]
+4.41437026e-03 * tc[1]
-1.10740702e-06 * tc[2]
+1.74496729e-10 * tc[3]
-1.18021041e-14 * tc[4]
+2.27163806e+00
;}

// species with midpoint at T=1800 kelvin
if (T < 1800) {
// species 9: NC16H34
species[9] =
-3.64057977e+00 * tc[0]
+1.96191732e-01 * tc[1]
-6.33793795e-05 * tc[2]
+1.37881825e-08 * tc[3]
-1.35142614e-12 * tc[4]
+6.45024954e+01
;} else {
// species 9: NC16H34
species[9] =
+4.98210237e+01 * tc[0]
+7.73881688e-02 * tc[1]
-1.38778948e-05 * tc[2]
+1.56559368e-09 * tc[3]
-7.82398108e-14 * tc[4]
-2.24842855e+02
;}
}

// Returns the mean specific heat at CP (Eq. 33)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPBL(const amrex::Real T, const amrex::Real x[], amrex::Real& cpbl)
{
amrex::Real result = 0; 
amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
amrex::Real cpor[12]; //  temporary storage
cp_R(cpor, tc);

// perform dot product
for (int id = 0; id < 12; ++id) {
result += x[id]*cpor[id];
}

cpbl = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CP (Eq. 34)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPBS(const amrex::Real T, const amrex::Real y[], amrex::Real& cpbs)
{
amrex::Real result = 0.0; 
const amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache

// compute Cp/R at the given temperature

// species with midpoint at T=1390 kelvin
if (tT < 1390) {
// species 0: NC5H12
result += y[0] * (
-7.36766122e-01
+6.07200973e-02 * tc[1]
-3.57592761e-05 * tc[2]
+1.04907042e-08 * tc[3]
-1.21487315e-12 * tc[4]
)* 0.0138598217626921
;
// species 1: NC6H14
result += y[1] * (
-9.69606184e-01
+7.29085608e-02 * tc[1]
-4.38853919e-05 * tc[2]
+1.32312807e-08 * tc[3]
-1.58437423e-12 * tc[4]
)* 0.0116038896238019
;} else {
// species 0: NC5H12
result += y[0] * (
+1.57257603e+01
+2.61086045e-02 * tc[1]
-8.90970996e-06 * tc[2]
+1.38102248e-09 * tc[3]
-8.00296536e-14 * tc[4]
)* 0.0138598217626921
;
// species 1: NC6H14
result += y[1] * (
+1.89634117e+01
+3.04480204e-02 * tc[1]
-1.03794829e-05 * tc[2]
+1.60775457e-09 * tc[3]
-9.31269728e-14 * tc[4]
)* 0.0116038896238019
;}

// species with midpoint at T=1391 kelvin
if (tT < 1391) {
// species 2: NC8H18
result += y[2] * (
-1.54218406e+00
+9.78112063e-02 * tc[1]
-6.09318358e-05 * tc[2]
+1.92005591e-08 * tc[3]
-2.42996250e-12 * tc[4]
)* 0.0087541144337839
;
// species 4: NC10H22
result += y[4] * (
-2.08416969e+00
+1.22535012e-01 * tc[1]
-7.76815739e-05 * tc[2]
+2.49834877e-08 * tc[3]
-3.23548038e-12 * tc[4]
)* 0.0070280983371519
;
// species 5: NC7H16
result += y[5] * (
-1.26836187e+00
+8.54355820e-02 * tc[1]
-5.25346786e-05 * tc[2]
+1.62945721e-08 * tc[3]
-2.02394925e-12 * tc[4]
)* 0.0099795419390250
;
// species 6: NC12H26
result += y[6] * (
-2.62181594e+00
+1.47237711e-01 * tc[1]
-9.43970271e-05 * tc[2]
+3.07441268e-08 * tc[3]
-4.03602230e-12 * tc[4]
)* 0.0058706117177410
;} else {
// species 2: NC8H18
result += y[2] * (
+2.54710194e+01
+3.90887037e-02 * tc[1]
-1.33038777e-05 * tc[2]
+2.05867527e-09 * tc[3]
-1.19167174e-13 * tc[4]
)* 0.0087541144337839
;
// species 4: NC10H22
result += y[4] * (
+3.19882239e+01
+4.77244922e-02 * tc[1]
-1.62276391e-05 * tc[2]
+2.50963259e-09 * tc[3]
-1.45215772e-13 * tc[4]
)* 0.0070280983371519
;
// species 5: NC7H16
result += y[5] * (
+2.22148969e+01
+3.47675750e-02 * tc[1]
-1.18407129e-05 * tc[2]
+1.83298478e-09 * tc[3]
-1.06130266e-13 * tc[4]
)* 0.0099795419390250
;
// species 6: NC12H26
result += y[6] * (
+3.85095037e+01
+5.63550048e-02 * tc[1]
-1.91493200e-05 * tc[2]
+2.96024862e-09 * tc[3]
-1.71244150e-13 * tc[4]
)* 0.0058706117177410
;}

// species with midpoint at T=1396 kelvin
if (tT < 1396) {
// species 3: IC8H18
result += y[3] * (
-4.20868893e+00
+1.11440581e-01 * tc[1]
-7.91346582e-05 * tc[2]
+2.92406242e-08 * tc[3]
-4.43743191e-12 * tc[4]
)* 0.0087541144337839
;} else {
// species 3: IC8H18
result += y[3] * (
+2.71373590e+01
+3.79004890e-02 * tc[1]
-1.29437358e-05 * tc[2]
+2.00760372e-09 * tc[3]
-1.16400580e-13 * tc[4]
)* 0.0087541144337839
;}

// species with midpoint at T=1000 kelvin
if (tT < 1000) {
// species 7: O2
result += y[7] * (
+3.78245636e+00
-2.99673416e-03 * tc[1]
+9.84730201e-06 * tc[2]
-9.68129509e-09 * tc[3]
+3.24372837e-12 * tc[4]
)* 0.0312519532470779
;
// species 8: N2
result += y[8] * (
+3.29867700e+00
+1.40824000e-03 * tc[1]
-3.96322200e-06 * tc[2]
+5.64151500e-09 * tc[3]
-2.44485500e-12 * tc[4]
)* 0.0356964374955379
;
// species 10: H2O
result += y[10] * (
+3.38684249e+00
+3.47498246e-03 * tc[1]
-6.35469633e-06 * tc[2]
+6.96858127e-09 * tc[3]
-2.50658847e-12 * tc[4]
)* 0.0555092978073827
;
// species 11: CO2
result += y[11] * (
+2.35677352e+00
+8.98459677e-03 * tc[1]
-7.12356269e-06 * tc[2]
+2.45919022e-09 * tc[3]
-1.43699548e-13 * tc[4]
)* 0.0227226249176305
;} else {
// species 7: O2
result += y[7] * (
+3.28253784e+00
+1.48308754e-03 * tc[1]
-7.57966669e-07 * tc[2]
+2.09470555e-10 * tc[3]
-2.16717794e-14 * tc[4]
)* 0.0312519532470779
;
// species 8: N2
result += y[8] * (
+2.92664000e+00
+1.48797700e-03 * tc[1]
-5.68476100e-07 * tc[2]
+1.00970400e-10 * tc[3]
-6.75335100e-15 * tc[4]
)* 0.0356964374955379
;
// species 10: H2O
result += y[10] * (
+2.67214561e+00
+3.05629289e-03 * tc[1]
-8.73026011e-07 * tc[2]
+1.20099639e-10 * tc[3]
-6.39161787e-15 * tc[4]
)* 0.0555092978073827
;
// species 11: CO2
result += y[11] * (
+3.85746029e+00
+4.41437026e-03 * tc[1]
-2.21481404e-06 * tc[2]
+5.23490188e-10 * tc[3]
-4.72084164e-14 * tc[4]
)* 0.0227226249176305
;}

// species with midpoint at T=1800 kelvin
if (tT < 1800) {
// species 9: NC16H34
result += y[9] * (
-3.64057977e+00
+1.96191732e-01 * tc[1]
-1.26758759e-04 * tc[2]
+4.13645475e-08 * tc[3]
-5.40570458e-12 * tc[4]
)* 0.0044160248710521
;} else {
// species 9: NC16H34
result += y[9] * (
+4.98210237e+01
+7.73881688e-02 * tc[1]
-2.77557896e-05 * tc[2]
+4.69678103e-09 * tc[3]
-3.12959243e-13 * tc[4]
)* 0.0044160248710521
;}

cpbs = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 35)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVBL(const amrex::Real T, const amrex::Real x[], amrex::Real& cvbl)
{
amrex::Real result = 0; 
amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
amrex::Real cvor[12]; //  temporary storage
cv_R(cvor, tc);

// perform dot product
for (int id = 0; id < 12; ++id) {
result += x[id]*cvor[id];
}

cvbl = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 36)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVBS(const amrex::Real T, const amrex::Real y[],  amrex::Real& cvbs)
{
amrex::Real result = 0.0; 
const amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
// compute Cv/R at the given temperature

// species with midpoint at T=1390 kelvin
if (tT < 1390) {
// species 0: NC5H12
result += y[0] * (
-1.73676612e+00
+6.07200973e-02 * tc[1]
-3.57592761e-05 * tc[2]
+1.04907042e-08 * tc[3]
-1.21487315e-12 * tc[4]
)* 0.0138598217626921
;
// species 1: NC6H14
result += y[1] * (
-1.96960618e+00
+7.29085608e-02 * tc[1]
-4.38853919e-05 * tc[2]
+1.32312807e-08 * tc[3]
-1.58437423e-12 * tc[4]
)* 0.0116038896238019
;} else {
// species 0: NC5H12
result += y[0] * (
+1.47257603e+01
+2.61086045e-02 * tc[1]
-8.90970996e-06 * tc[2]
+1.38102248e-09 * tc[3]
-8.00296536e-14 * tc[4]
)* 0.0138598217626921
;
// species 1: NC6H14
result += y[1] * (
+1.79634117e+01
+3.04480204e-02 * tc[1]
-1.03794829e-05 * tc[2]
+1.60775457e-09 * tc[3]
-9.31269728e-14 * tc[4]
)* 0.0116038896238019
;}

// species with midpoint at T=1391 kelvin
if (tT < 1391) {
// species 2: NC8H18
result += y[2] * (
-2.54218406e+00
+9.78112063e-02 * tc[1]
-6.09318358e-05 * tc[2]
+1.92005591e-08 * tc[3]
-2.42996250e-12 * tc[4]
)* 0.0087541144337839
;
// species 4: NC10H22
result += y[4] * (
-3.08416969e+00
+1.22535012e-01 * tc[1]
-7.76815739e-05 * tc[2]
+2.49834877e-08 * tc[3]
-3.23548038e-12 * tc[4]
)* 0.0070280983371519
;
// species 5: NC7H16
result += y[5] * (
-2.26836187e+00
+8.54355820e-02 * tc[1]
-5.25346786e-05 * tc[2]
+1.62945721e-08 * tc[3]
-2.02394925e-12 * tc[4]
)* 0.0099795419390250
;
// species 6: NC12H26
result += y[6] * (
-3.62181594e+00
+1.47237711e-01 * tc[1]
-9.43970271e-05 * tc[2]
+3.07441268e-08 * tc[3]
-4.03602230e-12 * tc[4]
)* 0.0058706117177410
;} else {
// species 2: NC8H18
result += y[2] * (
+2.44710194e+01
+3.90887037e-02 * tc[1]
-1.33038777e-05 * tc[2]
+2.05867527e-09 * tc[3]
-1.19167174e-13 * tc[4]
)* 0.0087541144337839
;
// species 4: NC10H22
result += y[4] * (
+3.09882239e+01
+4.77244922e-02 * tc[1]
-1.62276391e-05 * tc[2]
+2.50963259e-09 * tc[3]
-1.45215772e-13 * tc[4]
)* 0.0070280983371519
;
// species 5: NC7H16
result += y[5] * (
+2.12148969e+01
+3.47675750e-02 * tc[1]
-1.18407129e-05 * tc[2]
+1.83298478e-09 * tc[3]
-1.06130266e-13 * tc[4]
)* 0.0099795419390250
;
// species 6: NC12H26
result += y[6] * (
+3.75095037e+01
+5.63550048e-02 * tc[1]
-1.91493200e-05 * tc[2]
+2.96024862e-09 * tc[3]
-1.71244150e-13 * tc[4]
)* 0.0058706117177410
;}

// species with midpoint at T=1396 kelvin
if (tT < 1396) {
// species 3: IC8H18
result += y[3] * (
-5.20868893e+00
+1.11440581e-01 * tc[1]
-7.91346582e-05 * tc[2]
+2.92406242e-08 * tc[3]
-4.43743191e-12 * tc[4]
)* 0.0087541144337839
;} else {
// species 3: IC8H18
result += y[3] * (
+2.61373590e+01
+3.79004890e-02 * tc[1]
-1.29437358e-05 * tc[2]
+2.00760372e-09 * tc[3]
-1.16400580e-13 * tc[4]
)* 0.0087541144337839
;}

// species with midpoint at T=1000 kelvin
if (tT < 1000) {
// species 7: O2
result += y[7] * (
+2.78245636e+00
-2.99673416e-03 * tc[1]
+9.84730201e-06 * tc[2]
-9.68129509e-09 * tc[3]
+3.24372837e-12 * tc[4]
)* 0.0312519532470779
;
// species 8: N2
result += y[8] * (
+2.29867700e+00
+1.40824000e-03 * tc[1]
-3.96322200e-06 * tc[2]
+5.64151500e-09 * tc[3]
-2.44485500e-12 * tc[4]
)* 0.0356964374955379
;
// species 10: H2O
result += y[10] * (
+2.38684249e+00
+3.47498246e-03 * tc[1]
-6.35469633e-06 * tc[2]
+6.96858127e-09 * tc[3]
-2.50658847e-12 * tc[4]
)* 0.0555092978073827
;
// species 11: CO2
result += y[11] * (
+1.35677352e+00
+8.98459677e-03 * tc[1]
-7.12356269e-06 * tc[2]
+2.45919022e-09 * tc[3]
-1.43699548e-13 * tc[4]
)* 0.0227226249176305
;} else {
// species 7: O2
result += y[7] * (
+2.28253784e+00
+1.48308754e-03 * tc[1]
-7.57966669e-07 * tc[2]
+2.09470555e-10 * tc[3]
-2.16717794e-14 * tc[4]
)* 0.0312519532470779
;
// species 8: N2
result += y[8] * (
+1.92664000e+00
+1.48797700e-03 * tc[1]
-5.68476100e-07 * tc[2]
+1.00970400e-10 * tc[3]
-6.75335100e-15 * tc[4]
)* 0.0356964374955379
;
// species 10: H2O
result += y[10] * (
+1.67214561e+00
+3.05629289e-03 * tc[1]
-8.73026011e-07 * tc[2]
+1.20099639e-10 * tc[3]
-6.39161787e-15 * tc[4]
)* 0.0555092978073827
;
// species 11: CO2
result += y[11] * (
+2.85746029e+00
+4.41437026e-03 * tc[1]
-2.21481404e-06 * tc[2]
+5.23490188e-10 * tc[3]
-4.72084164e-14 * tc[4]
)* 0.0227226249176305
;}

// species with midpoint at T=1800 kelvin
if (tT < 1800) {
// species 9: NC16H34
result += y[9] * (
-4.64057977e+00
+1.96191732e-01 * tc[1]
-1.26758759e-04 * tc[2]
+4.13645475e-08 * tc[3]
-5.40570458e-12 * tc[4]
)* 0.0044160248710521
;} else {
// species 9: NC16H34
result += y[9] * (
+4.88210237e+01
+7.73881688e-02 * tc[1]
-2.77557896e-05 * tc[2]
+4.69678103e-09 * tc[3]
-3.12959243e-13 * tc[4]
)* 0.0044160248710521
;}

cvbs = result * 8.31446261815324e+07;
}

// Returns the mean enthalpy of the mixture in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHBML(const amrex::Real T, const amrex::Real x[], amrex::Real& hbml)
{
amrex::Real result = 0; 
amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
amrex::Real hml[12]; //  temporary storage
amrex::Real RT = 8.31446261815324e+07*tT; // R*T
speciesEnthalpy(hml, tc);

// perform dot product
for (int id = 0; id < 12; ++id) {
result += x[id]*hml[id];
}

hbml = result * RT;
}

// Returns mean enthalpy of mixture in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHBMS(const amrex::Real T, const amrex::Real y[],  amrex::Real& hbms)
{
amrex::Real result = 0.0;
const amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
const amrex::Real invT = 1.0 / T;

// species with midpoint at T=1390 kelvin
if (tT < 1390) {
// species 0: NC5H12
result += y[0] * (
-7.36766122e-01
+3.03600486e-02 * tc[1]
-1.19197587e-05 * tc[2]
+2.62267605e-09 * tc[3]
-2.42974630e-13 * tc[4]
-1.98934934e+04 * invT
)* 0.0138598217626921
;
// species 1: NC6H14
result += y[1] * (
-9.69606184e-01
+3.64542804e-02 * tc[1]
-1.46284640e-05 * tc[2]
+3.30782017e-09 * tc[3]
-3.16874846e-13 * tc[4]
-2.27803862e+04 * invT
)* 0.0116038896238019
;} else {
// species 0: NC5H12
result += y[0] * (
+1.57257603e+01
+1.30543023e-02 * tc[1]
-2.96990332e-06 * tc[2]
+3.45255620e-10 * tc[3]
-1.60059307e-14 * tc[4]
-2.60519543e+04 * invT
)* 0.0138598217626921
;
// species 1: NC6H14
result += y[1] * (
+1.89634117e+01
+1.52240102e-02 * tc[1]
-3.45982763e-06 * tc[2]
+4.01938643e-10 * tc[3]
-1.86253946e-14 * tc[4]
-3.01628739e+04 * invT
)* 0.0116038896238019
;}

// species with midpoint at T=1391 kelvin
if (tT < 1391) {
// species 2: NC8H18
result += y[2] * (
-1.54218406e+00
+4.89056031e-02 * tc[1]
-2.03106119e-05 * tc[2]
+4.80013978e-09 * tc[3]
-4.85992500e-13 * tc[4]
-2.85395641e+04 * invT
)* 0.0087541144337839
;
// species 4: NC10H22
result += y[4] * (
-2.08416969e+00
+6.12675060e-02 * tc[1]
-2.58938580e-05 * tc[2]
+6.24587193e-09 * tc[3]
-6.47096076e-13 * tc[4]
-3.43021863e+04 * invT
)* 0.0070280983371519
;
// species 5: NC7H16
result += y[5] * (
-1.26836187e+00
+4.27177910e-02 * tc[1]
-1.75115595e-05 * tc[2]
+4.07364302e-09 * tc[3]
-4.04789850e-13 * tc[4]
-2.56586565e+04 * invT
)* 0.0099795419390250
;
// species 6: NC12H26
result += y[6] * (
-2.62181594e+00
+7.36188555e-02 * tc[1]
-3.14656757e-05 * tc[2]
+7.68603170e-09 * tc[3]
-8.07204460e-13 * tc[4]
-4.00654253e+04 * invT
)* 0.0058706117177410
;} else {
// species 2: NC8H18
result += y[2] * (
+2.54710194e+01
+1.95443518e-02 * tc[1]
-4.43462590e-06 * tc[2]
+5.14668817e-10 * tc[3]
-2.38334348e-14 * tc[4]
-3.83962755e+04 * invT
)* 0.0087541144337839
;
// species 4: NC10H22
result += y[4] * (
+3.19882239e+01
+2.38622461e-02 * tc[1]
-5.40921303e-06 * tc[2]
+6.27408147e-10 * tc[3]
-2.90431544e-14 * tc[4]
-4.66392840e+04 * invT
)* 0.0070280983371519
;
// species 5: NC7H16
result += y[5] * (
+2.22148969e+01
+1.73837875e-02 * tc[1]
-3.94690430e-06 * tc[2]
+4.58246195e-10 * tc[3]
-2.12260532e-14 * tc[4]
-3.42760081e+04 * invT
)* 0.0099795419390250
;
// species 6: NC12H26
result += y[6] * (
+3.85095037e+01
+2.81775024e-02 * tc[1]
-6.38310667e-06 * tc[2]
+7.40062155e-10 * tc[3]
-3.42488300e-14 * tc[4]
-5.48843465e+04 * invT
)* 0.0058706117177410
;}

// species with midpoint at T=1396 kelvin
if (tT < 1396) {
// species 3: IC8H18
result += y[3] * (
-4.20868893e+00
+5.57202905e-02 * tc[1]
-2.63782194e-05 * tc[2]
+7.31015605e-09 * tc[3]
-8.87486382e-13 * tc[4]
-2.99446875e+04 * invT
)* 0.0087541144337839
;} else {
// species 3: IC8H18
result += y[3] * (
+2.71373590e+01
+1.89502445e-02 * tc[1]
-4.31457860e-06 * tc[2]
+5.01900930e-10 * tc[3]
-2.32801160e-14 * tc[4]
-4.07958177e+04 * invT
)* 0.0087541144337839
;}

// species with midpoint at T=1000 kelvin
if (tT < 1000) {
// species 7: O2
result += y[7] * (
+3.78245636e+00
-1.49836708e-03 * tc[1]
+3.28243400e-06 * tc[2]
-2.42032377e-09 * tc[3]
+6.48745674e-13 * tc[4]
-1.06394356e+03 * invT
)* 0.0312519532470779
;
// species 8: N2
result += y[8] * (
+3.29867700e+00
+7.04120000e-04 * tc[1]
-1.32107400e-06 * tc[2]
+1.41037875e-09 * tc[3]
-4.88971000e-13 * tc[4]
-1.02090000e+03 * invT
)* 0.0356964374955379
;
// species 10: H2O
result += y[10] * (
+3.38684249e+00
+1.73749123e-03 * tc[1]
-2.11823211e-06 * tc[2]
+1.74214532e-09 * tc[3]
-5.01317694e-13 * tc[4]
-3.02081133e+04 * invT
)* 0.0555092978073827
;
// species 11: CO2
result += y[11] * (
+2.35677352e+00
+4.49229839e-03 * tc[1]
-2.37452090e-06 * tc[2]
+6.14797555e-10 * tc[3]
-2.87399096e-14 * tc[4]
-4.83719697e+04 * invT
)* 0.0227226249176305
;} else {
// species 7: O2
result += y[7] * (
+3.28253784e+00
+7.41543770e-04 * tc[1]
-2.52655556e-07 * tc[2]
+5.23676387e-11 * tc[3]
-4.33435588e-15 * tc[4]
-1.08845772e+03 * invT
)* 0.0312519532470779
;
// species 8: N2
result += y[8] * (
+2.92664000e+00
+7.43988500e-04 * tc[1]
-1.89492033e-07 * tc[2]
+2.52426000e-11 * tc[3]
-1.35067020e-15 * tc[4]
-9.22797700e+02 * invT
)* 0.0356964374955379
;
// species 10: H2O
result += y[10] * (
+2.67214561e+00
+1.52814644e-03 * tc[1]
-2.91008670e-07 * tc[2]
+3.00249098e-11 * tc[3]
-1.27832357e-15 * tc[4]
-2.98992090e+04 * invT
)* 0.0555092978073827
;
// species 11: CO2
result += y[11] * (
+3.85746029e+00
+2.20718513e-03 * tc[1]
-7.38271347e-07 * tc[2]
+1.30872547e-10 * tc[3]
-9.44168328e-15 * tc[4]
-4.87591660e+04 * invT
)* 0.0227226249176305
;}

// species with midpoint at T=1800 kelvin
if (tT < 1800) {
// species 9: NC16H34
result += y[9] * (
-3.64057977e+00
+9.80958660e-02 * tc[1]
-4.22529197e-05 * tc[2]
+1.03411369e-08 * tc[3]
-1.08114092e-12 * tc[4]
-5.16052891e+04 * invT
)* 0.0044160248710521
;} else {
// species 9: NC16H34
result += y[9] * (
+4.98210237e+01
+3.86940844e-02 * tc[1]
-9.25192987e-06 * tc[2]
+1.17419526e-09 * tc[3]
-6.25918486e-14 * tc[4]
-7.08514663e+04 * invT
)* 0.0044160248710521
;}

const amrex::Real RT = 8.31446261815324e+07*tT; // R*T

hbms = result * RT;
}

// get mean internal energy in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUBML(const amrex::Real T, const amrex::Real x[], amrex::Real& ubml)
{
amrex::Real result = 0; 
amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
amrex::Real uml[12]; //  temporary energy array
amrex::Real RT = 8.31446261815324e+07*tT; // R*T
speciesInternalEnergy(uml, tc);

// perform dot product
for (int id = 0; id < 12; ++id) {
result += x[id]*uml[id];
}

ubml = result * RT;
}

// get mean internal energy in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUBMS(const amrex::Real T, const amrex::Real y[], amrex::Real& ubms)
{
amrex::Real result = 0.0;
const amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
const amrex::Real invT = 1.0 / T;

// species with midpoint at T=1390 kelvin
if (tT < 1390) {
// species 0: NC5H12
result += y[0] * (
-1.73676612e+00
+3.03600486e-02 * tc[1]
-1.19197587e-05 * tc[2]
+2.62267605e-09 * tc[3]
-2.42974630e-13 * tc[4]
-1.98934934e+04 * invT
)* 0.0138598217626921
;
// species 1: NC6H14
result += y[1] * (
-1.96960618e+00
+3.64542804e-02 * tc[1]
-1.46284640e-05 * tc[2]
+3.30782017e-09 * tc[3]
-3.16874846e-13 * tc[4]
-2.27803862e+04 * invT
)* 0.0116038896238019
;} else {
// species 0: NC5H12
result += y[0] * (
+1.47257603e+01
+1.30543023e-02 * tc[1]
-2.96990332e-06 * tc[2]
+3.45255620e-10 * tc[3]
-1.60059307e-14 * tc[4]
-2.60519543e+04 * invT
)* 0.0138598217626921
;
// species 1: NC6H14
result += y[1] * (
+1.79634117e+01
+1.52240102e-02 * tc[1]
-3.45982763e-06 * tc[2]
+4.01938643e-10 * tc[3]
-1.86253946e-14 * tc[4]
-3.01628739e+04 * invT
)* 0.0116038896238019
;}

// species with midpoint at T=1391 kelvin
if (tT < 1391) {
// species 2: NC8H18
result += y[2] * (
-2.54218406e+00
+4.89056031e-02 * tc[1]
-2.03106119e-05 * tc[2]
+4.80013978e-09 * tc[3]
-4.85992500e-13 * tc[4]
-2.85395641e+04 * invT
)* 0.0087541144337839
;
// species 4: NC10H22
result += y[4] * (
-3.08416969e+00
+6.12675060e-02 * tc[1]
-2.58938580e-05 * tc[2]
+6.24587193e-09 * tc[3]
-6.47096076e-13 * tc[4]
-3.43021863e+04 * invT
)* 0.0070280983371519
;
// species 5: NC7H16
result += y[5] * (
-2.26836187e+00
+4.27177910e-02 * tc[1]
-1.75115595e-05 * tc[2]
+4.07364302e-09 * tc[3]
-4.04789850e-13 * tc[4]
-2.56586565e+04 * invT
)* 0.0099795419390250
;
// species 6: NC12H26
result += y[6] * (
-3.62181594e+00
+7.36188555e-02 * tc[1]
-3.14656757e-05 * tc[2]
+7.68603170e-09 * tc[3]
-8.07204460e-13 * tc[4]
-4.00654253e+04 * invT
)* 0.0058706117177410
;} else {
// species 2: NC8H18
result += y[2] * (
+2.44710194e+01
+1.95443518e-02 * tc[1]
-4.43462590e-06 * tc[2]
+5.14668817e-10 * tc[3]
-2.38334348e-14 * tc[4]
-3.83962755e+04 * invT
)* 0.0087541144337839
;
// species 4: NC10H22
result += y[4] * (
+3.09882239e+01
+2.38622461e-02 * tc[1]
-5.40921303e-06 * tc[2]
+6.27408147e-10 * tc[3]
-2.90431544e-14 * tc[4]
-4.66392840e+04 * invT
)* 0.0070280983371519
;
// species 5: NC7H16
result += y[5] * (
+2.12148969e+01
+1.73837875e-02 * tc[1]
-3.94690430e-06 * tc[2]
+4.58246195e-10 * tc[3]
-2.12260532e-14 * tc[4]
-3.42760081e+04 * invT
)* 0.0099795419390250
;
// species 6: NC12H26
result += y[6] * (
+3.75095037e+01
+2.81775024e-02 * tc[1]
-6.38310667e-06 * tc[2]
+7.40062155e-10 * tc[3]
-3.42488300e-14 * tc[4]
-5.48843465e+04 * invT
)* 0.0058706117177410
;}

// species with midpoint at T=1396 kelvin
if (tT < 1396) {
// species 3: IC8H18
result += y[3] * (
-5.20868893e+00
+5.57202905e-02 * tc[1]
-2.63782194e-05 * tc[2]
+7.31015605e-09 * tc[3]
-8.87486382e-13 * tc[4]
-2.99446875e+04 * invT
)* 0.0087541144337839
;} else {
// species 3: IC8H18
result += y[3] * (
+2.61373590e+01
+1.89502445e-02 * tc[1]
-4.31457860e-06 * tc[2]
+5.01900930e-10 * tc[3]
-2.32801160e-14 * tc[4]
-4.07958177e+04 * invT
)* 0.0087541144337839
;}

// species with midpoint at T=1000 kelvin
if (tT < 1000) {
// species 7: O2
result += y[7] * (
+2.78245636e+00
-1.49836708e-03 * tc[1]
+3.28243400e-06 * tc[2]
-2.42032377e-09 * tc[3]
+6.48745674e-13 * tc[4]
-1.06394356e+03 * invT
)* 0.0312519532470779
;
// species 8: N2
result += y[8] * (
+2.29867700e+00
+7.04120000e-04 * tc[1]
-1.32107400e-06 * tc[2]
+1.41037875e-09 * tc[3]
-4.88971000e-13 * tc[4]
-1.02090000e+03 * invT
)* 0.0356964374955379
;
// species 10: H2O
result += y[10] * (
+2.38684249e+00
+1.73749123e-03 * tc[1]
-2.11823211e-06 * tc[2]
+1.74214532e-09 * tc[3]
-5.01317694e-13 * tc[4]
-3.02081133e+04 * invT
)* 0.0555092978073827
;
// species 11: CO2
result += y[11] * (
+1.35677352e+00
+4.49229839e-03 * tc[1]
-2.37452090e-06 * tc[2]
+6.14797555e-10 * tc[3]
-2.87399096e-14 * tc[4]
-4.83719697e+04 * invT
)* 0.0227226249176305
;} else {
// species 7: O2
result += y[7] * (
+2.28253784e+00
+7.41543770e-04 * tc[1]
-2.52655556e-07 * tc[2]
+5.23676387e-11 * tc[3]
-4.33435588e-15 * tc[4]
-1.08845772e+03 * invT
)* 0.0312519532470779
;
// species 8: N2
result += y[8] * (
+1.92664000e+00
+7.43988500e-04 * tc[1]
-1.89492033e-07 * tc[2]
+2.52426000e-11 * tc[3]
-1.35067020e-15 * tc[4]
-9.22797700e+02 * invT
)* 0.0356964374955379
;
// species 10: H2O
result += y[10] * (
+1.67214561e+00
+1.52814644e-03 * tc[1]
-2.91008670e-07 * tc[2]
+3.00249098e-11 * tc[3]
-1.27832357e-15 * tc[4]
-2.98992090e+04 * invT
)* 0.0555092978073827
;
// species 11: CO2
result += y[11] * (
+2.85746029e+00
+2.20718513e-03 * tc[1]
-7.38271347e-07 * tc[2]
+1.30872547e-10 * tc[3]
-9.44168328e-15 * tc[4]
-4.87591660e+04 * invT
)* 0.0227226249176305
;}

// species with midpoint at T=1800 kelvin
if (tT < 1800) {
// species 9: NC16H34
result += y[9] * (
-4.64057977e+00
+9.80958660e-02 * tc[1]
-4.22529197e-05 * tc[2]
+1.03411369e-08 * tc[3]
-1.08114092e-12 * tc[4]
-5.16052891e+04 * invT
)* 0.0044160248710521
;} else {
// species 9: NC16H34
result += y[9] * (
+4.88210237e+01
+3.86940844e-02 * tc[1]
-9.25192987e-06 * tc[2]
+1.17419526e-09 * tc[3]
-6.25918486e-14 * tc[4]
-7.08514663e+04 * invT
)* 0.0044160248710521
;}

const amrex::Real RT = 8.31446261815324e+07*tT; // R*T

ubms = result * RT;
}

// get mixture entropy in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSBML(const amrex::Real P, const amrex::Real T,const amrex::Real x[], amrex::Real& sbml)
{
amrex::Real result = 0; 
// Log of normalized pressure in cgs units dynes/cm^2 by Patm
amrex::Real logPratio = log ( P / 1013250.0 ); 
amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
amrex::Real sor[12]; //  temporary storage
speciesEntropy(sor, tc);

// Compute Eq 42
for (int id = 0; id < 12; ++id) {
result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
}

sbml = result * 8.31446261815324e+07;
}

// get mixture entropy in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSBMS(const amrex::Real P, const amrex::Real T,const amrex::Real y[], amrex::Real& sbms)
{
amrex::Real result = 0; 
// Log of normalized pressure in cgs units dynes/cm^2 by Patm
amrex::Real logPratio = log ( P / 1013250.0 ); 
amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
amrex::Real sor[12]; //  temporary storage
amrex::Real x[12]; //  need a ytx conversion
amrex::Real YOW = 0; // See Eq 4, 6 in CK Manual

// Compute inverse of mean molecular wt first
for (int i = 0; i < 12; i++)
{
YOW += y[i]* imw(i);
}
// Now compute y to x conversion
x[0] = y[0]/(72.151000*YOW); 
x[1] = y[1]/(86.178000*YOW); 
x[2] = y[2]/(114.232000*YOW); 
x[3] = y[3]/(114.232000*YOW); 
x[4] = y[4]/(142.286000*YOW); 
x[5] = y[5]/(100.205000*YOW); 
x[6] = y[6]/(170.340000*YOW); 
x[7] = y[7]/(31.998000*YOW); 
x[8] = y[8]/(28.014000*YOW); 
x[9] = y[9]/(226.448000*YOW); 
x[10] = y[10]/(18.015000*YOW); 
x[11] = y[11]/(44.009000*YOW); 
speciesEntropy(sor, tc);
// Perform computation in Eq 42 and 43
for (int i = 0; i < 12; i++)
{
result += x[i]*(sor[i]-log((x[i]+1e-100))-logPratio);
}
// Scale by R/W
sbms = result * 8.31446261815324e+07 * YOW;
}

//  get temperature given internal energy in mass units and mass fracs
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void GET_T_GIVEN_EY(const amrex::Real e, const amrex::Real y[], amrex::Real& t, int& ierr)
{
#ifdef CONVERGENCE
const int maxiter = 5000;
const amrex::Real tol = 1.e-12;
#else
const int maxiter = 200;
const amrex::Real tol = 1.e-6;
#endif
amrex::Real tmin = 90;// max lower bound for thermo def
amrex::Real tmax = 4000;// min upper bound for thermo def
amrex::Real e1,emin,emax,cv,t1,dt;
CKUBMS(tmin, y, emin);
CKUBMS(tmax, y, emax);
if (e < emin) {
// Linear Extrapolation below tmin
CKCVBS(tmin, y, cv);
t = tmin - (emin-e)/cv;
ierr = 1;
return;
}
if (e > emax) {
// Linear Extrapolation above tmax
CKCVBS(tmax, y, cv);
t = tmax - (emax-e)/cv;
ierr = 1;
return;
}
t1 = t;
if (t1 < tmin || t1 > tmax) {
t1 = tmin + (tmax-tmin)/(emax-emin)*(e-emin);
}
for (int i = 0; i < maxiter; ++i) {
CKUBMS(t1,y,e1);
CKCVBS(t1,y,cv);
dt = (e - e1) / cv;
if (dt > 100.) { dt = 100.; }
else if (dt < -100.) { dt = -100.; }
else if (fabs(dt) < tol) {break;}
t1 += dt;
}
t = t1;
ierr = 0;
}

//  get temperature given enthalpy in mass units and mass fracs
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void GET_T_GIVEN_HY(const amrex::Real h, const amrex::Real y[], amrex::Real& t, int& ierr)
{
#ifdef CONVERGENCE
const int maxiter = 5000;
const amrex::Real tol = 1.e-12;
#else
const int maxiter = 200;
const amrex::Real tol = 1.e-6;
#endif
amrex::Real tmin = 90;// max lower bound for thermo def
amrex::Real tmax = 4000;// min upper bound for thermo def
amrex::Real h1,hmin,hmax,cp,t1,dt;
CKHBMS(tmin, y, hmin);
CKHBMS(tmax, y, hmax);
if (h < hmin) {
// Linear Extrapolation below tmin
CKCPBS(tmin, y, cp);
t = tmin - (hmin-h)/cp;
ierr = 1;
return;
}
if (h > hmax) {
// Linear Extrapolation above tmax
CKCPBS(tmax, y, cp);
t = tmax - (hmax-h)/cp;
ierr = 1;
return;
}
t1 = t;
if (t1 < tmin || t1 > tmax) {
t1 = tmin + (tmax-tmin)/(hmax-hmin)*(h-hmin);
}
for (int i = 0; i < maxiter; ++i) {
CKHBMS(t1,y,h1);
CKCPBS(t1,y,cp);
dt = (h - h1) / cp;
if (dt > 100.) { dt = 100.; }
else if (dt < -100.) { dt = -100.; }
else if (fabs(dt) < tol) {break;}
t1 += dt;
}
t = t1;
ierr = 0;
}

// Compute P = rhoRT/W(x)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPX(const amrex::Real rho, const amrex::Real T, const amrex::Real x[], amrex::Real& P)
{
amrex::Real XW = 0;//  To hold mean molecular wt
XW += x[0]*72.151000; // NC5H12
XW += x[1]*86.178000; // NC6H14
XW += x[2]*114.232000; // NC8H18
XW += x[3]*114.232000; // IC8H18
XW += x[4]*142.286000; // NC10H22
XW += x[5]*100.205000; // NC7H16
XW += x[6]*170.340000; // NC12H26
XW += x[7]*31.998000; // O2
XW += x[8]*28.014000; // N2
XW += x[9]*226.448000; // NC16H34
XW += x[10]*18.015000; // H2O
XW += x[11]*44.009000; // CO2
P = rho * 8.31446261815324e+07 * T / XW; // P = rho*R*T/W

}

// Compute P = rhoRT/W(y)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPY(const amrex::Real rho, const amrex::Real T, const amrex::Real y[],  amrex::Real& P)
{
amrex::Real YOW = 0;//  for computing mean MW

for (int i = 0; i < 12; i++)
{
YOW += y[i]*imw(i);
}
P = rho * 8.31446261815324e+07 * T * YOW; // P = rho*R*T/W

}

// Compute P = rhoRT/W(c)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPC(const amrex::Real rho, const amrex::Real T, const amrex::Real c[],  amrex::Real& P)
{
// See Eq 5 in CK Manual
amrex::Real W = 0;
amrex::Real sumC = 0;
W += c[0]*72.151000; // NC5H12
W += c[1]*86.178000; // NC6H14
W += c[2]*114.232000; // NC8H18
W += c[3]*114.232000; // IC8H18
W += c[4]*142.286000; // NC10H22
W += c[5]*100.205000; // NC7H16
W += c[6]*170.340000; // NC12H26
W += c[7]*31.998000; // O2
W += c[8]*28.014000; // N2
W += c[9]*226.448000; // NC16H34
W += c[10]*18.015000; // H2O
W += c[11]*44.009000; // CO2

for (int id = 0; id < 12; ++id) {
sumC += c[id];
}
P = rho * 8.31446261815324e+07 * T * sumC / W; // P = rho*R*T/W

}

// Compute rho = PW(x)/RT
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOX(const amrex::Real P, const amrex::Real T,const amrex::Real x[], amrex::Real& rho)
{
amrex::Real XW = 0;//  To hold mean molecular wt
XW += x[0]*72.151000; // NC5H12
XW += x[1]*86.178000; // NC6H14
XW += x[2]*114.232000; // NC8H18
XW += x[3]*114.232000; // IC8H18
XW += x[4]*142.286000; // NC10H22
XW += x[5]*100.205000; // NC7H16
XW += x[6]*170.340000; // NC12H26
XW += x[7]*31.998000; // O2
XW += x[8]*28.014000; // N2
XW += x[9]*226.448000; // NC16H34
XW += x[10]*18.015000; // H2O
XW += x[11]*44.009000; // CO2
rho = P * XW / (8.31446261815324e+07 * T); // rho = P*W/(R*T)

}

// Compute rho = P*W(y)/RT
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOY(const amrex::Real P, const amrex::Real T,const amrex::Real y[], amrex::Real& rho)
{
amrex::Real YOW = 0;

for (int i = 0; i < 12; i++)
{
YOW += y[i]*imw(i);
}

rho = P / (8.31446261815324e+07 * T * YOW);// rho = P*W/(R*T)
}

// Compute rho = P*W(c)/(R*T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOC(const amrex::Real P, const amrex::Real T,const amrex::Real c[], amrex::Real& rho)
{
// See Eq 5 in CK Manual
amrex::Real W = 0;
amrex::Real sumC = 0;
W += c[0]*72.151000; // NC5H12
W += c[1]*86.178000; // NC6H14
W += c[2]*114.232000; // NC8H18
W += c[3]*114.232000; // IC8H18
W += c[4]*142.286000; // NC10H22
W += c[5]*100.205000; // NC7H16
W += c[6]*170.340000; // NC12H26
W += c[7]*31.998000; // O2
W += c[8]*28.014000; // N2
W += c[9]*226.448000; // NC16H34
W += c[10]*18.015000; // H2O
W += c[11]*44.009000; // CO2

for (int id = 0; id < 12; ++id) {
sumC += c[id];
}
rho = P * W / (sumC * T * 8.31446261815324e+07); // rho = PW/(R*T)

}

// get molecular weight for all species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWT( amrex::Real wt[])
{
get_mw(wt);
}

// given y[species]: mass fractions
// s mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWY(const amrex::Real y[], amrex::Real& wtm)
{
amrex::Real YOW = 0;

for (int i = 0; i < 12; i++)
{
YOW += y[i]*imw(i);
}

wtm = 1.0 / YOW;
}

// given x[species]: mole fractions
// returns mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWX(const amrex::Real x[],  amrex::Real& wtm)
{
amrex::Real XW = 0;//  see Eq 4 in CK Manual
XW += x[0]*72.151000; // NC5H12
XW += x[1]*86.178000; // NC6H14
XW += x[2]*114.232000; // NC8H18
XW += x[3]*114.232000; // IC8H18
XW += x[4]*142.286000; // NC10H22
XW += x[5]*100.205000; // NC7H16
XW += x[6]*170.340000; // NC12H26
XW += x[7]*31.998000; // O2
XW += x[8]*28.014000; // N2
XW += x[9]*226.448000; // NC16H34
XW += x[10]*18.015000; // H2O
XW += x[11]*44.009000; // CO2
wtm = XW;

}

// given c[species]: molar concentration
// returns mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWC(const amrex::Real c[],  amrex::Real& wtm)
{
// See Eq 5 in CK Manual
amrex::Real W = 0;
amrex::Real sumC = 0;
W += c[0]*72.151000; // NC5H12
W += c[1]*86.178000; // NC6H14
W += c[2]*114.232000; // NC8H18
W += c[3]*114.232000; // IC8H18
W += c[4]*142.286000; // NC10H22
W += c[5]*100.205000; // NC7H16
W += c[6]*170.340000; // NC12H26
W += c[7]*31.998000; // O2
W += c[8]*28.014000; // N2
W += c[9]*226.448000; // NC16H34
W += c[10]*18.015000; // H2O
W += c[11]*44.009000; // CO2

for (int id = 0; id < 12; ++id) {
sumC += c[id];
}
//  CK provides no guard against divison by zero
wtm = W/sumC;

}

// get Cp/R as a function of T 
// for all species (Eq 19)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPOR(const amrex::Real T, amrex::Real cpor[])
{
amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
cp_R(cpor, tc);
}

// get H/RT as a function of T 
// for all species (Eq 20)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHORT(const amrex::Real T, amrex::Real hort[])
{
amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
speciesEnthalpy(hort, tc);
}

// get S/R as a function of T 
// for all species (Eq 21)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSOR(const amrex::Real T, amrex::Real sor[])
{
amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
speciesEntropy(sor, tc);
}

// convert y[species] (mass fracs) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTX(const amrex::Real y[], amrex::Real x[])
{
amrex::Real YOW = 0;

for (int i = 0; i < 12; i++)
{
YOW += y[i]*imw(i);
}

amrex::Real YOWINV = 1.0/YOW;

for (int i = 0; i < 12; i++)
{
x[i] = y[i]*imw(i)*YOWINV;
}
}

// convert y[species] (mass fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTCP(const amrex::Real P, const amrex::Real T,const amrex::Real y[], amrex::Real c[])
{
amrex::Real YOW = 0;
amrex::Real PWORT;

// Compute inverse of mean molecular wt first
for (int i = 0; i < 12; i++)
{
c[i] = y[i]*imw(i);
}
for (int i = 0; i < 12; i++)
{
YOW += c[i];
}

// PW/RT (see Eq. 7)
PWORT = P/(YOW * 8.31446261815324e+07 * T); 
// Now compute conversion

for (int i = 0; i < 12; i++)
{
c[i] = PWORT * y[i] * imw(i);
}
}

// convert y[species] (mass fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTCR(const amrex::Real rho, amrex::Real /*T*/, const amrex::Real y[],  amrex::Real c[])
{

for (int i = 0; i < 12; i++)
{
c[i] = rho  * y[i] * imw(i);
}
}

// convert x[species] (mole fracs) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTY(const amrex::Real x[], amrex::Real y[])
{
amrex::Real XW = 0; // See Eq 4, 9 in CK Manual
// Compute mean molecular wt first
XW += x[0]*72.151000; // NC5H12
XW += x[1]*86.178000; // NC6H14
XW += x[2]*114.232000; // NC8H18
XW += x[3]*114.232000; // IC8H18
XW += x[4]*142.286000; // NC10H22
XW += x[5]*100.205000; // NC7H16
XW += x[6]*170.340000; // NC12H26
XW += x[7]*31.998000; // O2
XW += x[8]*28.014000; // N2
XW += x[9]*226.448000; // NC16H34
XW += x[10]*18.015000; // H2O
XW += x[11]*44.009000; // CO2
// Now compute conversion
amrex::Real XWinv = 1.0/XW;
y[0] = x[0]*72.151000*XWinv; 
y[1] = x[1]*86.178000*XWinv; 
y[2] = x[2]*114.232000*XWinv; 
y[3] = x[3]*114.232000*XWinv; 
y[4] = x[4]*142.286000*XWinv; 
y[5] = x[5]*100.205000*XWinv; 
y[6] = x[6]*170.340000*XWinv; 
y[7] = x[7]*31.998000*XWinv; 
y[8] = x[8]*28.014000*XWinv; 
y[9] = x[9]*226.448000*XWinv; 
y[10] = x[10]*18.015000*XWinv; 
y[11] = x[11]*44.009000*XWinv; 

}

// convert x[species] (mole fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTCP(const amrex::Real P, const amrex::Real T,const amrex::Real x[], amrex::Real c[])
{
amrex::Real PORT = P/(8.31446261815324e+07 * T); // P/RT

// Compute conversion, see Eq 10
for (int id = 0; id < 12; ++id) {
c[id] = x[id]*PORT;
}

}

// convert x[species] (mole fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTCR(const amrex::Real rho, const amrex::Real /*T*/, const amrex::Real x[], amrex::Real c[])
{
amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
amrex::Real ROW; 
// Compute mean molecular wt first
XW += x[0]*72.151000; // NC5H12
XW += x[1]*86.178000; // NC6H14
XW += x[2]*114.232000; // NC8H18
XW += x[3]*114.232000; // IC8H18
XW += x[4]*142.286000; // NC10H22
XW += x[5]*100.205000; // NC7H16
XW += x[6]*170.340000; // NC12H26
XW += x[7]*31.998000; // O2
XW += x[8]*28.014000; // N2
XW += x[9]*226.448000; // NC16H34
XW += x[10]*18.015000; // H2O
XW += x[11]*44.009000; // CO2
ROW = rho / XW;

// Compute conversion, see Eq 11
for (int id = 0; id < 12; ++id) {
c[id] = x[id]*ROW;
}

}

// convert c[species] (molar conc) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCTX(const amrex::Real c[], amrex::Real x[])
{
amrex::Real sumC = 0; 

// compute sum of c 
for (int id = 0; id < 12; ++id) {
sumC += c[id];
}

//  See Eq 13 
amrex::Real sumCinv = 1.0/sumC;
for (int id = 0; id < 12; ++id) {
x[id] = c[id]*sumCinv;
}

}

// convert c[species] (molar conc) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCTY(const amrex::Real c[], amrex::Real y[])
{
amrex::Real CW = 0; // See Eq 12 in CK Manual
// compute denominator in eq 12 first
CW += c[0]*72.151000; // NC5H12
CW += c[1]*86.178000; // NC6H14
CW += c[2]*114.232000; // NC8H18
CW += c[3]*114.232000; // IC8H18
CW += c[4]*142.286000; // NC10H22
CW += c[5]*100.205000; // NC7H16
CW += c[6]*170.340000; // NC12H26
CW += c[7]*31.998000; // O2
CW += c[8]*28.014000; // N2
CW += c[9]*226.448000; // NC16H34
CW += c[10]*18.015000; // H2O
CW += c[11]*44.009000; // CO2
// Now compute conversion
amrex::Real CWinv = 1.0/CW;
y[0] = c[0]*72.151000*CWinv; 
y[1] = c[1]*86.178000*CWinv; 
y[2] = c[2]*114.232000*CWinv; 
y[3] = c[3]*114.232000*CWinv; 
y[4] = c[4]*142.286000*CWinv; 
y[5] = c[5]*100.205000*CWinv; 
y[6] = c[6]*170.340000*CWinv; 
y[7] = c[7]*31.998000*CWinv; 
y[8] = c[8]*28.014000*CWinv; 
y[9] = c[9]*226.448000*CWinv; 
y[10] = c[10]*18.015000*CWinv; 
y[11] = c[11]*44.009000*CWinv; 

}

// get specific heat at constant volume as a function 
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVML(const amrex::Real T, amrex::Real cvml[])
{
amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
cv_R(cvml, tc);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
cvml[id] *= 8.31446261815324e+07;
}
}

// get specific heat at constant pressure as a 
// function of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPML(const amrex::Real T, amrex::Real cpml[])
{
amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
cp_R(cpml, tc);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
cpml[id] *= 8.31446261815324e+07;
}
}

// get internal energy as a function 
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUML(const amrex::Real T, amrex::Real uml[])
{
amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
amrex::Real RT = 8.31446261815324e+07*tT; // R*T
speciesInternalEnergy(uml, tc);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
uml[id] *= RT;
}
}

// get enthalpy as a function 
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHML(const amrex::Real T, amrex::Real hml[])
{
amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
amrex::Real RT = 8.31446261815324e+07*tT; // R*T
speciesEnthalpy(hml, tc);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
hml[id] *= RT;
}
}

// Returns the standard-state entropies in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSML(const amrex::Real T, amrex::Real sml[])
{
amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
speciesEntropy(sml, tc);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
sml[id] *= 8.31446261815324e+07;
}
}

// Returns the specific heats at constant volume
// in mass units (Eq. 29)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVMS(const amrex::Real T, amrex::Real cvms[])
{
amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
cv_R(cvms, tc);
// multiply by R/molecularweight
cvms[0] *= 1.152369699401705e+06; // NC5H12
cvms[1] *= 9.648010650227715e+05; // NC6H14
cvms[2] *= 7.278575721473177e+05; // NC8H18
cvms[3] *= 7.278575721473177e+05; // IC8H18
cvms[4] *= 5.843486090095469e+05; // NC10H22
cvms[5] *= 8.297452839831585e+05; // NC7H16
cvms[6] *= 4.881098167284983e+05; // NC12H26
cvms[7] *= 2.598431970171023e+06; // O2
cvms[8] *= 2.967966951578939e+06; // N2
cvms[9] *= 3.671687371119745e+05; // NC16H34
cvms[10] *= 4.615299815794193e+06; // H2O
cvms[11] *= 1.889264154639560e+06; // CO2
}

// Returns the specific heats at constant pressure
// in mass units (Eq. 26)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPMS(const amrex::Real T, amrex::Real cpms[])
{
amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
cp_R(cpms, tc);
// multiply by R/molecularweight
cpms[0] *= 1.152369699401705e+06; // NC5H12
cpms[1] *= 9.648010650227715e+05; // NC6H14
cpms[2] *= 7.278575721473177e+05; // NC8H18
cpms[3] *= 7.278575721473177e+05; // IC8H18
cpms[4] *= 5.843486090095469e+05; // NC10H22
cpms[5] *= 8.297452839831585e+05; // NC7H16
cpms[6] *= 4.881098167284983e+05; // NC12H26
cpms[7] *= 2.598431970171023e+06; // O2
cpms[8] *= 2.967966951578939e+06; // N2
cpms[9] *= 3.671687371119745e+05; // NC16H34
cpms[10] *= 4.615299815794193e+06; // H2O
cpms[11] *= 1.889264154639560e+06; // CO2
}

// Returns internal energy in mass units (Eq 30.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUMS(const amrex::Real T, amrex::Real ums[])
{
amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
amrex::Real RT = 8.31446261815324e+07*tT; // R*T

speciesInternalEnergy(ums, tc);

for (int i = 0; i < 12; i++)
{
ums[i] *= RT*imw(i);
}
}

// Returns enthalpy in mass units (Eq 27.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHMS(const amrex::Real T, amrex::Real hms[])
{
amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
amrex::Real RT = 8.31446261815324e+07*tT; // R*T

speciesEnthalpy(hms, tc);

for (int i = 0; i < 12; i++)
{
hms[i] *= RT*imw(i);
}
}

// Returns the entropies in mass units (Eq 28.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSMS(const amrex::Real T, amrex::Real sms[])
{
amrex::Real tT = T; // temporary temperature
const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; // temperature cache
speciesEntropy(sms, tc);
// multiply by R/molecularweight
sms[0] *= 1.152369699401705e+06; // NC5H12
sms[1] *= 9.648010650227715e+05; // NC6H14
sms[2] *= 7.278575721473177e+05; // NC8H18
sms[3] *= 7.278575721473177e+05; // IC8H18
sms[4] *= 5.843486090095469e+05; // NC10H22
sms[5] *= 8.297452839831585e+05; // NC7H16
sms[6] *= 4.881098167284983e+05; // NC12H26
sms[7] *= 2.598431970171023e+06; // O2
sms[8] *= 2.967966951578939e+06; // N2
sms[9] *= 3.671687371119745e+05; // NC16H34
sms[10] *= 4.615299815794193e+06; // H2O
sms[11] *= 1.889264154639560e+06; // CO2
}

// GPU version of productionRate: no more use of thermo namespace vectors 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_qfqr(amrex::Real * /*qf*/, amrex::Real * /*qr*/, const amrex::Real * /*sc*/, const amrex::Real * /*sc_qss*/, const amrex::Real * /*tc*/, const amrex::Real /*invT*/)
{

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void  productionRate(amrex::Real * wdot, const amrex::Real * /*sc*/, const amrex::Real /*T*/)
{


for (int i = 0; i < 12; ++i) {
wdot[i] = 0.0;
}


}


// compute the production rate for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWC(const amrex::Real T, amrex::Real C[], amrex::Real wdot[])
{

// convert to SI
for (int id = 0; id < 12; ++id) {
C[id] *= 1.0e6;
}

// convert to chemkin units
productionRate(wdot, C, T);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
C[id] *= 1.0e-6;
wdot[id] *= 1.0e-6;
}
}

// Returns the molar production rate of species
// Given P, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYP(const amrex::Real P, const amrex::Real T,const amrex::Real y[], amrex::Real wdot[])
{
amrex::Real c[12]; // temporary storage
amrex::Real YOW = 0; 
amrex::Real PWORT; 

// Compute inverse of mean molecular wt first
for (int i = 0; i < 12; i++)
{
YOW += y[i]*imw(i);
}
// PW/RT (see Eq. 7)
PWORT = P/(YOW * 8.31446261815324e+07 * T); 
// multiply by 1e6 so c goes to SI
PWORT *= 1e6; 
// Now compute conversion (and go to SI)
for (int i = 0; i < 12; i++)
{
c[i] = PWORT * y[i]*imw(i);
}

// convert to chemkin units
productionRate(wdot, c, T);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
wdot[id] *= 1.0e-6;
}
}

// Returns the molar production rate of species
// Given P, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXP(const amrex::Real P, const amrex::Real T,const amrex::Real x[], amrex::Real wdot[])
{
amrex::Real c[12]; // temporary storage
amrex::Real PORT = 1e6 * P/(8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

// Compute conversion, see Eq 10
for (int id = 0; id < 12; ++id) {
c[id] = x[id]*PORT;
}

// convert to chemkin units
productionRate(wdot, c, T);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
wdot[id] *= 1.0e-6;
}
}

// Returns the molar production rate of species
// Given rho, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYR(const amrex::Real rho, const amrex::Real T, const amrex::Real y[],  amrex::Real wdot[])
{
amrex::Real c[12]; // temporary storage

// See Eq 8 with an extra 1e6 so c goes to SI
for (int i = 0; i < 12; i++)
{
c[i] = 1e6 * rho * y[i]*imw(i);
}

// call productionRate
productionRate(wdot, c, T);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
wdot[id] *= 1.0e-6;
}
}

// Returns the molar production rate of species
// Given rho, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXR(const amrex::Real rho, const amrex::Real T, const amrex::Real x[],  amrex::Real wdot[])
{
amrex::Real c[12]; // temporary storage
amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
amrex::Real ROW; 
// Compute mean molecular wt first
XW += x[0]*72.151000; // NC5H12
XW += x[1]*86.178000; // NC6H14
XW += x[2]*114.232000; // NC8H18
XW += x[3]*114.232000; // IC8H18
XW += x[4]*142.286000; // NC10H22
XW += x[5]*100.205000; // NC7H16
XW += x[6]*170.340000; // NC12H26
XW += x[7]*31.998000; // O2
XW += x[8]*28.014000; // N2
XW += x[9]*226.448000; // NC16H34
XW += x[10]*18.015000; // H2O
XW += x[11]*44.009000; // CO2
// Extra 1e6 factor to take c to SI
ROW = 1e6*rho / XW;

// Compute conversion, see Eq 11
for (int id = 0; id < 12; ++id) {
c[id] = x[id]*ROW;
}

// convert to chemkin units
productionRate(wdot, c, T);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
wdot[id] *= 1.0e-6;
}
}

//  species unit charge number 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRG(int kcharge[])
{
kcharge[0] = 0;// NC5H12
kcharge[1] = 0;// NC6H14
kcharge[2] = 0;// NC8H18
kcharge[3] = 0;// IC8H18
kcharge[4] = 0;// NC10H22
kcharge[5] = 0;// NC7H16
kcharge[6] = 0;// NC12H26
kcharge[7] = 0;// O2
kcharge[8] = 0;// N2
kcharge[9] = 0;// NC16H34
kcharge[10] = 0;// H2O
kcharge[11] = 0;// CO2
}

//  species charge per unit mass 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRGMASS(amrex::Real zk[])
{

int kchrg[12];
CKCHRG(kchrg);

for (int id = 0; id < 12; ++id) {
zk[id] = 6.02214076e+23 * 1.60217663e-19 * kchrg[id] * imw(id);
}
}

// compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void dcvpRdT(amrex::Real * species, const amrex::Real *  tc)
{

// temperature
const amrex::Real T = tc[1];

// species with midpoint at T=1390 kelvin
if (T < 1390) {
// species 0: NC5H12
species[0] =
+6.07200973e-02
-7.15185522e-05 * tc[1]
+3.14721126e-08 * tc[2]
-4.85949260e-12 * tc[3]
;
// species 1: NC6H14
species[1] =
+7.29085608e-02
-8.77707838e-05 * tc[1]
+3.96938421e-08 * tc[2]
-6.33749692e-12 * tc[3]
;} else {
// species 0: NC5H12
species[0] =
+2.61086045e-02
-1.78194199e-05 * tc[1]
+4.14306744e-09 * tc[2]
-3.20118614e-13 * tc[3]
;
// species 1: NC6H14
species[1] =
+3.04480204e-02
-2.07589658e-05 * tc[1]
+4.82326371e-09 * tc[2]
-3.72507891e-13 * tc[3]
;}

// species with midpoint at T=1391 kelvin
if (T < 1391) {
// species 2: NC8H18
species[2] =
+9.78112063e-02
-1.21863672e-04 * tc[1]
+5.76016773e-08 * tc[2]
-9.71985000e-12 * tc[3]
;
// species 4: NC10H22
species[4] =
+1.22535012e-01
-1.55363148e-04 * tc[1]
+7.49504631e-08 * tc[2]
-1.29419215e-11 * tc[3]
;
// species 5: NC7H16
species[5] =
+8.54355820e-02
-1.05069357e-04 * tc[1]
+4.88837163e-08 * tc[2]
-8.09579700e-12 * tc[3]
;
// species 6: NC12H26
species[6] =
+1.47237711e-01
-1.88794054e-04 * tc[1]
+9.22323804e-08 * tc[2]
-1.61440892e-11 * tc[3]
;} else {
// species 2: NC8H18
species[2] =
+3.90887037e-02
-2.66077554e-05 * tc[1]
+6.17602581e-09 * tc[2]
-4.76668696e-13 * tc[3]
;
// species 4: NC10H22
species[4] =
+4.77244922e-02
-3.24552782e-05 * tc[1]
+7.52889777e-09 * tc[2]
-5.80863088e-13 * tc[3]
;
// species 5: NC7H16
species[5] =
+3.47675750e-02
-2.36814258e-05 * tc[1]
+5.49895434e-09 * tc[2]
-4.24521064e-13 * tc[3]
;
// species 6: NC12H26
species[6] =
+5.63550048e-02
-3.82986400e-05 * tc[1]
+8.88074586e-09 * tc[2]
-6.84976600e-13 * tc[3]
;}

// species with midpoint at T=1396 kelvin
if (T < 1396) {
// species 3: IC8H18
species[3] =
+1.11440581e-01
-1.58269316e-04 * tc[1]
+8.77218726e-08 * tc[2]
-1.77497276e-11 * tc[3]
;} else {
// species 3: IC8H18
species[3] =
+3.79004890e-02
-2.58874716e-05 * tc[1]
+6.02281116e-09 * tc[2]
-4.65602320e-13 * tc[3]
;}

// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 7: O2
species[7] =
-2.99673416e-03
+1.96946040e-05 * tc[1]
-2.90438853e-08 * tc[2]
+1.29749135e-11 * tc[3]
;
// species 8: N2
species[8] =
+1.40824000e-03
-7.92644400e-06 * tc[1]
+1.69245450e-08 * tc[2]
-9.77942000e-12 * tc[3]
;
// species 10: H2O
species[10] =
+3.47498246e-03
-1.27093927e-05 * tc[1]
+2.09057438e-08 * tc[2]
-1.00263539e-11 * tc[3]
;
// species 11: CO2
species[11] =
+8.98459677e-03
-1.42471254e-05 * tc[1]
+7.37757066e-09 * tc[2]
-5.74798192e-13 * tc[3]
;} else {
// species 7: O2
species[7] =
+1.48308754e-03
-1.51593334e-06 * tc[1]
+6.28411665e-10 * tc[2]
-8.66871176e-14 * tc[3]
;
// species 8: N2
species[8] =
+1.48797700e-03
-1.13695220e-06 * tc[1]
+3.02911200e-10 * tc[2]
-2.70134040e-14 * tc[3]
;
// species 10: H2O
species[10] =
+3.05629289e-03
-1.74605202e-06 * tc[1]
+3.60298917e-10 * tc[2]
-2.55664715e-14 * tc[3]
;
// species 11: CO2
species[11] =
+4.41437026e-03
-4.42962808e-06 * tc[1]
+1.57047056e-09 * tc[2]
-1.88833666e-13 * tc[3]
;}

// species with midpoint at T=1800 kelvin
if (T < 1800) {
// species 9: NC16H34
species[9] =
+1.96191732e-01
-2.53517518e-04 * tc[1]
+1.24093642e-07 * tc[2]
-2.16228183e-11 * tc[3]
;} else {
// species 9: NC16H34
species[9] =
+7.73881688e-02
-5.55115792e-05 * tc[1]
+1.40903431e-08 * tc[2]
-1.25183697e-12 * tc[3]
;}
}

// compute an approx to the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void aJacobian_precond(amrex::Real *  J, const amrex::Real *  /*sc*/, const amrex::Real /*T*/, const int /*HP*/)
{

#if defined(PELE_COMPILE_AJACOBIAN) || !defined(AMREX_USE_HIP)
for (int i=0; i<169; i++) {
J[i] = 0.0;
}
#else
amrex::Abort();
#endif
}

// compute an approx to the reaction Jacobian (for preconditioning)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void DWDOT_SIMPLIFIED(amrex::Real *  J, const amrex::Real *  sc, const amrex::Real *  Tp, const int * HP)
{
amrex::Real c[12];

for (int k=0; k<12; k++) {
c[k] = 1.e6 * sc[k];
}

aJacobian_precond(J, c, *Tp, *HP);

// dwdot[k]/dT
// dTdot/d[X]
for (int k=0; k<12; k++) {
J[156+k] *= 1.e-6;
J[k*13+12] *= 1.e6;
}

}


// compute the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void aJacobian(amrex::Real * J, const amrex::Real * /*sc*/, const amrex::Real /*T*/, const int /*consP*/)
{

#if defined(PELE_COMPILE_AJACOBIAN) || !defined(AMREX_USE_HIP)
for (int i=0; i<169; i++) {
J[i] = 0.0;
}
#else
amrex::Abort();
#endif
}

// compute the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void DWDOT(amrex::Real *  J, const amrex::Real *  sc, const amrex::Real *  Tp, const int * consP)
{
amrex::Real c[12];

for (int k=0; k<12; k++) {
c[k] = 1.e6 * sc[k];
}

aJacobian(J, c, *Tp, *consP);

// dwdot[k]/dT
// dTdot/d[X]
for (int k=0; k<12; k++) {
J[156+k] *= 1.e-6;
J[k*13+12] *= 1.e6;
}

}


// Transport function declarations 

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetLENIMC(int* LENIMC ) {
*LENIMC = 48;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetLENRMC(int* LENRMC ) {
*LENRMC = 3060;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNO(int* NO ) {
*NO = 4;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetKK(int* KK ) {
*KK = 12;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNLITE(int* NLITE ) {
*NLITE = 0;}


// Patm in ergs/cm3
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetPATM(amrex::Real* PATM) {
*PATM =   0.1013250000000000E+07;}

// the molecular weights in g/mol
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetWT(amrex::Real* WT ) {
WT[0] = 7.21510000E+01;
WT[1] = 8.61780000E+01;
WT[2] = 1.14232000E+02;
WT[3] = 1.14232000E+02;
WT[4] = 1.42286000E+02;
WT[5] = 1.00205000E+02;
WT[6] = 1.70340000E+02;
WT[7] = 3.19980000E+01;
WT[8] = 2.80140000E+01;
WT[9] = 2.26448000E+02;
WT[10] = 1.80150000E+01;
WT[11] = 4.40090000E+01;
}

// the lennard-jones potential well depth eps/kb in K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetEPS(amrex::Real* EPS ) {
EPS[0] = 4.58182000E+02;
EPS[1] = 5.12225000E+02;
EPS[2] = 6.13127000E+02;
EPS[3] = 4.58500000E+02;
EPS[4] = 7.04917000E+02;
EPS[5] = 4.59600000E+02;
EPS[6] = 7.89980000E+02;
EPS[7] = 1.07400000E+02;
EPS[8] = 9.75300000E+01;
EPS[9] = 6.50000000E+02;
EPS[10] = 5.72400000E+02;
EPS[11] = 2.44000000E+02;
}

// the lennard-jones collision diameter in Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetSIG(amrex::Real* SIG ) {
SIG[0] = 5.44500000E+00;
SIG[1] = 5.74200000E+00;
SIG[2] = 6.25000000E+00;
SIG[3] = 6.41400000E+00;
SIG[4] = 6.67500000E+00;
SIG[5] = 6.25300000E+00;
SIG[6] = 7.04700000E+00;
SIG[7] = 3.45800000E+00;
SIG[8] = 3.62100000E+00;
SIG[9] = 8.35000000E+00;
SIG[10] = 2.60500000E+00;
SIG[11] = 3.76300000E+00;
}

// the dipole moment in Debye
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetDIP(amrex::Real* DIP ) {
DIP[0] = 0.00000000E+00;
DIP[1] = 0.00000000E+00;
DIP[2] = 0.00000000E+00;
DIP[3] = 0.00000000E+00;
DIP[4] = 0.00000000E+00;
DIP[5] = 0.00000000E+00;
DIP[6] = 0.00000000E+00;
DIP[7] = 0.00000000E+00;
DIP[8] = 0.00000000E+00;
DIP[9] = 0.00000000E+00;
DIP[10] = 1.84400000E+00;
DIP[11] = 0.00000000E+00;
}

// the polarizability in cubic Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetPOL(amrex::Real* POL ) {
POL[0] = 0.00000000E+00;
POL[1] = 0.00000000E+00;
POL[2] = 0.00000000E+00;
POL[3] = 0.00000000E+00;
POL[4] = 0.00000000E+00;
POL[5] = 0.00000000E+00;
POL[6] = 0.00000000E+00;
POL[7] = 1.60000000E+00;
POL[8] = 1.76000000E+00;
POL[9] = 0.00000000E+00;
POL[10] = 0.00000000E+00;
POL[11] = 2.65000000E+00;
}

// the rotational relaxation collision number at 298 K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetZROT(amrex::Real* ZROT ) {
ZROT[0] = 1.00000000E+00;
ZROT[1] = 1.00000000E+00;
ZROT[2] = 1.00000000E+00;
ZROT[3] = 1.00000000E+00;
ZROT[4] = 1.00000000E+00;
ZROT[5] = 1.00000000E+00;
ZROT[6] = 1.00000000E+00;
ZROT[7] = 3.80000000E+00;
ZROT[8] = 4.00000000E+00;
ZROT[9] = 1.00000000E+00;
ZROT[10] = 4.00000000E+00;
ZROT[11] = 2.10000000E+00;
}

// 0: monoatomic, 1: linear, 2: nonlinear
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNLIN(int* NLIN) {
NLIN[0] = 2;
NLIN[1] = 2;
NLIN[2] = 2;
NLIN[3] = 2;
NLIN[4] = 2;
NLIN[5] = 2;
NLIN[6] = 2;
NLIN[7] = 1;
NLIN[8] = 1;
NLIN[9] = 2;
NLIN[10] = 2;
NLIN[11] = 1;
}


// Poly fits for the viscosities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFETA(amrex::Real* COFETA) {
COFETA[0] = -2.29005716E+01;
COFETA[1] = 3.92788495E+00;
COFETA[2] = -3.49780893E-01;
COFETA[3] = 1.21212758E-02;
COFETA[4] = -2.10271997E+01;
COFETA[5] = 3.03584039E+00;
COFETA[6] = -2.16607473E-01;
COFETA[7] = 5.68186005E-03;
COFETA[8] = -1.68951809E+01;
COFETA[9] = 1.14785929E+00;
COFETA[10] = 5.78369689E-02;
COFETA[11] = -7.31548437E-03;
COFETA[12] = -2.29889688E+01;
COFETA[13] = 3.92329359E+00;
COFETA[14] = -3.49082506E-01;
COFETA[15] = 1.20870516E-02;
COFETA[16] = -1.29489876E+01;
COFETA[17] = -6.04647944E-01;
COFETA[18] = 3.07320970E-01;
COFETA[19] = -1.89227139E-02;
COFETA[20] = -2.29712799E+01;
COFETA[21] = 3.90755399E+00;
COFETA[22] = -3.46688941E-01;
COFETA[23] = 1.19697999E-02;
COFETA[24] = -9.53395973E+00;
COFETA[25] = -2.09548278E+00;
COFETA[26] = 5.16271953E-01;
COFETA[27] = -2.85012704E-02;
COFETA[28] = -1.68118998E+01;
COFETA[29] = 2.52362554E+00;
COFETA[30] = -2.49309128E-01;
COFETA[31] = 1.10211025E-02;
COFETA[32] = -1.62526676E+01;
COFETA[33] = 2.24839597E+00;
COFETA[34] = -2.13428438E-01;
COFETA[35] = 9.46192413E-03;
COFETA[36] = -1.55685932E+01;
COFETA[37] = 4.51057467E-01;
COFETA[38] = 1.57706649E-01;
COFETA[39] = -1.19883346E-02;
COFETA[40] = -1.17770995E+01;
COFETA[41] = -8.26744826E-01;
COFETA[42] = 3.39009392E-01;
COFETA[43] = -2.00674472E-02;
COFETA[44] = -2.36749638E+01;
COFETA[45] = 4.99775518E+00;
COFETA[46] = -5.52687718E-01;
COFETA[47] = 2.34353338E-02;
}


// Poly fits for the conductivities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFLAM(amrex::Real* COFLAM) {
COFLAM[0] = -1.87743012E+01;
COFLAM[1] = 7.68887812E+00;
COFLAM[2] = -6.49527256E-01;
COFLAM[3] = 1.80038552E-02;
COFLAM[4] = -1.80140471E+01;
COFLAM[5] = 7.28886350E+00;
COFLAM[6] = -5.88065590E-01;
COFLAM[7] = 1.49879708E-02;
COFLAM[8] = -1.60169243E+01;
COFLAM[9] = 6.32682314E+00;
COFLAM[10] = -4.46230940E-01;
COFLAM[11] = 8.23381166E-03;
COFLAM[12] = -2.35634618E+01;
COFLAM[13] = 9.72687125E+00;
COFLAM[14] = -9.40839756E-01;
COFLAM[15] = 3.16881066E-02;
COFLAM[16] = -1.38317175E+01;
COFLAM[17] = 5.32397212E+00;
COFLAM[18] = -3.02248860E-01;
COFLAM[19] = 1.51838700E-03;
COFLAM[20] = -1.96434201E+01;
COFLAM[21] = 8.03403931E+00;
COFLAM[22] = -7.01574664E-01;
COFLAM[23] = 2.05526114E-02;
COFLAM[24] = -1.17923251E+01;
COFLAM[25] = 4.40896459E+00;
COFLAM[26] = -1.73024512E-01;
COFLAM[27] = -4.42005514E-03;
COFLAM[28] = -3.01283106E+00;
COFLAM[29] = 3.37554996E+00;
COFLAM[30] = -3.43353130E-01;
COFLAM[31] = 1.51043449E-02;
COFLAM[32] = 1.15507226E+01;
COFLAM[33] = -2.91453572E+00;
COFLAM[34] = 5.55045264E-01;
COFLAM[35] = -2.75173246E-02;
COFLAM[36] = -1.63275673E+01;
COFLAM[37] = 6.35061379E+00;
COFLAM[38] = -4.51182407E-01;
COFLAM[39] = 8.55345450E-03;
COFLAM[40] = 2.21730107E+01;
COFLAM[41] = -8.46933438E+00;
COFLAM[42] = 1.46153488E+00;
COFLAM[43] = -7.29500808E-02;
COFLAM[44] = -1.24047472E+01;
COFLAM[45] = 6.34783072E+00;
COFLAM[46] = -6.37857831E-01;
COFLAM[47] = 2.37613806E-02;
}


// Poly fits for the diffusion coefficients, dim NO*KK*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFD(amrex::Real* COFD) {
COFD[0] = -2.20396859E+01;
COFD[1] = 4.80541922E+00;
COFD[2] = -3.42753318E-01;
COFD[3] = 1.22750163E-02;
COFD[4] = -2.16073232E+01;
COFD[5] = 4.53847670E+00;
COFD[6] = -3.01502124E-01;
COFD[7] = 1.02366618E-02;
COFD[8] = -2.07682159E+01;
COFD[9] = 4.05314295E+00;
COFD[10] = -2.27706292E-01;
COFD[11] = 6.63337598E-03;
COFD[12] = -2.23094096E+01;
COFD[13] = 4.80406750E+00;
COFD[14] = -3.42539193E-01;
COFD[15] = 1.22642685E-02;
COFD[16] = -1.99004183E+01;
COFD[17] = 3.58882517E+00;
COFD[18] = -1.58571243E-01;
COFD[19] = 3.30756307E-03;
COFD[20] = -2.22468624E+01;
COFD[21] = 4.79938990E+00;
COFD[22] = -3.41798465E-01;
COFD[23] = 1.22270965E-02;
COFD[24] = -1.90130516E+01;
COFD[25] = 3.13800920E+00;
COFD[26] = -9.25506586E-02;
COFD[27] = 1.69413653E-04;
COFD[28] = -2.08991500E+01;
COFD[29] = 5.12777022E+00;
COFD[30] = -4.39040640E-01;
COFD[31] = 1.85511731E-02;
COFD[32] = -2.05259655E+01;
COFD[33] = 5.00332689E+00;
COFD[34] = -4.24470438E-01;
COFD[35] = 1.79814505E-02;
COFD[36] = -2.08372255E+01;
COFD[37] = 3.87614972E+00;
COFD[38] = -2.01138578E-01;
COFD[39] = 5.34840295E-03;
COFD[40] = -1.98446554E+01;
COFD[41] = 4.24085889E+00;
COFD[42] = -2.56126131E-01;
COFD[43] = 8.01556198E-03;
COFD[44] = -2.26675871E+01;
COFD[45] = 5.51808537E+00;
COFD[46] = -4.65712484E-01;
COFD[47] = 1.87719102E-02;
COFD[48] = -2.16073232E+01;
COFD[49] = 4.53847670E+00;
COFD[50] = -3.01502124E-01;
COFD[51] = 1.02366618E-02;
COFD[52] = -2.11005299E+01;
COFD[53] = 4.23978495E+00;
COFD[54] = -2.55962392E-01;
COFD[55] = 8.00756055E-03;
COFD[56] = -2.01139357E+01;
COFD[57] = 3.69246340E+00;
COFD[58] = -1.73863180E-01;
COFD[59] = 4.03859662E-03;
COFD[60] = -2.18815094E+01;
COFD[61] = 4.53658146E+00;
COFD[62] = -3.01212880E-01;
COFD[63] = 1.02224813E-02;
COFD[64] = -1.90458224E+01;
COFD[65] = 3.14563115E+00;
COFD[66] = -9.36748910E-02;
COFD[67] = 2.22927330E-04;
COFD[68] = -2.18123243E+01;
COFD[69] = 4.53003044E+00;
COFD[70] = -3.00213235E-01;
COFD[71] = 1.01734776E-02;
COFD[72] = -1.83188596E+01;
COFD[73] = 2.76235018E+00;
COFD[74] = -3.72816184E-02;
COFD[75] = -2.45837985E-03;
COFD[76] = -2.13316775E+01;
COFD[77] = 5.23707202E+00;
COFD[78] = -4.50870877E-01;
COFD[79] = 1.89710559E-02;
COFD[80] = -2.10312389E+01;
COFD[81] = 5.14647415E+00;
COFD[82] = -4.41208890E-01;
COFD[83] = 1.86350422E-02;
COFD[84] = -2.00936890E+01;
COFD[85] = 3.47661962E+00;
COFD[86] = -1.42053106E-01;
COFD[87] = 2.51950342E-03;
COFD[88] = -1.92911861E+01;
COFD[89] = 3.92801999E+00;
COFD[90] = -2.08870639E-01;
COFD[91] = 5.72067207E-03;
COFD[92] = -2.28316747E+01;
COFD[93] = 5.50566151E+00;
COFD[94] = -4.60020840E-01;
COFD[95] = 1.83701806E-02;
COFD[96] = -2.07682159E+01;
COFD[97] = 4.05314295E+00;
COFD[98] = -2.27706292E-01;
COFD[99] = 6.63337598E-03;
COFD[100] = -2.01139357E+01;
COFD[101] = 3.69246340E+00;
COFD[102] = -1.73863180E-01;
COFD[103] = 4.03859662E-03;
COFD[104] = -1.88215225E+01;
COFD[105] = 3.01899523E+00;
COFD[106] = -7.49981344E-02;
COFD[107] = -6.65959618E-04;
COFD[108] = -2.10509667E+01;
COFD[109] = 4.05101803E+00;
COFD[110] = -2.27386254E-01;
COFD[111] = 6.61784537E-03;
COFD[112] = -1.77893546E+01;
COFD[113] = 2.48796358E+00;
COFD[114] = 2.47414627E-03;
COFD[115] = -4.32916021E-03;
COFD[116] = -2.09766491E+01;
COFD[117] = 4.04380664E+00;
COFD[118] = -2.26298990E-01;
COFD[119] = 6.56504901E-03;
COFD[120] = -1.67927381E+01;
COFD[121] = 1.99564934E+00;
COFD[122] = 7.32168356E-02;
COFD[123] = -7.63548809E-03;
COFD[124] = -2.19478238E+01;
COFD[125] = 5.37866722E+00;
COFD[126] = -4.64853773E-01;
COFD[127] = 1.94041759E-02;
COFD[128] = -2.16548852E+01;
COFD[129] = 5.29519757E+00;
COFD[130] = -4.56363642E-01;
COFD[131] = 1.91281929E-02;
COFD[132] = -1.88607802E+01;
COFD[133] = 2.82475565E+00;
COFD[134] = -4.63864261E-02;
COFD[135] = -2.02729765E-03;
COFD[136] = -1.80014844E+01;
COFD[137] = 3.26089295E+00;
COFD[138] = -1.10490863E-01;
COFD[139] = 1.01966774E-03;
COFD[140] = -2.27782552E+01;
COFD[141] = 5.34930975E+00;
COFD[142] = -4.31287593E-01;
COFD[143] = 1.67983303E-02;
COFD[144] = -2.23094096E+01;
COFD[145] = 4.80406750E+00;
COFD[146] = -3.42539193E-01;
COFD[147] = 1.22642685E-02;
COFD[148] = -2.18815094E+01;
COFD[149] = 4.53658146E+00;
COFD[150] = -3.01212880E-01;
COFD[151] = 1.02224813E-02;
COFD[152] = -2.10509667E+01;
COFD[153] = 4.05101803E+00;
COFD[154] = -2.27386254E-01;
COFD[155] = 6.61784537E-03;
COFD[156] = -2.25918928E+01;
COFD[157] = 4.80271460E+00;
COFD[158] = -3.42324914E-01;
COFD[159] = 1.22535140E-02;
COFD[160] = -2.01889639E+01;
COFD[161] = 3.58629280E+00;
COFD[162] = -1.58198171E-01;
COFD[163] = 3.28975875E-03;
COFD[164] = -2.25240439E+01;
COFD[165] = 4.79806454E+00;
COFD[166] = -3.41588375E-01;
COFD[167] = 1.22165519E-02;
COFD[168] = -1.93075642E+01;
COFD[169] = 3.13582885E+00;
COFD[170] = -9.22291119E-02;
COFD[171] = 1.54109619E-04;
COFD[172] = -2.11685652E+01;
COFD[173] = 5.12871609E+00;
COFD[174] = -4.39152449E-01;
COFD[175] = 1.85555947E-02;
COFD[176] = -2.07853786E+01;
COFD[177] = 5.00395099E+00;
COFD[178] = -4.24537954E-01;
COFD[179] = 1.79838424E-02;
COFD[180] = -2.11319430E+01;
COFD[181] = 3.87383924E+00;
COFD[182] = -2.00794381E-01;
COFD[183] = 5.33183840E-03;
COFD[184] = -2.01066507E+01;
COFD[185] = 4.23909888E+00;
COFD[186] = -2.55857796E-01;
COFD[187] = 8.00244944E-03;
COFD[188] = -2.29441492E+01;
COFD[189] = 5.51837255E+00;
COFD[190] = -4.65728192E-01;
COFD[191] = 1.87717847E-02;
COFD[192] = -1.99004183E+01;
COFD[193] = 3.58882517E+00;
COFD[194] = -1.58571243E-01;
COFD[195] = 3.30756307E-03;
COFD[196] = -1.90458224E+01;
COFD[197] = 3.14563115E+00;
COFD[198] = -9.36748910E-02;
COFD[199] = 2.22927330E-04;
COFD[200] = -1.77893546E+01;
COFD[201] = 2.48796358E+00;
COFD[202] = 2.47414627E-03;
COFD[203] = -4.32916021E-03;
COFD[204] = -2.01889639E+01;
COFD[205] = 3.58629280E+00;
COFD[206] = -1.58198171E-01;
COFD[207] = 3.28975875E-03;
COFD[208] = -1.65247183E+01;
COFD[209] = 1.86376261E+00;
COFD[210] = 9.20701178E-02;
COFD[211] = -8.51470516E-03;
COFD[212] = -2.01080776E+01;
COFD[213] = 3.57754981E+00;
COFD[214] = -1.56910225E-01;
COFD[215] = 3.22829712E-03;
COFD[216] = -1.53428459E+01;
COFD[217] = 1.29341827E+00;
COFD[218] = 1.73454228E-01;
COFD[219] = -1.23007501E-02;
COFD[220] = -2.24759923E+01;
COFD[221] = 5.50822096E+00;
COFD[222] = -4.78234659E-01;
COFD[223] = 1.98512785E-02;
COFD[224] = -2.21643795E+01;
COFD[225] = 5.41977231E+00;
COFD[226] = -4.69168499E-01;
COFD[227] = 1.95518107E-02;
COFD[228] = -1.77303626E+01;
COFD[229] = 2.25150083E+00;
COFD[230] = 3.65779467E-02;
COFD[231] = -5.92687944E-03;
COFD[232] = -1.70666161E+01;
COFD[233] = 2.77459742E+00;
COFD[234] = -3.90660125E-02;
COFD[235] = -2.37392826E-03;
COFD[236] = -2.25761704E+01;
COFD[237] = 5.16216553E+00;
COFD[238] = -3.99975925E-01;
COFD[239] = 1.51719252E-02;
COFD[240] = -2.22468624E+01;
COFD[241] = 4.79938990E+00;
COFD[242] = -3.41798465E-01;
COFD[243] = 1.22270965E-02;
COFD[244] = -2.18123243E+01;
COFD[245] = 4.53003044E+00;
COFD[246] = -3.00213235E-01;
COFD[247] = 1.01734776E-02;
COFD[248] = -2.09766491E+01;
COFD[249] = 4.04380664E+00;
COFD[250] = -2.26298990E-01;
COFD[251] = 6.56504901E-03;
COFD[252] = -2.25240439E+01;
COFD[253] = 4.79806454E+00;
COFD[254] = -3.41588375E-01;
COFD[255] = 1.22165519E-02;
COFD[256] = -2.01080776E+01;
COFD[257] = 3.57754981E+00;
COFD[258] = -1.56910225E-01;
COFD[259] = 3.22829712E-03;
COFD[260] = -2.24581406E+01;
COFD[261] = 4.79347331E+00;
COFD[262] = -3.40860948E-01;
COFD[263] = 1.21800530E-02;
COFD[264] = -1.92271958E+01;
COFD[265] = 3.12826317E+00;
COFD[266] = -9.11138226E-02;
COFD[267] = 1.01042801E-04;
COFD[268] = -2.11294239E+01;
COFD[269] = 5.13172698E+00;
COFD[270] = -4.39502623E-01;
COFD[271] = 1.85691915E-02;
COFD[272] = -2.07461896E+01;
COFD[273] = 5.00611298E+00;
COFD[274] = -4.24772056E-01;
COFD[275] = 1.79921460E-02;
COFD[276] = -2.10486881E+01;
COFD[277] = 3.86584757E+00;
COFD[278] = -1.99604029E-01;
COFD[279] = 5.27455854E-03;
COFD[280] = -2.00488777E+01;
COFD[281] = 4.23300736E+00;
COFD[282] = -2.54929309E-01;
COFD[283] = 7.95708569E-03;
COFD[284] = -2.28975481E+01;
COFD[285] = 5.51934005E+00;
COFD[286] = -4.65779010E-01;
COFD[287] = 1.87711908E-02;
COFD[288] = -1.90130516E+01;
COFD[289] = 3.13800920E+00;
COFD[290] = -9.25506586E-02;
COFD[291] = 1.69413653E-04;
COFD[292] = -1.83188596E+01;
COFD[293] = 2.76235018E+00;
COFD[294] = -3.72816184E-02;
COFD[295] = -2.45837985E-03;
COFD[296] = -1.67927381E+01;
COFD[297] = 1.99564934E+00;
COFD[298] = 7.32168356E-02;
COFD[299] = -7.63548809E-03;
COFD[300] = -1.93075642E+01;
COFD[301] = 3.13582885E+00;
COFD[302] = -9.22291119E-02;
COFD[303] = 1.54109619E-04;
COFD[304] = -1.53428459E+01;
COFD[305] = 1.29341827E+00;
COFD[306] = 1.73454228E-01;
COFD[307] = -1.23007501E-02;
COFD[308] = -1.92271958E+01;
COFD[309] = 3.12826317E+00;
COFD[310] = -9.11138226E-02;
COFD[311] = 1.01042801E-04;
COFD[312] = -1.44331280E+01;
COFD[313] = 8.41942943E-01;
COFD[314] = 2.37470782E-01;
COFD[315] = -1.52471999E-02;
COFD[316] = -2.26775238E+01;
COFD[317] = 5.51592833E+00;
COFD[318] = -4.75325322E-01;
COFD[319] = 1.95715263E-02;
COFD[320] = -2.25658626E+01;
COFD[321] = 5.51389289E+00;
COFD[322] = -4.78413122E-01;
COFD[323] = 1.98371302E-02;
COFD[324] = -1.65777457E+01;
COFD[325] = 1.69078901E+00;
COFD[326] = 1.16782972E-01;
COFD[327] = -9.66685173E-03;
COFD[328] = -1.61225000E+01;
COFD[329] = 2.30448095E+00;
COFD[330] = 2.89601588E-02;
COFD[331] = -5.57107525E-03;
COFD[332] = -2.23094683E+01;
COFD[333] = 4.96580377E+00;
COFD[334] = -3.68253958E-01;
COFD[335] = 1.35580688E-02;
COFD[336] = -2.08991500E+01;
COFD[337] = 5.12777022E+00;
COFD[338] = -4.39040640E-01;
COFD[339] = 1.85511731E-02;
COFD[340] = -2.13316775E+01;
COFD[341] = 5.23707202E+00;
COFD[342] = -4.50870877E-01;
COFD[343] = 1.89710559E-02;
COFD[344] = -2.19478238E+01;
COFD[345] = 5.37866722E+00;
COFD[346] = -4.64853773E-01;
COFD[347] = 1.94041759E-02;
COFD[348] = -2.11685652E+01;
COFD[349] = 5.12871609E+00;
COFD[350] = -4.39152449E-01;
COFD[351] = 1.85555947E-02;
COFD[352] = -2.24759923E+01;
COFD[353] = 5.50822096E+00;
COFD[354] = -4.78234659E-01;
COFD[355] = 1.98512785E-02;
COFD[356] = -2.11294239E+01;
COFD[357] = 5.13172698E+00;
COFD[358] = -4.39502623E-01;
COFD[359] = 1.85691915E-02;
COFD[360] = -2.26775238E+01;
COFD[361] = 5.51592833E+00;
COFD[362] = -4.75325322E-01;
COFD[363] = 1.95715263E-02;
COFD[364] = -1.53110597E+01;
COFD[365] = 3.37317428E+00;
COFD[366] = -2.24900439E-01;
COFD[367] = 9.81228151E-03;
COFD[368] = -1.50096252E+01;
COFD[369] = 3.25515933E+00;
COFD[370] = -2.09710110E-01;
COFD[371] = 9.15941830E-03;
COFD[372] = -2.25795432E+01;
COFD[373] = 5.43595045E+00;
COFD[374] = -4.70905585E-01;
COFD[375] = 1.96130544E-02;
COFD[376] = -2.10639926E+01;
COFD[377] = 5.50980703E+00;
COFD[378] = -4.78335492E-01;
COFD[379] = 1.98515434E-02;
COFD[380] = -1.81197250E+01;
COFD[381] = 4.33684042E+00;
COFD[382] = -3.44981265E-01;
COFD[383] = 1.48142449E-02;
COFD[384] = -2.05259655E+01;
COFD[385] = 5.00332689E+00;
COFD[386] = -4.24470438E-01;
COFD[387] = 1.79814505E-02;
COFD[388] = -2.10312389E+01;
COFD[389] = 5.14647415E+00;
COFD[390] = -4.41208890E-01;
COFD[391] = 1.86350422E-02;
COFD[392] = -2.16548852E+01;
COFD[393] = 5.29519757E+00;
COFD[394] = -4.56363642E-01;
COFD[395] = 1.91281929E-02;
COFD[396] = -2.07853786E+01;
COFD[397] = 5.00395099E+00;
COFD[398] = -4.24537954E-01;
COFD[399] = 1.79838424E-02;
COFD[400] = -2.21643795E+01;
COFD[401] = 5.41977231E+00;
COFD[402] = -4.69168499E-01;
COFD[403] = 1.95518107E-02;
COFD[404] = -2.07461896E+01;
COFD[405] = 5.00611298E+00;
COFD[406] = -4.24772056E-01;
COFD[407] = 1.79921460E-02;
COFD[408] = -2.25658626E+01;
COFD[409] = 5.51389289E+00;
COFD[410] = -4.78413122E-01;
COFD[411] = 1.98371302E-02;
COFD[412] = -1.50096252E+01;
COFD[413] = 3.25515933E+00;
COFD[414] = -2.09710110E-01;
COFD[415] = 9.15941830E-03;
COFD[416] = -1.47639411E+01;
COFD[417] = 3.15955654E+00;
COFD[418] = -1.97590757E-01;
COFD[419] = 8.64692156E-03;
COFD[420] = -2.22525051E+01;
COFD[421] = 5.34409219E+00;
COFD[422] = -4.61334919E-01;
COFD[423] = 1.92897986E-02;
COFD[424] = -2.08123331E+01;
COFD[425] = 5.42470175E+00;
COFD[426] = -4.69700438E-01;
COFD[427] = 1.95706912E-02;
COFD[428] = -1.77350630E+01;
COFD[429] = 4.19328271E+00;
COFD[430] = -3.26911461E-01;
COFD[431] = 1.40520357E-02;
COFD[432] = -2.08372255E+01;
COFD[433] = 3.87614972E+00;
COFD[434] = -2.01138578E-01;
COFD[435] = 5.34840295E-03;
COFD[436] = -2.00936890E+01;
COFD[437] = 3.47661962E+00;
COFD[438] = -1.42053106E-01;
COFD[439] = 2.51950342E-03;
COFD[440] = -1.88607802E+01;
COFD[441] = 2.82475565E+00;
COFD[442] = -4.63864261E-02;
COFD[443] = -2.02729765E-03;
COFD[444] = -2.11319430E+01;
COFD[445] = 3.87383924E+00;
COFD[446] = -2.00794381E-01;
COFD[447] = 5.33183840E-03;
COFD[448] = -1.77303626E+01;
COFD[449] = 2.25150083E+00;
COFD[450] = 3.65779467E-02;
COFD[451] = -5.92687944E-03;
COFD[452] = -2.10486881E+01;
COFD[453] = 3.86584757E+00;
COFD[454] = -1.99604029E-01;
COFD[455] = 5.27455854E-03;
COFD[456] = -1.65777457E+01;
COFD[457] = 1.69078901E+00;
COFD[458] = 1.16782972E-01;
COFD[459] = -9.66685173E-03;
COFD[460] = -2.25795432E+01;
COFD[461] = 5.43595045E+00;
COFD[462] = -4.70905585E-01;
COFD[463] = 1.96130544E-02;
COFD[464] = -2.22525051E+01;
COFD[465] = 5.34409219E+00;
COFD[466] = -4.61334919E-01;
COFD[467] = 1.92897986E-02;
COFD[468] = -1.87949917E+01;
COFD[469] = 2.58040365E+00;
COFD[470] = -1.08977961E-02;
COFD[471] = -3.70103387E-03;
COFD[472] = -1.80095597E+01;
COFD[473] = 3.05076642E+00;
COFD[474] = -7.96958601E-02;
COFD[475] = -4.41991430E-04;
COFD[476] = -2.31273254E+01;
COFD[477] = 5.28177733E+00;
COFD[478] = -4.19694755E-01;
COFD[479] = 1.61874765E-02;
COFD[480] = -1.98446554E+01;
COFD[481] = 4.24085889E+00;
COFD[482] = -2.56126131E-01;
COFD[483] = 8.01556198E-03;
COFD[484] = -1.92911861E+01;
COFD[485] = 3.92801999E+00;
COFD[486] = -2.08870639E-01;
COFD[487] = 5.72067207E-03;
COFD[488] = -1.80014844E+01;
COFD[489] = 3.26089295E+00;
COFD[490] = -1.10490863E-01;
COFD[491] = 1.01966774E-03;
COFD[492] = -2.01066507E+01;
COFD[493] = 4.23909888E+00;
COFD[494] = -2.55857796E-01;
COFD[495] = 8.00244944E-03;
COFD[496] = -1.70666161E+01;
COFD[497] = 2.77459742E+00;
COFD[498] = -3.90660125E-02;
COFD[499] = -2.37392826E-03;
COFD[500] = -2.00488777E+01;
COFD[501] = 4.23300736E+00;
COFD[502] = -2.54929309E-01;
COFD[503] = 7.95708569E-03;
COFD[504] = -1.61225000E+01;
COFD[505] = 2.30448095E+00;
COFD[506] = 2.89601588E-02;
COFD[507] = -5.57107525E-03;
COFD[508] = -2.10639926E+01;
COFD[509] = 5.50980703E+00;
COFD[510] = -4.78335492E-01;
COFD[511] = 1.98515434E-02;
COFD[512] = -2.08123331E+01;
COFD[513] = 5.42470175E+00;
COFD[514] = -4.69700438E-01;
COFD[515] = 1.95706912E-02;
COFD[516] = -1.80095597E+01;
COFD[517] = 3.05076642E+00;
COFD[518] = -7.96958601E-02;
COFD[519] = -4.41991430E-04;
COFD[520] = -1.31492509E+01;
COFD[521] = 1.48004045E+00;
COFD[522] = 1.60499956E-01;
COFD[523] = -1.19765871E-02;
COFD[524] = -2.12021420E+01;
COFD[525] = 5.20775027E+00;
COFD[526] = -4.07348285E-01;
COFD[527] = 1.55473262E-02;
COFD[528] = -2.26675871E+01;
COFD[529] = 5.51808537E+00;
COFD[530] = -4.65712484E-01;
COFD[531] = 1.87719102E-02;
COFD[532] = -2.28316747E+01;
COFD[533] = 5.50566151E+00;
COFD[534] = -4.60020840E-01;
COFD[535] = 1.83701806E-02;
COFD[536] = -2.27782552E+01;
COFD[537] = 5.34930975E+00;
COFD[538] = -4.31287593E-01;
COFD[539] = 1.67983303E-02;
COFD[540] = -2.29441492E+01;
COFD[541] = 5.51837255E+00;
COFD[542] = -4.65728192E-01;
COFD[543] = 1.87717847E-02;
COFD[544] = -2.25761704E+01;
COFD[545] = 5.16216553E+00;
COFD[546] = -3.99975925E-01;
COFD[547] = 1.51719252E-02;
COFD[548] = -2.28975481E+01;
COFD[549] = 5.51934005E+00;
COFD[550] = -4.65779010E-01;
COFD[551] = 1.87711908E-02;
COFD[552] = -2.23094683E+01;
COFD[553] = 4.96580377E+00;
COFD[554] = -3.68253958E-01;
COFD[555] = 1.35580688E-02;
COFD[556] = -1.81197250E+01;
COFD[557] = 4.33684042E+00;
COFD[558] = -3.44981265E-01;
COFD[559] = 1.48142449E-02;
COFD[560] = -1.77350630E+01;
COFD[561] = 4.19328271E+00;
COFD[562] = -3.26911461E-01;
COFD[563] = 1.40520357E-02;
COFD[564] = -2.31273254E+01;
COFD[565] = 5.28177733E+00;
COFD[566] = -4.19694755E-01;
COFD[567] = 1.61874765E-02;
COFD[568] = -2.12021420E+01;
COFD[569] = 5.20775027E+00;
COFD[570] = -4.07348285E-01;
COFD[571] = 1.55473262E-02;
COFD[572] = -2.10907633E+01;
COFD[573] = 5.29211327E+00;
COFD[574] = -4.56068366E-01;
COFD[575] = 1.91195062E-02;
}


// List of specs with small weight, dim NLITE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetKTDIF(int* /*KTDIF*/) {
}


// Poly fits for thermal diff ratios, dim NO*NLITE*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFTD(amrex::Real* /*COFTD*/) {
}


// compute the critical parameters for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void GET_CRITPARAMS(amrex::Real *  Tci, amrex::Real *  ai, amrex::Real *  bi, amrex::Real *  acentric_i)
{

amrex::Real   EPS[12];
amrex::Real   SIG[12];
amrex::Real    wt[12];
amrex::Real Rcst = 83.144598;// in bar [CGS] !
amrex::Real avogadro = 6.02214199e23;
amrex::Real boltzmann = 1.3806503e-16;// we work in CGS

egtransetEPS(EPS);
egtransetSIG(SIG);
get_mw(wt);

// species 0: NC5H12
Tci[0] = 1.316 * EPS[0] ; 
ai[0] = (5.55 * avogadro * avogadro * EPS[0]*boltzmann * 1e-24 * SIG[0] * SIG[0] * SIG[0] ) / (wt[0] * wt[0]); 
bi[0] = 0.855 * avogadro * 1e-24 * SIG[0] * SIG[0] * SIG[0] / (wt[0]); 
acentric_i[0] = 0.0 ;

// species 1: NC6H14
Tci[1] = 1.316 * EPS[1] ; 
ai[1] = (5.55 * avogadro * avogadro * EPS[1]*boltzmann * 1e-24 * SIG[1] * SIG[1] * SIG[1] ) / (wt[1] * wt[1]); 
bi[1] = 0.855 * avogadro * 1e-24 * SIG[1] * SIG[1] * SIG[1] / (wt[1]); 
acentric_i[1] = 0.0 ;

// species 2: NC8H18
Tci[2] = 1.316 * EPS[2] ; 
ai[2] = (5.55 * avogadro * avogadro * EPS[2]*boltzmann * 1e-24 * SIG[2] * SIG[2] * SIG[2] ) / (wt[2] * wt[2]); 
bi[2] = 0.855 * avogadro * 1e-24 * SIG[2] * SIG[2] * SIG[2] / (wt[2]); 
acentric_i[2] = 0.0 ;

// species 3: IC8H18
Tci[3] = 1.316 * EPS[3] ; 
ai[3] = (5.55 * avogadro * avogadro * EPS[3]*boltzmann * 1e-24 * SIG[3] * SIG[3] * SIG[3] ) / (wt[3] * wt[3]); 
bi[3] = 0.855 * avogadro * 1e-24 * SIG[3] * SIG[3] * SIG[3] / (wt[3]); 
acentric_i[3] = 0.0 ;

// species 4: NC10H22
Tci[4] = 1.316 * EPS[4] ; 
ai[4] = (5.55 * avogadro * avogadro * EPS[4]*boltzmann * 1e-24 * SIG[4] * SIG[4] * SIG[4] ) / (wt[4] * wt[4]); 
bi[4] = 0.855 * avogadro * 1e-24 * SIG[4] * SIG[4] * SIG[4] / (wt[4]); 
acentric_i[4] = 0.0 ;

// species 5: NC7H16
Tci[5] = 1.316 * EPS[5] ; 
ai[5] = (5.55 * avogadro * avogadro * EPS[5]*boltzmann * 1e-24 * SIG[5] * SIG[5] * SIG[5] ) / (wt[5] * wt[5]); 
bi[5] = 0.855 * avogadro * 1e-24 * SIG[5] * SIG[5] * SIG[5] / (wt[5]); 
acentric_i[5] = 0.0 ;

// species 6: NC12H26
Tci[6] = 1.316 * EPS[6] ; 
ai[6] = (5.55 * avogadro * avogadro * EPS[6]*boltzmann * 1e-24 * SIG[6] * SIG[6] * SIG[6] ) / (wt[6] * wt[6]); 
bi[6] = 0.855 * avogadro * 1e-24 * SIG[6] * SIG[6] * SIG[6] / (wt[6]); 
acentric_i[6] = 0.0 ;

// species 7: O2
// Imported from NIST
Tci[7] = 154.581000 ; 
ai[7] = 1e6 * 0.42748 * Rcst * Rcst * Tci[7] * Tci[7] / (31.998800 * 31.998800 * 50.430466); 
bi[7] = 0.08664 * Rcst * Tci[7] / (31.998800 * 50.430466); 
acentric_i[7] = 0.022200 ;

// species 8: N2
// Imported from NIST
Tci[8] = 126.192000 ; 
ai[8] = 1e6 * 0.42748 * Rcst * Rcst * Tci[8] * Tci[8] / (28.013400 * 28.013400 * 33.958000); 
bi[8] = 0.08664 * Rcst * Tci[8] / (28.013400 * 33.958000); 
acentric_i[8] = 0.037200 ;

// species 9: NC16H34
Tci[9] = 1.316 * EPS[9] ; 
ai[9] = (5.55 * avogadro * avogadro * EPS[9]*boltzmann * 1e-24 * SIG[9] * SIG[9] * SIG[9] ) / (wt[9] * wt[9]); 
bi[9] = 0.855 * avogadro * 1e-24 * SIG[9] * SIG[9] * SIG[9] / (wt[9]); 
acentric_i[9] = 0.0 ;

// species 10: H2O
// Imported from NIST
Tci[10] = 647.096000 ; 
ai[10] = 1e6 * 0.42748 * Rcst * Rcst * Tci[10] * Tci[10] / (18.015340 * 18.015340 * 220.640000); 
bi[10] = 0.08664 * Rcst * Tci[10] / (18.015340 * 220.640000); 
acentric_i[10] = 0.344300 ;

// species 11: CO2
// Imported from NIST
Tci[11] = 304.120000 ; 
ai[11] = 1e6 * 0.42748 * Rcst * Rcst * Tci[11] * Tci[11] / (44.009950 * 44.009950 * 73.740000); 
bi[11] = 0.08664 * Rcst * Tci[11] / (44.009950 * 73.740000); 
acentric_i[11] = 0.225000 ;

}

// gauss-jordan solver external routine
// Replace this routine with the one generated by the Gauss Jordan solver of DW
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void sgjsolve(amrex::Real* /*A*/, amrex::Real* /*x*/, amrex::Real* /*b*/) {
amrex::Abort("sgjsolve not implemented, choose a different solver ");
}

// Replace this routine with the one generated by the Gauss Jordan solver of DW
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void sgjsolve_simplified(amrex::Real* /*A*/, amrex::Real* /*x*/, amrex::Real* /*b*/) {
amrex::Abort("sgjsolve_simplified not implemented, choose a different solver ");
}

#endif
