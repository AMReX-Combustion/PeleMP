#include "mechanism.H"

// Returns 0-based map of reaction order
void GET_RMAP
(int * /*_rmap*/)
{
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void CKINU(const int i, int& nspec, int* /*ki*/, int* /*nu*/)
{
if (i < 1) {
// Return max num species per reaction
nspec = 0;
} else {
nspec = -1;
}
}

// Returns the progress rates of each reactions
// Given P, T, and mole fractions
void CKKFKR(const amrex::Real P, const amrex::Real T, const amrex::Real x[], amrex::Real q_f[], amrex::Real q_r[])
{
amrex::Real c[12]; // temporary storage
amrex::Real PORT = 1e6 * P/(8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

// Compute conversion, see Eq 10
for (int id = 0; id < 12; ++id) {
c[id] = x[id]*PORT;
}

// convert to chemkin units
progressRateFR(q_f, q_r, c, T);

}

// compute the progress rate for each reaction
// USES progressRate : todo switch to GPU
void progressRateFR(amrex::Real *  q_f, amrex::Real *  q_r, amrex::Real *  sc, amrex::Real T)
{
}

// save atomic weights into array
void atomicWeight(amrex::Real *  awt)
{
awt[0] = 12.011000; // C
awt[1] = 1.008000; // H
awt[2] = 15.999000; // O
awt[3] = 14.007000; // N
}

// get atomic weight for all elements
void CKAWT( amrex::Real *  awt)
{
atomicWeight(awt);
}

// Returns the elemental composition 
// of the speciesi (mdim is num of elements)
void CKNCF(int * ncf)
{
int kd = 4; 
// Zero ncf
for (int id = 0; id < kd * 12; ++ id) {
 ncf[id] = 0; 
}

// NC5H12
ncf[ 0 * kd + 0 ] = 5; // C
ncf[ 0 * kd + 1 ] = 12; // H

// NC6H14
ncf[ 1 * kd + 0 ] = 6; // C
ncf[ 1 * kd + 1 ] = 14; // H

// NC8H18
ncf[ 2 * kd + 0 ] = 8; // C
ncf[ 2 * kd + 1 ] = 18; // H

// IC8H18
ncf[ 3 * kd + 0 ] = 8; // C
ncf[ 3 * kd + 1 ] = 18; // H

// NC10H22
ncf[ 4 * kd + 0 ] = 10; // C
ncf[ 4 * kd + 1 ] = 22; // H

// NC7H16
ncf[ 5 * kd + 0 ] = 7; // C
ncf[ 5 * kd + 1 ] = 16; // H

// NC12H26
ncf[ 6 * kd + 0 ] = 12; // C
ncf[ 6 * kd + 1 ] = 26; // H

// O2
ncf[ 7 * kd + 2 ] = 2; // O

// N2
ncf[ 8 * kd + 3 ] = 2; // N

// NC16H34
ncf[ 9 * kd + 0 ] = 16; // C
ncf[ 9 * kd + 1 ] = 34; // H

// H2O
ncf[ 10 * kd + 1 ] = 2; // H
ncf[ 10 * kd + 2 ] = 1; // O

// CO2
ncf[ 11 * kd + 0 ] = 1; // C
ncf[ 11 * kd + 2 ] = 2; // O

}

// Returns the vector of strings of element names
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
ename.resize(4);
ename[0] = "C";
ename[1] = "H";
ename[2] = "O";
ename[3] = "N";
}

// Returns the vector of strings of species names
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
kname.resize(12);
kname[0] = "NC5H12";
kname[1] = "NC6H14";
kname[2] = "NC8H18";
kname[3] = "IC8H18";
kname[4] = "NC10H22";
kname[5] = "NC7H16";
kname[6] = "NC12H26";
kname[7] = "O2";
kname[8] = "N2";
kname[9] = "NC16H34";
kname[10] = "H2O";
kname[11] = "CO2";
}

// compute the sparsity pattern of the chemistry Jacobian
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,169> Jac = {0.0};
amrex::GpuArray<amrex::Real,12> conc = {0.0};
for (int n=0; n<12; n++) {
    conc[n] = 1.0/ 12.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<13; k++) {
for (int l=0; l<13; l++) {
if(Jac[ 13 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}

*nJdata = NCELLS * nJdata_tmp;
}



// compute the sparsity pattern of the system Jacobian
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,169> Jac = {0.0};
amrex::GpuArray<amrex::Real,12> conc = {0.0};
for (int n=0; n<12; n++) {
    conc[n] = 1.0/ 12.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<13; k++) {
for (int l=0; l<13; l++) {
if(k == l){
nJdata_tmp = nJdata_tmp + 1;
} else {
if(Jac[ 13 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}
}

*nJdata = NCELLS * nJdata_tmp;
}



// compute the sparsity pattern of the simplified (for preconditioning) system Jacobian
void SPARSITY_INFO_SYST_SIMPLIFIED( int * nJdata, const int * consP)
{
amrex::GpuArray<amrex::Real,169> Jac = {0.0};
amrex::GpuArray<amrex::Real,12> conc = {0.0};
for (int n=0; n<12; n++) {
    conc[n] = 1.0/ 12.000000 ;
}
aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<13; k++) {
for (int l=0; l<13; l++) {
if(k == l){
nJdata_tmp = nJdata_tmp + 1;
} else {
if(Jac[ 13 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}
}

nJdata[0] = nJdata_tmp;
}


// compute the sparsity pattern of the chemistry Jacobian in CSC format -- base 0
void SPARSITY_PREPROC_CSC(int *  rowVals, int *  colPtrs, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,169> Jac = {0.0};
amrex::GpuArray<amrex::Real,12> conc = {0.0};
for (int n=0; n<12; n++) {
    conc[n] = 1.0/ 12.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

colPtrs[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset_row = nc * 13;
int offset_col = nc * 13;
for (int k=0; k<13; k++) {
for (int l=0; l<13; l++) {
if(Jac[13*k + l] != 0.0) {
rowVals[nJdata_tmp] = l + offset_row; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
colPtrs[offset_col + (k + 1)] = nJdata_tmp;
}
}
}

// compute the sparsity pattern of the chemistry Jacobian in CSR format -- base 0
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, const int * consP, int NCELLS, int base)
{
amrex::GpuArray<amrex::Real,169> Jac = {0.0};
amrex::GpuArray<amrex::Real,12> conc = {0.0};
for (int n=0; n<12; n++) {
    conc[n] = 1.0/ 12.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

if (base == 1) {
rowPtrs[0] = 1;
int nJdata_tmp = 1;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 13;
for (int l=0; l<13; l++) {
for (int k=0; k<13; k++) {
if(Jac[13*k + l] != 0.0) {
colVals[nJdata_tmp-1] = k+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
rowPtrs[offset + (l + 1)] = nJdata_tmp;
}
}
} else {
rowPtrs[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 13;
for (int l=0; l<13; l++) {
for (int k=0; k<13; k++) {
if(Jac[13*k + l] != 0.0) {
colVals[nJdata_tmp] = k + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
rowPtrs[offset + (l + 1)] = nJdata_tmp;
}
}
}
}

// compute the sparsity pattern of the system Jacobian
// CSR format BASE is user choice
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtr, const int * consP, int NCELLS, int base)
{
amrex::GpuArray<amrex::Real,169> Jac = {0.0};
amrex::GpuArray<amrex::Real,12> conc = {0.0};
for (int n=0; n<12; n++) {
    conc[n] = 1.0/ 12.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

if (base == 1) {
rowPtr[0] = 1;
int nJdata_tmp = 1;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 13;
for (int l=0; l<13; l++) {
for (int k=0; k<13; k++) {
if (k == l) {
colVals[nJdata_tmp-1] = l+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[13*k + l] != 0.0) {
colVals[nJdata_tmp-1] = k+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[offset + (l + 1)] = nJdata_tmp;
}
}
} else {
rowPtr[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 13;
for (int l=0; l<13; l++) {
for (int k=0; k<13; k++) {
if (k == l) {
colVals[nJdata_tmp] = l + offset; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[13*k + l] != 0.0) {
colVals[nJdata_tmp] = k + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[offset + (l + 1)] = nJdata_tmp;
}
}
}
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian on CPU
// BASE 0
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, const int * consP)
{
amrex::GpuArray<amrex::Real,169> Jac = {0.0};
amrex::GpuArray<amrex::Real,12> conc = {0.0};
for (int n=0; n<12; n++) {
    conc[n] = 1.0/ 12.000000 ;
}
aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

colPtrs[0] = 0;
int nJdata_tmp = 0;
for (int k=0; k<13; k++) {
for (int l=0; l<13; l++) {
if (k == l) {
rowVals[nJdata_tmp] = l; 
indx[nJdata_tmp] = 13*k + l;
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[13*k + l] != 0.0) {
rowVals[nJdata_tmp] = l; 
indx[nJdata_tmp] = 13*k + l;
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
colPtrs[k+1] = nJdata_tmp;
}
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// CSR format BASE is under choice
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, const int * consP, int base)
{
amrex::GpuArray<amrex::Real,169> Jac = {0.0};
amrex::GpuArray<amrex::Real,12> conc = {0.0};
for (int n=0; n<12; n++) {
    conc[n] = 1.0/ 12.000000 ;
}
aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

if (base == 1) {
rowPtr[0] = 1;
int nJdata_tmp = 1;
for (int l=0; l<13; l++) {
for (int k=0; k<13; k++) {
if (k == l) {
colVals[nJdata_tmp-1] = l+1; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[13*k + l] != 0.0) {
colVals[nJdata_tmp-1] = k+1; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[l+1] = nJdata_tmp;
}
} else {
rowPtr[0] = 0;
int nJdata_tmp = 0;
for (int l=0; l<13; l++) {
for (int k=0; k<13; k++) {
if (k == l) {
colVals[nJdata_tmp] = l; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[13*k + l] != 0.0) {
colVals[nJdata_tmp] = k; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[l+1] = nJdata_tmp;
}
}
}
