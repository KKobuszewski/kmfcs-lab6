#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
// #include <stdint.h>
#include <cmath>
#include <complex>
#include <vector>
// #include <algorithm>


// #include <Eigen/Core>
// #include <Eigen/Dense>
// #include <Eigen/Eigenvalues>


using namespace std::literals::complex_literals;


// ============================================= PROGRAM PARAMETERS =============================================
#define NX   16384/4

#define ZMAX ((double) 250.0)
#define ZA   ((double)  75.0)
#define ZB   ((double) 175.0)
#define Z0   ((double) 125.0)


// ============================================= NUMEROV ALGORITHM SETTINGS ======================================

#define HBAR2_OVER_2M0 ((double) (197.3269*197.3269) / (0.51099 * 20))  // hbar^2 / 2m_e
#define gamma_eff(x) (( (x > ZA) && (x < ZB) ) ? 0.92 : 1.00) 

#include <numerov.h>


// electric potential
#define ELECTRIC 0.1


inline double electric_potential(const double x, const double x0 = 125.0, const double eF = 0.0)
{
    return (x - x0)*eF;
}

inline double potential(const double x, const double V0 = -200, const double x1 = ZA, const double x2 = ZB)
{
    return electric_potential(x,Z0,ELECTRIC) + (( (x > x1) && (x < x2) ) ? V0 : 0.0);
}


// definition in numerov.h
inline double numerov_potential(double x, double* params)
{
    return potential(x);
}


// ===============================================  SHOOTING METHOD ==============================================1

// 
inline double diff_norm(std::complex<double>* psi_fwd, std::complex<double>* psi_bck, const unsigned index, const int sgn)
{
    return  std::real(psi_fwd[index]) - sgn*std::real(psi_bck[index]);
    //return  std::real(psi_fwd[index] - psi_bck[index]);
}


// buffer for wavefunctions
std::complex<double> psi_fwd[NX] = {0};
std::complex<double> psi_bck[NX] = {0};

std::complex<double> psi_fwd_0[2] = {0.0,0.00001};
std::complex<double> psi_bck_0[2] = {0.0,0.00001};


#include <RootFinder.hpp>


double diff_for_E( double E, void *params )
{
    double xmin = 0.0;
    double xmax = ZMAX;
    double dx   = (xmax - xmin)/((double) NX);
    
    int *p = (int *)params;
    
    numerov1D<NUMEROV_FORWARD>(psi_fwd, psi_fwd_0, E, NX, xmin, xmax);
    numerov1D<NUMEROV_BACKWARD>(psi_bck, psi_bck_0, E, NX, xmin, xmax);
    
    return diff_norm(psi_fwd, psi_bck, p[0], p[1]); 
}


inline void extract_wavefunction(std::complex<double>* wf, double trunc, const double dx, const int sgn)
{
    const unsigned index1 = floor(trunc*((double)NX));
    const unsigned index2 = NX - index1 - 1;
    
    unsigned it = 0;
    for (it; it < index1; it++)
        wf[it] = psi_fwd[it];
    
    for (it; it <= index2; it++)
        wf[it] = 0.5*(psi_fwd[it] + ((double)sgn)*(psi_bck[it]));
    
    for (it; it < NX; it++)
        wf[it] = psi_bck[it];
    
    normalize(wf, NX, dx);
}


/*
 * This function looks for roots
 */
int shooting_method( std::vector<std::complex<double>*>* wavefunctions, 
                      std::vector<double>* energies,
                      double Emin, double Emax, 
                      const unsigned it, const unsigned index_to_compare = 3*NX/4,
                      const double accurancy=1e-8)
{   
    const double xmin = 0.0;
    const double xmax = ZMAX;
    const double dx   = (xmax - xmin)/((double) NX);
    const int sgn = (it%2) ? -1 : 1;
    int int_params[2] = {(int)index_to_compare,sgn};
#ifdef DEBUG
    std::cout << "x = [" << xmin << "," << xmax << "]" << std::endl;
#endif
    
    // for gsl root solver
    gsl_function F;
    F.function = &diff_for_E;
    F.params = (void*) int_params;
    RootFinder rf;
    rf.set_function(&F);
    
    // find differences on egdes of the energy interval
#ifdef DEBUG
    std::cout << "NUMEROV FORWARD" << std::endl;
#endif
    numerov1D<NUMEROV_FORWARD>(psi_fwd, psi_fwd_0, Emin, NX, xmin, xmax);
#ifdef DEBUG
    std::cout << "NUMEROV BACKWARD" << std::endl;
#endif
    numerov1D<NUMEROV_BACKWARD>(psi_bck, psi_bck_0, Emin, NX, xmin, xmax);
    double diff_min = diff_norm(psi_fwd, psi_bck, index_to_compare,sgn);
    
    
    numerov1D<NUMEROV_FORWARD>(psi_fwd, psi_fwd_0, Emax, NX, xmin, xmax);
    numerov1D<NUMEROV_BACKWARD>(psi_bck, psi_bck_0, Emax, NX, xmin, xmax);
    double diff_max = diff_norm(psi_fwd, psi_bck, index_to_compare,sgn);
    
    // if differences have different signs it means root is inside
    if ( diff_min * diff_max < 0 )
    {
        std::cout << "[" << diff_min << "," << diff_max << "]\t\tPossible root in interval E=[" << Emin << "," << Emax << "] !" << std::endl;
    std::cout << "SIGN: " << sgn << std::endl;
        
#ifdef DEBUG
        std::ostringstream oss; oss << "data/well_fwd_E" << Emax << ".bin";
        FILE *file_fwd = fopen (oss.str().c_str(), "wb");
        fwrite(psi_fwd,sizeof(std::complex<double>),NX,file_fwd);
        
        std::ostringstream osb; osb << "data/well_bck_E" << Emax << ".bin";
        FILE *file_bck = fopen (osb.str().c_str(), "wb");
        fwrite(psi_bck,sizeof(std::complex<double>),NX,file_bck);
#endif
        // find root
        double root;
        rf.find_root(&root,Emin,Emax);
        
        // save energies
        energies->push_back(root);
        //std::cout << "Root for E=" << root << std::endl;
        
        // save wavefunction
        numerov1D<NUMEROV_FORWARD>(psi_fwd, psi_fwd_0, root, NX, xmin, xmax);
        extract_wavefunction(wavefunctions->at(it),0.75,dx,sgn);
        
        return 1;
    }
    else
    {
        //std::cout << "[" << diff_min << "," << diff_max << "]\t\tNo root found in interval E=[" << Emin << "," << Emax << "] !" << std::endl;
        return 0;
    }
    
    
    // metoda bisekcji przesiac kazdy przedzial bo osiagnac rzadana dokladnosc
    
    
    // porownac wartosc i pochodna w punkcie zszycia ...
    
    // zp = 3/4 (nie srodek studni...)
}




int main(int argc, char* argv[])
{
    // najnizsza wartosc wlasna
    // E0 = -197.37   (okolo)
    
    
    double xmin = 0.0;
    double xmax = ZMAX;
    double dx   = (xmax - xmin)/((double) NX);
    
    std::ofstream file; file.open("potential_well_eF.dat");
    for (unsigned ix=0; ix < NX; ix++)
        file << std::fixed << (xmin + ix*dx) << "\t" << numerov_potential(xmin + ix*dx) << "\t" << compute_k2(xmin + ix*dx,0) << std::endl;
        
    
    
    
    unsigned max_wavefunctions = 10;
    std::vector<std::complex<double>*> wavefunctions(max_wavefunctions);
    std::vector<double> energies;
    
    std::complex<double>* wf_buffer = new std::complex<double> [NX*max_wavefunctions];
    for (unsigned it=0; it < wavefunctions.size(); it++)
        wavefunctions[it] =  wf_buffer + it*NX;
    
    double dE = 6.0;
    //double E = -1.965604e+02;
    unsigned it=0;
    for (double E = -220; E <= -40.0; E += dE)
    {
        it += shooting_method(&wavefunctions,&energies,E,E+dE,it,3*NX/4);
    }
    
    std::cout << std::endl;
    FILE* output_wf = fopen( "wavefunctions_eF.bin", "w" );
    if (!output_wf) { fprintf(stderr,"Error opening file wavefunction.bin!\n"); exit(EXIT_FAILURE); }
    fwrite(wf_buffer,1,sizeof(std::complex<double>)*NX*energies.size(),output_wf);
    
    std::ofstream output;
    output.open("energies_eF.dat");
    for (it=0; it<energies.size(); it++)
    {
        std::cout << "Root for E="   << energies[it] << std::endl;
        output    << std::scientific << energies[it] << "\t";
    }
    output << std::endl;
    fclose(output_wf);
    
    delete wf_buffer;
    
    return EXIT_SUCCESS;
}