#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <stdint.h>

// numerov algorithm settings and include
#define ZA   (-100.0)
#define ZB   (100.0)
#define HBAR2_OVER_2M0 0.5//((double) 197.3269*197.3269 / (0.51099 * 2))  // hbar^2 / 2m_e
#define gamma_eff(x)  1.0 
#include <numerov.h>


#define NX 100000

using namespace std::literals::complex_literals;


// definition in numerov.h
inline double numerov_potential(double x, double* params)
{
    return 0.5*x*x;
}

inline std::complex<double> harmonic_groundstate(const double x)
{
    return  pow( 1./M_PI , 0.25) * exp ( - 0.5* x*x) + 1i*0.0;
}

int main()
{
    double xmin = -5.0;
    double xmax =  5.0;
    double dx   = (xmax - xmin)/((double) NX);
    
    std::complex<double>* psi = new std::complex<double> [NX];
    std::complex<double> psi_0[2] = {harmonic_groundstate(xmin), harmonic_groundstate(xmin+dx)};
    
    double energy;
    for (energy = 0.49; energy <= 0.51; energy += 0.01)
    {
        std::cout << "Energy: " << energy << std::endl;
        numerov1D<NUMEROV_FORWARD>(psi, psi_0, energy, NX, xmin, xmax);
        
        //psi_0[0] = harmonic_groundstate(xmax);
        //psi_0[1] = harmonic_groundstate(xmax-dx);
        //numerov1D<NUMEROV_BACKWARD>(psi+NX/2, psi_0, 0.5, NX/2, 0.0, xmax);
        //::cout << psi[0].real() << " + " << psi[0].imag() << "j" << std::endl;
        //std::cout << psi[1].real() << " + " << psi[1].imag() << "j" << std::endl;
        
        std::ostringstream oss;
        oss << "harmonic_E" << energy << ".dat";
        std::cout << oss.str() << std::endl;
        
        std::ofstream file;
        file.open(oss.str());
        
        for (unsigned ix=0; ix < NX; ix++)
        {
            file      << std::setprecision(15) << std::fixed;
            file      << xmin + ix*dx << "\t" << psi[ix].real() << "\t" << psi[ix].imag() << "\t" << harmonic_groundstate(xmin + ix*dx).real() << std::endl;
        }
        
        file.close();
    }
    
    std::cout << std::endl;
    std::cout << "Numerov backwardd scheme" << std::endl;
    
    
    psi_0[0] = harmonic_groundstate(xmax);
    psi_0[1] = harmonic_groundstate(xmax-dx);
    numerov1D<NUMEROV_BACKWARD>(psi, psi_0, 0.5, NX, xmin, xmax);
    
    std::ofstream file;
    file.open("harmonic_backward_E0.5.dat");
    
    for (unsigned ix=0; ix < NX; ix++)
    {
        file      << std::setprecision(15) << std::fixed;
        file      << xmax - ix*dx << "\t" << psi[ix].real() << "\t" << psi[ix].imag() << "\t" << harmonic_groundstate(xmax - ix*dx).real() << std::endl;
    }
    
    file.close();
    
    delete psi;
    
    return EXIT_SUCCESS;
}