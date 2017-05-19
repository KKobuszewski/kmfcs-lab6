#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <stdint.h>







// to be declared somewhere else
inline double numerov_potential(double x, double* params = NULL);



inline double normalize(std::complex<double> *psi, const size_t Nx, const double dx)
{
    double norm = 0.0;
    
    #pragma omp simd
    for (size_t ix=0; ix < Nx; ix++)
    {
        norm += psi[ix].real()*psi[ix].real() + psi[ix].imag()*psi[ix].imag();
    }
    
    norm = sqrt(norm*dx);
    #pragma omp simd
    for (size_t ix=0; ix < Nx; ix++)
    {
        psi[ix] /= norm;
    }
    
    return norm*norm;
}



inline double compute_k2(const double x, const double E,double* params = NULL)
{
    return  (E - numerov_potential(x,params))/HBAR2_OVER_2M0*gamma_eff(x);
}


#define NUMEROV_FORWARD  ( 1)
#define NUMEROV_BACKWARD (-1)

template <int8_t direction> // direction = {-1,1}
void numerov1D(std::complex<double> *psi, std::complex<double> *psi_0, const double E, const size_t Nx, const double xmin, const double xmax)
{
    const double dx = (xmax - xmin)/((double) Nx);
    
    switch(direction)
    {
        case NUMEROV_FORWARD:
        {
            psi[0] = psi_0[0];
            psi[1] = psi_0[1];
            double x = xmin + 2.0*dx;
#ifdef DEBUG
    std::cout << std::fixed << xmin    << "\t" << psi[0].real() << " + " << psi[0].imag() << "j" << std::endl;
    std::cout << std::fixed << xmin+dx << "\t" << psi[1].real() << " + " << psi[1].imag() << "j" << std::endl;
#endif
            //std::cout << "dx: " << dx <<std::endl;
            //std::cout << " x: " <<  x <<std::endl;
            
            for (size_t ix=2; ix < Nx; ix++)
            {
                psi[ix] = ( 2.0*(1.0 - 5.0/12.0 * dx*dx * compute_k2(x -     dx,E)) * psi[ix-1]  - 
                                (1.0 + 1.0/12.0 * dx*dx * compute_k2(x - 2.0*dx,E)) * psi[ix-2] ) /
                                (1.0 + 1.0/12.0 * dx*dx * compute_k2(x         ,E));
#ifdef DEBUG
                 std::cout << std::fixed << x << "\t" << psi[ix].real() << " + " << psi[ix].imag() << "j\t" << compute_k2(x,E) << std::endl;
#endif
                x += dx;
            }
            
    //std::cout << psi[2].real() << " + " << psi[2].imag() << "j" << std::endl;
    //std::cout << psi[Nx-1].real() << " + " << psi[Nx-1].imag() << "j" << std::endl;
            
            break;
        }
        
        case NUMEROV_BACKWARD:
        {
            psi[Nx-1] = psi_0[0];
            psi[Nx-2] = psi_0[1];
            double x = xmax - 2.0*dx;
#ifdef DEBUG
    std::cout << std::fixed << xmax    << "\t" << psi[0].real() << " + " << psi[0].imag() << "j" << std::endl;
    std::cout << std::fixed << xmax-dx << "\t" << psi[1].real() << " + " << psi[1].imag() << "j" << std::endl;
#endif
            
            for (size_t ix=Nx-3; x >= 0.0; ix--)
            {
                psi[ix] = ( 2.0*(1.0 - 5.0/12.0 * dx*dx * compute_k2(x +     dx,E)) * psi[ix+1]  - 
                                (1.0 + 1.0/12.0 * dx*dx * compute_k2(x + 2.0*dx,E)) * psi[ix+2] ) /
                                (1.0 + 1.0/12.0 * dx*dx * compute_k2(x         ,E));
#ifdef DEBUG
                 std::cout << std::fixed << x << "\t" << psi[ix].real() << " + " << psi[ix].imag() << "j\t" << compute_k2(x,E) << std::endl;
#endif
                x -= dx;
            }
            break;
        }
    }
    
    normalize(psi,Nx,dx);
}
