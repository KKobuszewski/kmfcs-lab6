#ifndef __ROOT_FINDER_HPP__
#define __ROOT_FINDER_HPP__

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>


class RootFinder
{
private:
    const gsl_root_fsolver_type *_T;
    gsl_root_fsolver *_solver;
    gsl_function* _F;
    double _root;
    
public:
    int status;
    
    RootFinder(const gsl_root_fsolver_type *T = gsl_root_fsolver_brent);
    ~RootFinder();
    
    void set_function(gsl_function* func);
    int find_root(double* r, double x_lo, double x_hi, const double accurancy = 1e-08, const unsigned max_iter = 100);
};

RootFinder::RootFinder(const gsl_root_fsolver_type* T):_T(T)
{
    _solver = gsl_root_fsolver_alloc(_T);
}

RootFinder::~RootFinder()
{
    gsl_root_fsolver_free(_solver);
}

inline void RootFinder::set_function(gsl_function* func)
{
    _F = func;
}

inline int RootFinder::find_root(double* r, double x_lo, double x_hi, const double accurancy, const unsigned max_iter)
{
    status = gsl_root_fsolver_set(_solver, _F, x_lo, x_hi);
    if ( status != GSL_SUCCESS) return status;

    printf("Using %s method\n", gsl_root_fsolver_name(_solver));
    printf("%5s [%9s, %9s] %9s %10s\n", "iter", "lower", "upper", "root","err(est)");
    
    unsigned iter = 0;
    do
    {
        iter++;
        status = gsl_root_fsolver_iterate(_solver);
        _root  = gsl_root_fsolver_root(_solver);
        x_lo   = gsl_root_fsolver_x_lower(_solver);
        x_hi   = gsl_root_fsolver_x_upper(_solver);
        status = gsl_root_test_interval(x_lo, x_hi, 0, accurancy);
        
        if (status == GSL_SUCCESS) printf ("Converged:\n");
        printf ("%5d [%.7f, %.7f] %.7f %.7f\n", iter, x_lo, x_hi, _root, x_hi - x_lo);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    
    *r = _root;
    return status;
}












#endif