#include <cstdlib>
#include <iostream>

#include <RootFinder.hpp>
#include "demo_fn.h"


double nonpolynomial(double x, void *params)
{
  struct quadratic_params *p 
    = (struct quadratic_params *) params;

  double a = p->a;
  double b = p->b;
  double c = p->c;

  return (a * x + b) * sin(x) + c;
}


int main()
{
    gsl_function F;
    struct quadratic_params params = {1.0, 0.0, -5.0};
    F.function = &quadratic;
    F.params = &params;
    
    double root = 0.0;
    int status;
    
    RootFinder rf;
    rf.set_function(&F);
    status = rf.find_root(&root,1.0,5.0);
    std::cout << "GSL status: " << std::string( gsl_strerror(status) ) << std::endl;
    
    std::cout << std::endl;
    std::cout << "Second root" << std::endl;
    
    
    rf.set_function(&F);
    status = rf.find_root(&root,-5.0,1.0);
    std::cout << "GSL status: " << std::string( gsl_strerror(status) ) << std::endl;
    
    std::cout << std::endl;
    std::cout << "Another function" << std::endl;
    
    F.function = &nonpolynomial;
    rf.set_function(&F);
    status = rf.find_root(&root,-10.0,10.0);
    std::cout << "GSL status: " << std::string( gsl_strerror(status) ) << std::endl;
    
    return EXIT_SUCCESS;
}