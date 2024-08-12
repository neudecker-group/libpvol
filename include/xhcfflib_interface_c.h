#ifndef XHCFFLIB_INTERFACE_C_H
#define XHCFFLIB_INTERFACE_C_H

#include <stdbool.h> // For booleans
#include <stdio.h>   // For FILE*
#include <stdlib.h>  // For standard library functions


#ifdef __cplusplus
extern "C" {
#endif

// Declare the Fortran structure
typedef struct {
  void *ptr;
} c_xhcfflib_calculator;

// Declate the initializer
extern c_xhcfflib_calculator
c_xhcfflib_calculator_init(int nat, int *at, double (*xyz)[3], double pressure,
                           int model, int gridpts, double proberad,
                           bool verbose, int printlevel, int vdwSet);

// Declare the deallocator
extern void c_xhcfflib_calculator_deallocate(c_xhcfflib_calculator *calculator);

// Declate the singlepoint calculator
extern void c_xhcfflib_calculator_singlepoint(c_xhcfflib_calculator *calculator,
                                              int nat, int *at,
                                              double (*xyz)[3], double *energy,
                                              double (*gradient)[3],
                                              int *iostat);

// Declate the print routine
extern void c_xhcfflib_calculator_info(c_xhcfflib_calculator *calculator,
                                       int iunit);


#ifdef __cplusplus
}
#endif

#endif /* XHCFFLIB_INTERFACE_C_H */
