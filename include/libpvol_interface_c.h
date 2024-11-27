#ifndef LIBPVOL_INTERFACE_C_H
#define LIBPVOL_INTERFACE_C_H

#include <stdbool.h> // For booleans
#include <stdio.h>   // For FILE*
#include <stdlib.h>  // For standard library functions


#ifdef __cplusplus
extern "C" {
#endif

// Declare the Fortran structure
typedef struct {
  void *ptr;
} c_libpvol_calculator;

// Declate the initializer
extern c_libpvol_calculator
c_libpvol_calculator_init(int nat, int *at, double (*xyz)[3], double pressure,
                           int model, int gridpts, double proberad,
                           bool verbose, int printlevel, int vdwSet);

// Declare the deallocator
extern void c_libpvol_calculator_deallocate(c_libpvol_calculator *calculator);

// Declate the singlepoint calculator
extern void c_libpvol_calculator_singlepoint(c_libpvol_calculator *calculator,
                                              int nat, int *at,
                                              double (*xyz)[3], double *energy,
                                              double (*gradient)[3],
                                              int *iostat);

// Declate the print routine
extern void c_libpvol_calculator_info(c_libpvol_calculator *calculator,
                                       int iunit);


#ifdef __cplusplus
}
#endif

#endif /* LIBPVOL_INTERFACE_C_H */
