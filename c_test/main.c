#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

// Declare the Fortran structure
typedef struct {
  void *ptr;
} c_xhcfflib_calculator;

// Declate the initializer
extern c_xhcfflib_calculator
c_xhcfflib_calculator_init(int nat, int *at, double xyz[3][24], double pressure,
                           int model, int gridpts, double proberad,
                           bool verbose, int printlevel, int vdwSet);

// Declare the deallocator
extern void c_xhcfflib_calculator_deallocate(c_xhcfflib_calculator *calculator);

// Declate the singlepoint calculator
extern void c_xhcfflib_calculator_singlepoint(c_xhcfflib_calculator *calculator,
                                              int nat, int *at,
                                              double xyz[3][24], double *energy,
                                              double gradient[3][24],
                                              int *iostat);

// Declate the print routine
extern void c_xhcfflib_calculator_info(c_xhcfflib_calculator* calculator, int iunit);

// The test implementation program
int main() {
  // Test molecule: caffeine
  int nat = 24; // Number of atoms
  int at[24] = {6, 7, 6, 7, 6, 6, 6, 8, 7, 6, 8, 7,
                6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // Atom types
  double xyz[3][24] = {
      {2.02799738646442,  4.75011007621000,  6.33434307654413,
       8.72860718071825,  8.65318821103610,  6.23857175648671,
       5.63266886875962,  3.44931709749015,  7.77508917214346,
       10.30229550927022, 12.07410272485492, 10.70038521493902,
       13.24597858727017, 7.40891694074004,  1.38702118184179,
       1.34622199478497,  1.34624089204623,  5.65596919189118,
       14.67430918222276, 13.50897177220290, 13.50780014200488,
       5.41408424778406,  8.31919801555568,  8.31511620712388},
      {0.09231312124713,  0.02373496014051,  2.07098865582721,
       1.38002919517619,  -1.19324866489847, -2.08353643730276,
       -4.69950321056008, -5.48092386085491, -6.24427872938674,
       -5.39739796609292, -6.91573621641911, -2.79078533715849,
       -1.76969072232377, -8.95905928176407, 2.05575746325296,
       -0.86356704498496, -0.86133716815647, 4.00172183859480,
       -3.26230980007732, -0.60815166181684, -0.60614855212345,
       -9.49239668625902, -9.74947502841788, -9.76854236502758},
      {-0.14310895950963, -0.14324124033844, -0.14235306905930,
       -0.14265542523943, -0.14231527453678, -0.14218299370797,
       -0.13940509630299, -0.14318454855466, -0.13107140408805,
       -0.13672168520430, -0.13666499342053, -0.14148379504141,
       -0.14218299370797, -0.11636933482904, -0.14178615122154,
       1.55590600570783,  -1.84340893849267, -0.14131371969009,
       -0.14344911021228, 1.54898960808727,  -1.83214617078268,
       -0.11022772492007, 1.56539243085954,  -1.79108242206824}};

  // Test variables
  double pressure = 1.0;
  int model = 1; // 0 for XHCFF, 1 for PV
  int gridpts = 2030;
  double proberad = 1.4;
  bool verbose = false; // we test the info routine separately below
  int vdwSet = 0;
  int printlevel=2;

  // Call the Fortran function
  c_xhcfflib_calculator calc = c_xhcfflib_calculator_init(
      nat, at, xyz, pressure, model, gridpts, proberad, verbose, printlevel, vdwSet);

  if (calc.ptr == NULL) {
    printf("Error initializing xhcfflib calculator.\n");
    return 1;
  }

  printf("xhcfflib calculator initialized successfully.\n");

  // Use the calculator...
  double energy;
  double gradient[3][nat]; // Adjust the size to match nat
  int iostat;

  // Call the singlepoint function
  c_xhcfflib_calculator_singlepoint(&calc, nat, at, xyz, &energy, gradient,
                                    &iostat);

  // Check the result and print it
  if (iostat == 0) {
    printf("Singlepoint calculation successful.\n");
    printf("Energy: %f\n", energy);

    // Print the gradient for the first atoms(optional)
    for (int i = 0; i < 3; i++) {
      int j=0 ;
      printf("Gradient[%d][%d] = %e\n", i, j, gradient[i][j]);
    }
  } else {
    printf("Singlepoint calculation failed with iostat = %d\n", iostat);
  }

  // Also test the intrinsic print routine
  int iunit = 6;  // Use 6 to get Fortran's STDOUT

  // Call the print function
  c_xhcfflib_calculator_info(&calc, iunit);


  // Deallocate the Fortran object
  c_xhcfflib_calculator_deallocate(&calc);

  if (calc.ptr == NULL) {
    printf("xhcfflib calculator deallocated successfully.\n");
  } else {
    printf("Error deallocating xhcfflib calculator.\n");
  }

  return 0;
}
