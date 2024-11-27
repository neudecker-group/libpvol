# pylibpvol/calculator.py

import numpy as np
from .bindings import initialize_calculator

class LibpvolCalculator:
    """
    A high-level Python interface for the libpvol calculator, 
    providing methods to perform calculations and manage resources.

    This class wraps the lower-level pybind11 bindings and offers
    a more Pythonic API for users.
    """

    def __init__(self, nat, at, xyz, pressure=1.0, model=1, gridpts=1202, proberad=1.5,
                 verbose=True, printlevel=2, vdwSet=0):
        """
        Initialize the libpvol calculator with optional arguments.

        Parameters:
            nat (int): Number of atoms in the system.
            at (numpy.ndarray): Atomic numbers of the atoms (1D array).
            xyz (numpy.ndarray): Atomic coordinates (2D array with shape (nat, 3)) in [Bohr].
            pressure (float, optional): Pressure in the system in [GPa].
            model (int, optional): Model index specifying the calculation model.
                                   Options are 0=XHCFF, 1=PV
            gridpts (int, optional): Number of angular (Lebedev) grid points for calculations.
                                     Default is 1202.
            proberad (float, optional): Probe radius used in calculations in [Ang]. 
                                        Default is 1.5.
            verbose (bool, optional): If True, provides detailed output *during the setup*
            printlevel (int, optional): Level of detail for printed output. Default is 2.
            vdwSet (int, optional): Index specifying the Van der Waals set to use. Default is 0.
        """
        self.calculator = initialize_calculator(nat, at, xyz, pressure,
                                                model, gridpts, proberad,
                                                verbose, printlevel, vdwSet)

    def calculate_single_point(self, nat, at, xyz):
        """
        Perform a single-point energy and force calculation.

        Parameters:
            nat (int): Number of atoms in the system.
            at (numpy.ndarray): Atomic numbers of the atoms (1D array).
            xyz (numpy.ndarray): Atomic coordinates (2D array with shape (nat, 3)).

        Returns:
           energy (float): Calculated energy of the system in [Hartree].
           forces (numpy.ndarray): Calculated forces on each atom (shape (nat, 3)) in [Hartree/Bohr].
           iostat (int): I/O status indicating success (0) or failure (non-zero).
        """
        energy = np.zeros(1, dtype=np.float64)  # Assuming energy is a scalar
        forces = np.zeros((nat, 3), dtype=np.float64)  # Force array of size (nat, 3)
        iostat = np.zeros(1, dtype=np.int32)  # Assuming iostat is a scalar or single value

        # Call the C++ function through the Python binding
        self.calculator.singlepoint(nat, at, xyz, energy, forces, iostat)
        
        # Return the results
        return energy[0], forces, iostat[0]


    def print_info(self, iunit=6):
        """
        Print detailed information about the calculator's internal state.

        Parameters:
            iunit (int, optional): Output unit number where the information
                                   will be printed. Default is 6, which is stdout in Fortran
        """
        self.calculator.info(iunit)

    def finalize(self):
        """
        Deallocate resources used by the libpvol calculator.
        """
        self.calculator.deallocate()

