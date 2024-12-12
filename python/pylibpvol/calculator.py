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

    def __init__(self, nat: int, at: np.ndarray, xyz: np.ndarray, pressure: float,
                  model=1, gridpts=1202, proberad=1.5,verbose=False, printlevel=2, vdwSet=0):
        """
        Initialize the libpvol calculator with optional arguments. \n
            :param int nat: Number of atoms in the system.
            :param numpy.ndarray at: Atomic numbers of the atoms (1D array).
            :param numpy.ndarray xyz: Atomic coordinates (2D array with shape (nat, 3)) in [Bohr].
            pressure (float): Pressure in the system in [GPa].
            :param float pressure: Pressure in the system in [GPa].
            :param int, optional model: Model index specifying the calculation model. \n
                                Options: 0=XHCFF, 1=PV.
            :param int, optional gridpts: Number of angular (Lebedev) grid points for calculations. \n
                                    Default is 1202.
            :param float, otional proberad: Probe radius used in calculations in [Ang].\n
                                        Default is 1.5.
            :param bool, optional verbose: If True, provides detailed output *during the setup*
            :param int, optional printlevel: Level of detail for printed output. \n
                                            Default is 2.
            :param int, optional vdwSet: Index specifying the Van der Waals set to use. \n
                                        Options: 0=D3 (H-Pu), 1=Bondi (H-Ar)  \n
                                        Default is 0.
        """
        self.calculator = initialize_calculator(nat, at, xyz, pressure,
                                                model, gridpts, proberad,
                                                verbose, printlevel, vdwSet)

    def calculate_single_point(self, nat: int, at: np.ndarray, xyz: np.ndarray):
        """
        Perform a single-point energy and force calculation. \n
            :param int nat: Number of atoms in the system.
            :param numpy.ndarray at: Atomic numbers of the atoms (1D array).
            :param numpy.ndarray xyz: Atomic coordinates (2D array with shape (nat, 3)).
            :returns: energy -- Calculated energy of the system in [Hartree]. \n
                      forces -- Calculated forces on each atom (shape (nat, 3)) in [Hartree/Bohr].\n
                      iostat -- I/O status indicating success (0) or failure (non-zero).
            :rtype: (energy: int, forces: np.ndarray, iostat: int)
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
        :param int iunit: Output unit number where the information will be printed.\n
          Default is 6, which is stdout in Fortran
        """
        self.calculator.info(iunit)

    def finalize(self):
        """
        Deallocate resources used by the libpvol calculator.
        """
        self.calculator.deallocate()

