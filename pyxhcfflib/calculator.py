# pyxhcfflib/calculator.py

from .bindings import initialize_calculator

class XHCFFLibCalculator:
    """
    A high-level Python interface for the XHCFFLib calculator, 
    providing methods to perform calculations and manage resources.

    This class wraps the lower-level pybind11 bindings and offers
    a more Pythonic API for users.
    """

    def __init__(self, nat, at, xyz, pressure=1.0, model=1, gridpts=2030, proberad=1.2,
                 verbose=True, printlevel=2, vdwSet=0):
        """
        Initialize the XHCFFLib calculator with optional arguments.

        Parameters:
            nat (int): Number of atoms in the system.
            at (numpy.ndarray): Atomic numbers of the atoms (1D array).
            xyz (numpy.ndarray): Atomic coordinates (2D array with shape (nat, 3)) in [Bohr].
            pressure (float, optional): Pressure in the system in [GPa].
            model (int, optional): Model index specifying the calculation model.
                                   Options are 0=XHCFF, 1=PV
            gridpts (int, optional): Number of angular (Lebedev) grid points for calculations.
                                     Default is 2030.
            proberad (float, optional): Probe radius used in calculations in [Ang]. 
                                        Default is 1.4.
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
        energy, forces, iostat = self.calculator.singlepoint(nat, at, xyz)
        return energy, forces, iostat

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
        Deallocate resources used by the XHCFFLib calculator.
        """
        self.calculator.deallocate()

