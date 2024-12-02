# pylibpvol/ase_calculator.py

from ase.calculators.calculator import Calculator, all_changes
from ase.units import Bohr,Hartree
from .calculator import LibpvolCalculator
import numpy as np

class LibpvolASECalculator(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, pressure: float, model=1, gridpts=1202, proberad=1.5,
                 verbose=False, printlevel=2, vdwSet=1, attached_calculator=None, **kwargs):
        """
        Initialize the libpvol ASE Calculator. \n
        One can attach another ASE calculator to it as the lib provides additive
        energies and gradients, rather than a full potential. \n
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
        Calculator.__init__(self, **kwargs)
        self.pressure = pressure
        self.model = model
        self.gridpts = gridpts
        self.proberad = proberad
        self.verbose = verbose
        self.printlevel = printlevel
        self.vdwSet = vdwSet
        self.calculator = None
        self.attached_calculator = attached_calculator

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        """
        Perform a calculation.
        """
        Calculator.calculate(self, atoms, properties, system_changes)


        # Extract information from the atoms object
        nat = len(atoms)
        at = atoms.get_atomic_numbers()
        xyz = atoms.get_positions() / Bohr  # pylibpvol expects Bohr!

        # use the attached calculator for the underlying potential
        if self.attached_calculator:
            self.ase_calculator.calculate(atoms, properties, system_changes)
            attached_energy = self.attached_calculator.get_potential_energy()
            attached_forces = self.attached_calculator.get_forces()
        else:
            attached_energy = 0.0
            attached_forces = 0.0


        # Initialize the custom calculator if it hasn’t been done yet
        if self.calculator is None:
            self.calculator = LibpvolCalculator(nat, at, xyz, self.pressure,
                                                 self.model, self.gridpts, self.proberad,
                                                 self.verbose, self.printlevel, self.vdwSet)
        
        # Perform the calculation
        energy, forces, iostat = self.calculator.calculate_single_point(nat, at, xyz)
        
        # Store the results in the ASE calculator's results dictionary
        # but make sure to do so in eV and eV/Ang units!
        self.results['energy'] = energy*Hartree + attached_energy
        self.results['forces'] = forces*Hartree/Bohr + attached_forces

    def clean_up(self):
        """Deallocate the calculator when done."""
        if self.calculator is not None:
            self.calculator.finalize()

