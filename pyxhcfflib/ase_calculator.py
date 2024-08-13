# pyxhcfflib/ase_calculator.py

from ase.calculators.calculator import Calculator, all_changes
from .calculator import XHCFFLibCalculator
import numpy as np

class XHCFFLibASECalculator(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, pressure=1.0, model=1, gridpts=2030, proberad=1.4,
                 verbose=True, printlevel=2, vdwSet=0, **kwargs):
        """
        Initialize the XHCFFLib ASE Calculator.
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

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        """
        Perform a calculation.
        """
        Calculator.calculate(self, atoms, properties, system_changes)

        # Extract information from the atoms object
        nat = len(atoms)
        at = atoms.get_atomic_numbers()
        xyz = atoms.get_positions()

        # Initialize the custom calculator if it hasn’t been done yet
        if self.calculator is None:
            self.calculator = XHCFFLibCalculator(nat, at, xyz, self.pressure,
                                                 self.model, self.gridpts, self.proberad,
                                                 self.verbose, self.printlevel, self.vdwSet)
        
        # Perform the calculation
        energy, forces, iostat = self.calculator.calculate_single_point(nat, at, xyz)
        
        # Store the results in the ASE calculator's results dictionary
        self.results['energy'] = energy
        self.results['forces'] = forces

    def clean_up(self):
        """Deallocate the calculator when done."""
        if self.calculator is not None:
            self.calculator.finalize()

