# pyxhcfflib/ase_calculator.py

from ase.calculators.calculator import Calculator, all_changes
from ase.units import Bohr,Hartree
from .calculator import XHCFFLibCalculator
import numpy as np

class XHCFFLibASECalculator(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, pressure=1.0, model=1, gridpts=2030, proberad=1.2,
                 verbose=False, printlevel=2, vdwSet=1, attached_calculator=None, **kwargs):
        """
        Initialize the XHCFFLib ASE Calculator.
        One can attach another ASE calculator to it as the lib provides additive
        energies and gradients, rather than a full potential.
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
        xyz = atoms.get_positions() / Bohr  # pyxhcfflib expects Bohr!

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
            self.calculator = XHCFFLibCalculator(nat, at, xyz, self.pressure,
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

