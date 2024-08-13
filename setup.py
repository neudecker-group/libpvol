from skbuild import setup
from setuptools import find_packages

# Define CMake arguments
cmake_args = [
    '-DPYTHON_BINDINGS=ON',  # Pass the custom option to CMake
]

setup(
    name='pyxhcfflib',
    version='0.0.1',
    author='Felix Zeller, Philipp Pracht, Tim Neudecker',
    author_email='zellerf@uni-bremen.de, research@philipp-pracht.de',
    description='A Python interface to the xhcfflib Fortran/C++ code',
    long_description=open('README.md').read(),
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    cmake_args=cmake_args,  # Pass the CMake arguments here
    install_requires=[
        'ase>=3.22.0',  # Ensure ASE is installed with a version that meets your needs
        # Add other dependencies if necessary
    ],
)

