from skbuild import setup
from setuptools import find_packages

# Define CMake arguments
cmake_args = [
    '-DPYTHON_BINDINGS=ON',  # Pass the custom option to CMake
]

setup(
    name='pyxhcfflib',
    version='0.0.1',
    author='Author1, Author2, Author3',
    author_email='email1@example.com, email2@example.com',
    description='A Python interface to the xhcfflib Fortran/C++ code',
    long_description=open('README.md').read(),
    packages=find_packages(),
    cmake_install_dir='_build_pyxhcfflib',
    include_package_data=True,
    zip_safe=False,
    cmake_args=cmake_args,  # Pass the CMake arguments here
)

