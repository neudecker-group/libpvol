
# PVol standalone library

This repository contains a standalone implementation of
the PV and XHCFF methods (https://doi.org/10.1063/5.0024671). 

The project provides a statically linked library (`libpvol.a`) with Fortran and C++ interfaces that can be linked in other projects. A Python implementation via the C++ interface is also provided.

`test_pvol.F90` and `main.cpp` in `test/` demonstrate the in-code usage.


## Building the Project

Make sure you have the following dependencies installed:

- CMake and `make` (or meson and ninja build systems)
- Fortran and C compilers (e.g., `gfortran`/`gcc` or `ifort`/`icc`)


### Installation

While some portion the program backend is written in Fortran, the Python and C++ bindings enable to  install the package via `pip`, which will automatically compile the code. Use the `-v` option to     keep track of the compilation process:
```bash
git clone git@github.com:neudecker-group/libpvol.git
cd libpvol
pip install . -v
```


<details>
<summary><h4>Fortran build instructions (without `pip`)</h4></summary>

Follow these steps to build the project:

1. Create a build directory and navigate to it
   ```bash
   mkdir _build
   cd _build
   ```

2. Export the compilers (here for example `ifort`/`icc`) and depending on your chosen build system set up the build:
   - generate the build files using CMake:
     ```bash
     FC=ifort CC=icc cmake ..
     ```
   - generate the build files using meson:
     ```bash
     FC=ifort CC=icc meson ..
     ```
   I you wish to build the test-binary, add `-Dbuild_exe=true` to either the `cmake` or `meson` setup command.


3. Depending on your chosen build system, build the project. If you have multiple cores/processors, you can speed up the build process by specifying the number of cores to use with the `-j` option. For example, to use 4 cores:
   - With CMake/`make`:
     ```shell
     make -j4
     ```
   - With meson/`ninja`:
     ```shell
     ninja -j4
     ```

### Cleaning the Build

To clean the build files, simply delete the `build` directory:

```shell
rm -rf _build
```
</details>



