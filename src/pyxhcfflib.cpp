#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "xhcfflib_interface_c.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

// a minimal test function
int add(int i, int j) {
    return i + j;
}

// Wrapper for the c_xhcfflib_calculator_init function
c_xhcfflib_calculator calculator_init(int nat, py::array_t<int> at, py::array_t<double> xyz, double pressure,
                                      int model, int gridpts, double proberad,
                                      bool verbose, int printlevel, int vdwSet) {
    // Ensure the arrays are contiguous
    auto at_buf = at.request();
    auto xyz_buf = xyz.request();

    if (at_buf.ndim != 1 || xyz_buf.ndim != 2 || xyz_buf.shape[1] != 3) {
        throw std::runtime_error("Invalid array dimensions");
    }

    int *at_ptr = static_cast<int *>(at_buf.ptr);
    double (*xyz_ptr)[3] = reinterpret_cast<double (*)[3]>(xyz_buf.ptr);

    return c_xhcfflib_calculator_init(nat, at_ptr, xyz_ptr, pressure, model, gridpts, proberad, verbose, printlevel, vdwSet);
}

// Wrapper for the c_xhcfflib_calculator_deallocate function
void calculator_deallocate(c_xhcfflib_calculator *calculator) {
    c_xhcfflib_calculator_deallocate(calculator);
}

// Wrapper for the c_xhcfflib_calculator_singlepoint function
void calculator_singlepoint(c_xhcfflib_calculator *calculator,
                            int nat, py::array_t<int> at,
                            py::array_t<double> xyz, py::array_t<double> energy,
                            py::array_t<double> gradient, py::array_t<int> iostat) {
    // Ensure the arrays are contiguous
    auto at_buf = at.request();
    auto xyz_buf = xyz.request();
    auto energy_buf = energy.request();
    auto gradient_buf = gradient.request();
    auto iostat_buf = iostat.request();

    if (at_buf.ndim != 1 || xyz_buf.ndim != 2 || xyz_buf.shape[1] != 3 ||
        energy_buf.ndim != 1 || gradient_buf.ndim != 2 || gradient_buf.shape[1] != 3 ||
        iostat_buf.ndim != 1) {
        throw std::runtime_error("Invalid array dimensions");
    }

    int *at_ptr = static_cast<int *>(at_buf.ptr);
    double (*xyz_ptr)[3] = reinterpret_cast<double (*)[3]>(xyz_buf.ptr);
    double *energy_ptr = static_cast<double *>(energy_buf.ptr);
    double (*gradient_ptr)[3] = reinterpret_cast<double (*)[3]>(gradient_buf.ptr);
    int *iostat_ptr = static_cast<int *>(iostat_buf.ptr);

    c_xhcfflib_calculator_singlepoint(calculator, nat, at_ptr, xyz_ptr, energy_ptr, gradient_ptr, iostat_ptr);
}

// Wrapper for the c_xhcfflib_calculator_info function
void calculator_info(c_xhcfflib_calculator *calculator, int iunit) {
    c_xhcfflib_calculator_info(calculator, iunit);
}

PYBIND11_MODULE(_xhcfflib, m) {
    m.doc() = "Python interface for xhcfflib C/Fortran code";

    py::class_<c_xhcfflib_calculator>(m, "Calculator")
        .def(py::init(&calculator_init), "Initialize the calculator")
        .def("deallocate", &calculator_deallocate, "Deallocate the calculator")
        .def("singlepoint", &calculator_singlepoint, "Perform single-point calculation")
        .def("info", &calculator_info, "Print calculator information");

    m.def("add", &add, "small test function");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}

