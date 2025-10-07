#include "PySTS.h"

namespace py = pybind11;

PYBIND11_MODULE(PySTS, m)
{
    m.doc() = "ScrewTheorySolvers Python wrapper"; // optional module docstring
    m.attr("__version__") = "1.0.0";
    init_solvers(m);
}
