#include <boost/python.hpp>
#include "../breakupKernels/CoulaloglouTavlarides.H"

BOOST_PYTHON_MODULE(kernels)
{
    using namespace boost::python;
    using Foam::scalar;

    class_<Foam::breakupKernels::CoulaloglouTavlaridesImp>("CoulaloglouTavlaridesImp"
                                     , init<scalar, scalar, scalar, scalar>())
            .def("S", &Foam::breakupKernels::CoulaloglouTavlaridesImp::S);
}
