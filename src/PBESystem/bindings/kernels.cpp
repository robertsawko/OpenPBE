#include <boost/python.hpp>
#include "../breakupKernels/CoulaloglouTavlarides.H"
#include "../breakupKernels/binaryBreakup.H"
#include "../breakupKernels/noBreakup.H"
#include "../coalescenceKernels/CoulaloglouTavlarides.H"
#include "../coalescenceKernels/constantCoalescence.H"
#include "../coalescenceKernels/noCoalescence.H"

BOOST_PYTHON_MODULE(kernels)
{
    using namespace boost::python;
    using Foam::scalar;

    class_<Foam::breakupKernels::CoulaloglouTavlaridesImp>(
                "CoulaloglouTavlaridesBreakup",
                init<scalar, scalar, scalar, scalar>())
            .def("S", &Foam::breakupKernels::CoulaloglouTavlaridesImp::S);

    class_<Foam::breakupKernels::binaryBreakupImpl>(
                "BinaryBreakup")
            .def("S", &Foam::breakupKernels::binaryBreakupImpl::S);

    class_<Foam::breakupKernels::noBreakupImpl>(
                "NoBreakup")
            .def("S", &Foam::breakupKernels::noBreakupImpl::S);

    class_<Foam::coalescenceKernels::CoulaloglouTavlaridesCImpl>(
                "CoulaloglouTavlaridesCoalescence",
                init<scalar, scalar, scalar, scalar>())
            .def("S", &Foam::coalescenceKernels::CoulaloglouTavlaridesCImpl::S);

    class_<Foam::coalescenceKernels::constantCoalescenceImpl>(
                "ConstantCoalescence",
                init<scalar>())
            .def("S", &Foam::coalescenceKernels::constantCoalescenceImpl::S);

    class_<Foam::coalescenceKernels::noCoalescenceImpl>(
                "NoCoalescence")
            .def("S", &Foam::coalescenceKernels::noCoalescenceImpl::S);
}
