#include "config.h" // TODO

#include "CoSTAModule.hpp"

#include <pybind11/pybind11.h>

#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/TTagFlowProblemTPFA.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

namespace Opm::Properties {

template<class TypeTag>
struct Linearizer<TypeTag, TTag::FlowProblemTPFA> { using type = TpfaLinearizer<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FlowProblemTPFA> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowProblemTPFA> { static constexpr bool value = false; };

}

void export_BlackOil(pybind11::module& m)
{
    CoSTAModule<Opm::Properties::TTag::FlowProblemTPFA>::pyExport(m, "BlackOil");
}
