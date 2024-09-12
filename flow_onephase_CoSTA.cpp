#include "config.h" // TODO

#include "CoSTAModule.hpp"

#include <pybind11/pybind11.h>

#include <opm/simulators/flow/Main.hpp>
#include <opm/models/blackoil/blackoilonephaseindices.hh>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

namespace Opm::Properties {

namespace TTag {
struct FlowWaterOnlyProblem {
    using InheritsFrom = std::tuple<FlowProblem>;
};
}

template<class TypeTag>
struct Linearizer<TypeTag, TTag::FlowWaterOnlyProblem> { using type = TpfaLinearizer<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FlowWaterOnlyProblem> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowWaterOnlyProblem> { static constexpr bool value = false; };

//! The indices required by the model
template<class TypeTag>
struct Indices<TypeTag, TTag::FlowWaterOnlyProblem>
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because this leads
    // to cyclic definitions of some properties. if this happens the compiler error
    // messages unfortunately are *really* confusing and not really helpful.
    using BaseTypeTag = TTag::FlowProblem;
    using FluidSystem = GetPropType<BaseTypeTag, Properties::FluidSystem>;

public:
    using type = BlackOilOnePhaseIndices<getPropValue<TypeTag, Properties::EnableSolvent>(),
                                         getPropValue<TypeTag, Properties::EnableExtbo>(),
                                         getPropValue<TypeTag, Properties::EnablePolymer>(),
                                         getPropValue<TypeTag, Properties::EnableEnergy>(),
                                         getPropValue<TypeTag, Properties::EnableFoam>(),
                                         getPropValue<TypeTag, Properties::EnableBrine>(),
                                         /*PVOffset=*/0,
                                         /*enabledCompIdx=*/FluidSystem::waterCompIdx,
                                         getPropValue<TypeTag, Properties::EnableMICP>()>;
};

} // namespace Opm::Properties

void export_OnePhase(pybind11::module& m)
{
    CoSTAModule<Opm::Properties::TTag::FlowWaterOnlyProblem>::pyExport(m, "OnePhase");
}
