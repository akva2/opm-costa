// $Id$
//==============================================================================
//!
//! \file opm_CoSTA.cpp
//!
//! \date Sep 12 2024
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Exports the opm_costa python module.
//!
//==============================================================================

#include <pybind11/numpy.h>
#include <pybind11/stl.h>

void export_BlackOil(pybind11::module& m);
void export_OnePhase(pybind11::module& m);

//! \brief Exports the opm_CoSTA python module.
PYBIND11_MODULE(opm_CoSTA, m)
{
    export_BlackOil(m);
    export_OnePhase(m);
}
