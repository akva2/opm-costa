#ifndef COSTA_MODULE_HPP
#define COSTA_MODULE_HPP

#include <opm/simulators/flow/Main.hpp>

#include <map>
#include <pybind11/pybind11.h>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

class CoSTAMain : public Opm::Main
{
public:
    using Opm::Main::Main;

    template<class TypeTag>
    int initialize()
    {
        int exitCode = EXIT_SUCCESS;
        if (initialize_<TypeTag>(exitCode)) {
            // TODO: check that this deck really represents a blackoil
            // case. E.g. check that number of phases == 3
            this->setupVanguard();
        } else {
            throw std::runtime_error("Error initializing simulator");
        }
    }
};

/*!
 \brief Class adding a CoSTA interface to a simulator.
*/

template<class TypeTag>
class CoSTAModule
{
public:
    using ParameterMap = std::map<std::string, std::variant<double,std::vector<double>>>; //!< Map for parameters

    //! \brief Constructor.
    //! \param infile Input file to parse
    //! \param verbose If false, suppress stdout output
    CoSTAModule(const std::string& infile, bool verbose)
    {
        const char* argv[] = {"opm_CoSTA", infile.c_str(), nullptr};
        main = std::make_unique<CoSTAMain>(2, const_cast<char**>(argv));
        main->initialize<TypeTag>();
    }

    //! \brief Perform a prediction step in a CoSTA loop.
    //! \param mu Model parameters
    //! \param uprev State to make a time-step from
    std::vector<double> predict(const std::vector<double>& uprev)
    {
        return {};
    }

    //! \brief Perform a residual calculation step in a CoSTA loop.
    //! \param mu Model parameters
    //! \param uprev State to make a time-step from
    //! \param unext State to calculate residual for
    std::vector<double>
    residual(const ParameterMap& mu,
             const std::vector<double>& uprev, const std::vector<double>& unext)
    {
        return {};
    }

    //! \brief Perform a correction step in a CoSTA loop.
    //! \param mu Model parameters
    //! \param uprev State to make a time-step from
    //! \param sigma Right-hand-side correction to use
    std::vector<double>
    correct(const ParameterMap& mu,
            const std::vector<double>& uprev, const std::vector<double>& sigma)
    {
        return {};
    }

    //! \brief Get the IDs of all Dirichlet DoFs.
    std::vector<int> dirichletDofs()
    {
        return {};
    }

    //! \brief Returns initial condition for solution.
    std::vector<double> initialCondition(const ParameterMap& mu)
    {
        return {};
    }

    //! \brief Returns the analytical solution.
    //! \param mu Map of parameters
    std::map<std::string,std::vector<double>> anaSol(const ParameterMap& mu)
    {
        return {};
    }

    //! \brief Returns a quantity of interest from a solution vector.
    //! \param mu Map of parameters
    //! \param u Solution vector to extract QI from
    //! \param qi Name of quantity of interest
    std::vector<double>
    getQI(const ParameterMap& mu,
          const std::vector<double>& u,
          const std::string& qi)
    {
        return {};
    }

    size_t ndof; //!< Number of degrees of freedom in simulator

    //! \brief Static helper to export to python.
    //! \param m Module to export to
    //! \param name Name of python class
    static void pyExport(pybind11::module& m, const char* name)
    {
      pybind11::class_<CoSTAModule<TypeTag>>(m, name)
          .def(pybind11::init<const std::string&, bool>(),
               pybind11::arg("infile"),
               pybind11::arg("verbose") = true)
          .def("correct", &CoSTAModule<TypeTag>::correct)
          .def("predict", &CoSTAModule<TypeTag>::predict)
          .def("residual", &CoSTAModule<TypeTag>::residual)
          .def("dirichlet_dofs", &CoSTAModule<TypeTag>::dirichletDofs)
          .def("initial_condition", &CoSTAModule<TypeTag>::initialCondition)
          .def("anasol", &CoSTAModule<TypeTag>::anaSol)
          .def("qi", &CoSTAModule<TypeTag>::getQI)
          .def_readonly("ndof", &CoSTAModule<TypeTag>::ndof);
    }

protected:
    //! \brief Get a scalar parameter from map.
    //! \param map Map with parameters
    //! \param key Name of parameter
    double getScalarParameter(const ParameterMap& map, const std::string& key) const
    {
        const auto it = map.find(key);
        if (it == map.end())
            throw std::runtime_error("Need "+key+" in parameters");
        if (!std::holds_alternative<double>(it->second))
            throw std::runtime_error(key+" needs to be a double");
        return std::get<double>(it->second);
    }

    //! \brief Helper function to set a parameter in simulator.
    //! \param name Name of parameter
    //! \param value Value of parameter
    void setParam(const std::string& name, double value)
    {
    }

    //! \brief Set parameters from map in simulator.
    //! \param map Map of parameters
    void setParameters(const ParameterMap& map)
    {
        static const std::vector<std::string> blackList = {"dt", "t"};
        for (const auto& entry : map)
            if (std::holds_alternative<double>(entry.second))
                if (std::find(blackList.begin(), blackList.end(), entry.first) == blackList.end())
                    this->setParam(entry.first, std::get<double>(entry.second));
    }

    //! \brief Extract time stepping parameters from parameter map.
    std::tuple<double, double> getTimeParams(const ParameterMap& map) const
    {
        double dt = this->getScalarParameter(map, "dt");
        double t = dt;
        if (map.find("t") != map.end())
            t = this->getScalarParameter(map, "t");

        return std::make_tuple(dt,t);
    }

    std::unique_ptr<CoSTAMain> main;
};

#endif
