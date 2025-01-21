#ifndef COSTA_MODULE_HPP
#define COSTA_MODULE_HPP

#include <opm/models/parallel/threadmanager.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoil.hpp>

#include <map>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

class CoSTAMain : public Opm::Main
{
public:
    using Opm::Main::Main;

    template<class TypeTag>
    bool initialize()
    {
        int exitCode = EXIT_SUCCESS;
        if (initialize_<TypeTag>(exitCode)) {
            // TODO: check that this deck really represents a blackoil
            // case. E.g. check that number of phases == 3
            this->setupVanguard();
        } else {
            throw std::runtime_error("Error initializing simulator");
        }

        return exitCode == EXIT_SUCCESS;
    }
};

 template<class TypeTag>
 class CoSTASimulator : public Opm::SimulatorFullyImplicitBlackoil<TypeTag>
 {
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    using GlobalEqVec = Opm::GetPropType<TypeTag, Opm::Properties::GlobalEqVector>;
    using Indices = Opm::GetPropType<TypeTag, Opm::Properties::Indices>;
    using Model = typename Opm::SimulatorFullyImplicitBlackoil<TypeTag>::Model;
    using PrimaryVariables = Opm::GetPropType<TypeTag, Opm::Properties::PrimaryVariables>;
    using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;

    using BVector = typename Model::BVector;
    static constexpr auto numEq = Opm::getPropValue<TypeTag, Opm::Properties::NumEq>();

 public:
    explicit CoSTASimulator(Simulator& simulator)
        : Opm::SimulatorFullyImplicitBlackoil<TypeTag>(simulator)
    {}

    ~CoSTASimulator()
    {
        this->finalize();
    }

    void create()
    {
        this->solver_ = this->createSolver(this->wellModel_());
    }

    void updateSolution(const std::vector<double>& vec, std::size_t idx = 0)
    {
        auto& result = this->solver_->model().simulator().model().solution(idx);
        auto it = vec.begin();
        for (std::size_t i = 0; i < result.size(); ++i, it += numEq) {
            auto& pvVars = result[i];
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) &&
                pvVars.primaryVarsMeaningGas() != PrimaryVariables::GasMeaning::Disabled)
            {
                pvVars[Indices::compositionSwitchIdx] = *(it + Indices::compositionSwitchIdx);
            }
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) &&
                pvVars.primaryVarsMeaningWater() != PrimaryVariables::WaterMeaning::Disabled)
            {
                pvVars[Indices::waterSwitchIdx] = *(it + Indices::waterSwitchIdx);
            }
            pvVars[Indices::pressureSwitchIdx] = *(it + Indices::pressureSwitchIdx);
        }

        // unsigned nc = this->simulator_.model().numGridDof();
        // BVector x(nc);
        // this->solver_->model().updateSolution(x);
    }

    bool predictStep(Opm::SimulatorTimer& timer)
    {
        if (this->schedule().exitStatus().has_value()) {
            if (this->terminalOutput_) {
                Opm::OpmLog::info("Stopping simulation since EXIT was triggered by an action keyword.");
            }
            this->report_.success.exit_status = this->schedule().exitStatus().value();
            return false;
        }

        // write the inital state at the report stage
        if (timer.initialStep()) {
            Dune::Timer perfTimer;
            perfTimer.start();

            this->simulator_.setEpisodeIndex(-1);
            this->simulator_.setEpisodeLength(0.0);
            this->simulator_.setTimeStepSize(0.0);
            this->wellModel_().beginReportStep(timer.currentStepNum());
            this->simulator_.problem().writeOutput(false);

            this->report_.success.output_write_time += perfTimer.stop();
        }

        this->simulator_.startNextEpisode(
            this->simulator_.startTime()
               + this->schedule().seconds(timer.currentStepNum()),
            timer.currentStepLength());
        this->simulator_.setEpisodeIndex(timer.currentStepNum());

        this->solver_->model().beginReportStep();

        // solve for complete report step
        this->solver_->step(timer);

        return true;
    }

    auto& vanguard()
    {
        return this->simulator_.vanguard();
    }

    const auto& solution(std::size_t idx)
    {
        return this->solver_->model().simulator().model().solution(idx);
    }

    void linearize(const Opm::SimulatorTimer& timer)
    {
        Opm::SimulatorReportSingle dummy;
        this->solver_->model().initialLinearization(dummy, 0, 1, 2, timer);
    }

    void correctStep(Opm::SimulatorTimer& timer)
    {
        this->solver_->model().nonlinearIteration(2, timer, *this->solver_);
        this->solver_->model().afterStep(timer);
        this->simulator_.problem().writeOutput(false);
        this->solver_->model().endReportStep();
        ++timer;
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
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    using GlobalEqVec = Opm::GetPropType<TypeTag, Opm::Properties::GlobalEqVector>;
    using ModelSimulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;
    using PrimaryVariables = Opm::GetPropType<TypeTag, Opm::Properties::PrimaryVariables>;
    using SolutionVec = Opm::GetPropType<TypeTag, Opm::Properties::SolutionVector>;

    static constexpr auto numEq = Opm::getPropValue<TypeTag, Opm::Properties::NumEq>();

    //! \brief Constructor.
    //! \param infile Input file to parse
    //! \param verbose If false, suppress stdout output
    CoSTAModule(const std::string& infile, bool verbose)
    {
        const std::string verb = verbose ? "--enable-terminal-output=true"
                                         : "--enable-terminal-output=false";
        const std::string output_mode = verbose ? "--output-mode=all"
                                                : "--output-mode=none";
        const char* arg[] = {"opm_CoSTA", infile.c_str(),
                             verb.c_str(), output_mode.c_str(),
                            "--matrix-add-well-contributions=true",
                            "--linear-solver=ilu0",
                            "--enable-adaptive-time-stepping=0",
                            nullptr};
        char** argv = const_cast<char**>(arg);
        int argc = 7;

        // read deck and initialize eclipse structures
        main = std::make_unique<CoSTAMain>(argc, argv);
        if (!main->initialize<TypeTag>()) {
            throw std::runtime_error("Error initializing simulator");
        }

        Opm::FlowMain<TypeTag>::setupParameters_(argc, argv, Opm::FlowGenericVanguard::comm());

        // setup the model simulator
        msim = std::make_unique<ModelSimulator>(Opm::FlowGenericVanguard::comm(), false);
        msim->model().applyInitialSolution();
        sim = std::make_unique<CoSTASimulator<TypeTag>>(*msim);

        timer.init(msim->vanguard().schedule(), 0);
        msim->setEpisodeIndex(timer.currentStepNum());
        msim->model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
        sim->init(timer, 0, nullptr);
        sim->create();

        ndof = sim->model().simulator().model().numGridDof() * numEq;
    }

    ~CoSTAModule()
    {
        sim.reset();
        msim.reset();
        main.reset();
    }

    //! \brief Perform a prediction step in a CoSTA loop.
    //! \param mu Model parameters
    //! \param uprev State to make a time-step from
    std::vector<double> predict(const std::vector<double>& uprev)
    {
        if (!timer.initialStep()) {
            sim->updateSolution(uprev, 1);
        }
        sim->predictStep(timer);

        return priVarsToVec(sim->model().simulator().model().solution(0));
    }

    //! \brief Perform a residual calculation step in a CoSTA loop.
    //! \param mu Model parameters
    //! \param uprev State to make a time-step from
    //! \param unext State to calculate residual for
    std::vector<double>
    residual(const ParameterMap&,
             const std::vector<double>& uprev, const std::vector<double>& unext)
    {
        sim->updateSolution(uprev, 1);
        sim->updateSolution(unext, 0);
        sim->linearize(timer);

        return eqVecToVec(sim->model().simulator().model().linearizer().residual());
    }

    //! \brief Perform a correction step in a CoSTA loop.
    //! \param mu Model parameters
    //! \param uprev State to make a time-step from
    //! \param sigma Right-hand-side correction to use
    std::vector<double>
    correct(const ParameterMap&,
            const std::vector<double>& uprev, const std::vector<double>& sigma)
    {
        using Indices = Opm::GetPropType<TypeTag, Opm::Properties::Indices>;
        sim->updateSolution(uprev, 0);
        Opm::Source source;
        for (std::size_t i = 0; i < ndof; ++i) {
            const auto& pvVars = sim->solution(/*timeIdx=*/0)[i];
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) &&
                pvVars.primaryVarsMeaningGas() != PrimaryVariables::GasMeaning::Disabled)
            {
                Opm::Source::SourceCell cell;
                sim->vanguard().cartesianCoordinate(i, cell.ijk);
                cell.component = Opm::SourceComponent::GAS;
                cell.rate = sigma[i * numEq + Indices::compositionSwitchIdx];
                source.addSourceCell(cell);
            }
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) &&
                pvVars.primaryVarsMeaningWater() != PrimaryVariables::WaterMeaning::Disabled)
            {
                Opm::Source::SourceCell cell;
                sim->vanguard().cartesianCoordinate(i, cell.ijk);
                cell.component = Opm::SourceComponent::WATER;
                cell.rate = sigma[i * numEq + Indices::waterSwitchIdx];
                source.addSourceCell(cell);
            }
        }
        auto& origSource = const_cast<Opm::Source&>(sim->vanguard().schedule()[timer.currentStepNum()].source());
        std::swap(source, origSource);
        sim->correctStep(timer);
        std::swap(source, origSource);

        return priVarsToVec(sim->solution(0));
    }

    //! \brief Get the IDs of all Dirichlet DoFs.
    std::vector<int> dirichletDofs()
    {
        return {};
    }

    //! \brief Returns initial condition for solution.
    std::vector<double> initialCondition(const ParameterMap&)
    {
        return {};
    }

    //! \brief Returns the analytical solution.
    //! \param mu Map of parameters
    std::map<std::string,std::vector<double>> anaSol(const ParameterMap&)
    {
        return {};
    }

    //! \brief Returns a quantity of interest from a solution vector.
    //! \param mu Map of parameters
    //! \param u Solution vector to extract QI from
    //! \param qi Name of quantity of interest
    std::vector<double>
    getQI(const ParameterMap& ,
          const std::vector<double>&,
          const std::string&)
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
    std::vector<double> priVarsToVec(const SolutionVec& vec)
    {
        using Indices = Opm::GetPropType<TypeTag, Opm::Properties::Indices>;
        std::vector<double> result(vec.size() * numEq);
        auto it = result.begin();
        for (std::size_t i = 0; i < vec.size(); ++i, it += numEq) {
            const auto& pVars = vec[i];
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) &&
                pVars.primaryVarsMeaningGas() != PrimaryVariables::GasMeaning::Disabled)
            {
                *(it + Indices::compositionSwitchIdx) = pVars[Indices::compositionSwitchIdx];
            }
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) &&
                pVars.primaryVarsMeaningWater() != PrimaryVariables::WaterMeaning::Disabled)
            {
                *(it + Indices::waterSwitchIdx) = pVars[Indices::waterSwitchIdx];
            }
            *(it + Indices::pressureSwitchIdx) = pVars[Indices::pressureSwitchIdx];
        }
        return result;
    }

    std::vector<double> eqVecToVec(const GlobalEqVec& vec)
    {
        using Indices = Opm::GetPropType<TypeTag, Opm::Properties::Indices>;
        std::vector<double> result(vec.size() * numEq);
        const auto& pVars = sim->solution(0)[0];
        auto it = result.begin();
        for (std::size_t i = 0; i < vec.size(); ++i, it += numEq) {
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) &&
                pVars.primaryVarsMeaningGas() != PrimaryVariables::GasMeaning::Disabled)
            {
                *(it + Indices::compositionSwitchIdx) = vec[i][Indices::compositionSwitchIdx];
            }
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) &&
                pVars.primaryVarsMeaningWater() != PrimaryVariables::WaterMeaning::Disabled)
            {
                *(it + Indices::waterSwitchIdx) = vec[i][Indices::waterSwitchIdx];
            }
            *(it + Indices::pressureSwitchIdx) = vec[i][Indices::pressureSwitchIdx];
        }
        return result;
    }

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
    void setParam(const std::string&, double)
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
    std::unique_ptr<ModelSimulator> msim;
    std::unique_ptr<CoSTASimulator<TypeTag>> sim;
    Opm::SimulatorTimer timer;
};

#endif
