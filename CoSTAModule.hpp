#ifndef COSTA_MODULE_HPP
#define COSTA_MODULE_HPP

#include <opm/models/parallel/threadmanager.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicit.hpp>

#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
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
 class CoSTASimulator : public Opm::SimulatorFullyImplicit<TypeTag>
 {
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    using GlobalEqVec = Opm::GetPropType<TypeTag, Opm::Properties::GlobalEqVector>;
    using Indices = Opm::GetPropType<TypeTag, Opm::Properties::Indices>;
    using Model = typename Opm::SimulatorFullyImplicit<TypeTag>::Model;
    using PrimaryVariables = Opm::GetPropType<TypeTag, Opm::Properties::PrimaryVariables>;
    using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;

    using BVector = typename Model::BVector;
    static constexpr auto numEq = Opm::getPropValue<TypeTag, Opm::Properties::NumEq>();
    static constexpr auto numPhases = FluidSystem::numPhases;
    static constexpr auto waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static constexpr auto oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr auto gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr bool enableDissolution =
        Indices::compositionSwitchIdx != std::numeric_limits<unsigned>::max();
    static constexpr bool enableVapwat =
        Opm::getPropValue<TypeTag, Opm::Properties::EnableVapwat>();
    static constexpr bool enableDisgasInWater =
        Opm::getPropValue<TypeTag, Opm::Properties::EnableDisgasInWater>();

    //! \brief Whether a dissolution quantity is part of the state.
    static bool hasRs() { return enableDissolution && FluidSystem::enableDissolvedGas(); }
    static bool hasRv() { return enableDissolution && FluidSystem::enableVaporizedOil(); }
    static bool hasRvw() { return enableVapwat && FluidSystem::enableVaporizedWater(); }
    static bool hasRsw() { return enableDisgasInWater && FluidSystem::enableDissolvedGasInWater(); }

    //! \brief The phase whose pressure is used in the physical state.
    static unsigned referencePhase()
    {
        if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
            return oilPhaseIdx;
        }
        return FluidSystem::phaseIsActive(gasPhaseIdx) ? gasPhaseIdx : waterPhaseIdx;
    }

 public:
    explicit CoSTASimulator(Simulator& simulator)
        : Opm::SimulatorFullyImplicit<TypeTag>(simulator)
    {}

    ~CoSTASimulator()
    {
        this->finalize();
    }

    void create()
    {
        this->solver_ = this->createSolver(this->wellModel_());
    }

    //! \brief Names of the physical state variables of a cell.
    //! \details The state is independent of the primary variables used by
    //! the simulator, whose meaning switches per cell (e.g. Sg, Rs or Rv).
    //! The pressure is the pressure of the oil phase if active, otherwise
    //! of the gas phase if active, otherwise of the water phase.
    static std::vector<std::string> stateVars()
    {
        std::vector<std::string> result{"p"};
        if (FluidSystem::phaseIsActive(waterPhaseIdx) && FluidSystem::numActivePhases() > 1) {
            result.push_back("sw");
        }
        if (FluidSystem::phaseIsActive(gasPhaseIdx) && FluidSystem::phaseIsActive(oilPhaseIdx)) {
            result.push_back("sg");
        }
        if (hasRs()) {
            result.push_back("rs");
        }
        if (hasRv()) {
            result.push_back("rv");
        }
        if (hasRvw()) {
            result.push_back("rvw");
        }
        if (hasRsw()) {
            result.push_back("rsw");
        }
        return result;
    }

    //! \brief Set the solution at a time level from a physical state.
    //! \param vec Physical state, stateVars() values per cell
    //! \param idx Time level
    //! \details Primary variables and their meanings are assigned by the
    //! simulator from the resulting fluid state, as for explicit initial
    //! conditions.
    void updateSolution(const std::vector<double>& vec, std::size_t idx = 0)
    {
        using MaterialLaw = Opm::GetPropType<TypeTag, Opm::Properties::MaterialLaw>;
        using Scalar = Opm::GetPropType<TypeTag, Opm::Properties::Scalar>;

        const auto vars = stateVars();
        const auto& problem = this->simulator_.problem();
        auto& solution = this->simulator_.model().solution(idx);
        if (vec.size() != solution.size() * vars.size()) {
            throw std::runtime_error("State vector has wrong size");
        }

        const auto refPhase = referencePhase();
        auto it = vec.begin();
        for (std::size_t i = 0; i < solution.size(); ++i) {
            // The initial fluid state provides temperature, salt etc.
            auto fs = problem.initialFluidState(i);
            const Scalar p = *it++;
            Scalar sw = 0.0, sg = 0.0;
            if (FluidSystem::phaseIsActive(waterPhaseIdx) && FluidSystem::numActivePhases() > 1) {
                sw = *it++;
            }
            if (FluidSystem::phaseIsActive(gasPhaseIdx) && FluidSystem::phaseIsActive(oilPhaseIdx)) {
                sg = *it++;
            }
            else if (FluidSystem::phaseIsActive(gasPhaseIdx) &&
                     FluidSystem::phaseIsActive(waterPhaseIdx)) {
                sg = 1.0 - sw;
            }
            if (FluidSystem::numActivePhases() == 1) {
                sw = sg = 1.0;
            }
            if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
                fs.setSaturation(waterPhaseIdx, sw);
            }
            if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
                fs.setSaturation(gasPhaseIdx, sg);
            }
            if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                fs.setSaturation(oilPhaseIdx, FluidSystem::numActivePhases() == 1
                                                  ? 1.0 : 1.0 - sw - sg);
            }

            std::array<Scalar, numPhases> pc{};
            MaterialLaw::capillaryPressures(pc, problem.materialLawParams(i), fs);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (FluidSystem::phaseIsActive(phaseIdx)) {
                    fs.setPressure(phaseIdx, p + (pc[phaseIdx] - pc[refPhase]));
                }
            }

            // Storage for these only exist in the fluid state if compiled in
            if constexpr (enableDissolution) {
                if (hasRs()) {
                    fs.setRs(*it++);
                }
                if (hasRv()) {
                    fs.setRv(*it++);
                }
            }
            if constexpr (enableVapwat) {
                if (hasRvw()) {
                    fs.setRvw(*it++);
                }
            }
            if constexpr (enableDisgasInWater) {
                if (hasRsw()) {
                    fs.setRsw(*it++);
                }
            }

            solution[i].assignNaive(fs);
        }

        // The linearizer uses cached intensive quantities, which the
        // simulator only refreshes after its own Newton updates.
        if (idx == 0) {
            this->simulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
        }
    }

    //! \brief Get the physical state of the current solution.
    //! \return stateVars() values per cell
    std::vector<double> state()
    {
        auto& model = this->simulator_.model();
        model.invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);

        const auto refPhase = referencePhase();
        std::vector<double> result;
        result.reserve(model.numGridDof() * stateVars().size());
        for (unsigned i = 0; i < model.numGridDof(); ++i) {
            const auto& fs = model.intensiveQuantities(i, /*timeIdx=*/0).fluidState();
            result.push_back(Opm::getValue(fs.pressure(refPhase)));
            if (FluidSystem::phaseIsActive(waterPhaseIdx) && FluidSystem::numActivePhases() > 1) {
                result.push_back(Opm::getValue(fs.saturation(waterPhaseIdx)));
            }
            if (FluidSystem::phaseIsActive(gasPhaseIdx) && FluidSystem::phaseIsActive(oilPhaseIdx)) {
                result.push_back(Opm::getValue(fs.saturation(gasPhaseIdx)));
            }
            if (hasRs()) {
                result.push_back(Opm::getValue(fs.Rs()));
            }
            if (hasRv()) {
                result.push_back(Opm::getValue(fs.Rv()));
            }
            if (hasRvw()) {
                result.push_back(Opm::getValue(fs.Rvw()));
            }
            if (hasRsw()) {
                result.push_back(Opm::getValue(fs.Rsw()));
            }
        }
        return result;
    }

    //! \brief Component amounts in each cell for a physical state.
    //! \param vec Physical state, stateVars() values per cell
    //! \return Storage times cell volume, numEq values per cell, in the
    //! units of the residual times time (surface volumes by default)
    std::vector<double> storage(const std::vector<double>& vec)
    {
        using EqVector = Opm::GetPropType<TypeTag, Opm::Properties::EqVector>;
        using LocalResidual = Opm::GetPropType<TypeTag, Opm::Properties::LocalResidual>;
        using Scalar = Opm::GetPropType<TypeTag, Opm::Properties::Scalar>;

        auto& model = this->simulator_.model();
        const auto current = model.solution(0);
        this->updateSolution(vec, 0);
        std::vector<double> result;
        result.reserve(model.numGridDof() * numEq);
        for (unsigned i = 0; i < model.numGridDof(); ++i) {
            EqVector cell;
            LocalResidual::template computeStorage<Scalar>(cell, model.intensiveQuantities(i, 0));
            for (unsigned eq = 0; eq < numEq; ++eq) {
                result.push_back(cell[eq] * model.dofTotalVolume(i));
            }
        }
        model.solution(0) = current;
        model.invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
        return result;
    }

    //! \brief Find the physical state holding given component amounts.
    //! \param amounts Component amounts as returned by storage()
    //! \param guess Initial guess for the state
    //! \details Solves storage(state) = amounts cell by cell with Newton's
    //! method, including pressure. Primary variable switching is handled as
    //! in the simulator.
    std::vector<double> stateFromStorage(const std::vector<double>& amounts,
                                         const std::vector<double>& guess)
    {
        using Evaluation = Opm::GetPropType<TypeTag, Opm::Properties::Evaluation>;
        using GlobalEqVector = Opm::GetPropType<TypeTag, Opm::Properties::GlobalEqVector>;
        using LocalResidual = Opm::GetPropType<TypeTag, Opm::Properties::LocalResidual>;
        using Scalar = Opm::GetPropType<TypeTag, Opm::Properties::Scalar>;
        constexpr int maxIter = 50;
        constexpr Scalar tolerance = 1e-12;

        auto& model = this->simulator_.model();
        const unsigned numCells = model.numGridDof();
        if (amounts.size() != numCells * numEq) {
            throw std::runtime_error("Component amounts have wrong size");
        }
        std::array<Scalar, numEq> scale{};
        for (unsigned i = 0; i < numCells; ++i) {
            for (unsigned eq = 0; eq < numEq; ++eq) {
                scale[eq] = std::max(scale[eq], std::abs(amounts[i * numEq + eq]));
            }
        }

        const auto current = model.solution(0);
        this->updateSolution(guess, 0);
        GlobalEqVector residual(numCells), update(numCells);
        bool converged = false;
        for (int iter = 0; iter < maxIter && !converged; ++iter) {
            converged = true;
            for (unsigned i = 0; i < numCells; ++i) {
                Dune::FieldVector<Evaluation, numEq> cell;
                LocalResidual::template computeStorage<Evaluation>(cell,
                                                                   model.intensiveQuantities(i, 0));
                Dune::FieldMatrix<Scalar, numEq, numEq> jacobian;
                const Scalar volume = model.dofTotalVolume(i);
                for (unsigned eq = 0; eq < numEq; ++eq) {
                    residual[i][eq] = cell[eq].value() * volume - amounts[i * numEq + eq];
                    for (unsigned pv = 0; pv < numEq; ++pv) {
                        jacobian[eq][pv] = cell[eq].derivative(pv) * volume;
                    }
                    if (std::abs(residual[i][eq]) > tolerance * scale[eq]) {
                        converged = false;
                    }
                }
                jacobian.solve(update[i], residual[i]);
            }
            if (!converged) {
                model.newtonMethod().applyUpdate(model.solution(0), model.solution(0),
                                                 update, residual);
                model.invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
            }
        }
        if (!converged) {
            model.solution(0) = current;
            model.invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
            throw std::runtime_error("Finding the state for given component amounts did not converge");
        }

        auto result = this->state();
        model.solution(0) = current;
        model.invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
        return result;
    }

    //! \brief Set the state at the start of the time step.
    //! \details Besides the primary variables, this updates the cached
    //! storage term of the old time level. The linearizer otherwise
    //! recycles the storage of the iterate in the first Newton iteration,
    //! which assumes that the iterate equals the start-of-step state.
    void setOldState(const std::vector<double>& vec)
    {
        using EqVector = Opm::GetPropType<TypeTag, Opm::Properties::EqVector>;
        using LocalResidual = Opm::GetPropType<TypeTag, Opm::Properties::LocalResidual>;
        using Scalar = Opm::GetPropType<TypeTag, Opm::Properties::Scalar>;

        auto& model = this->simulator_.model();
        this->updateSolution(vec, 1);
        if (!model.enableStorageCache()) {
            model.invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/1);
            return;
        }

        // Evaluate the storage of the old state using the time level 0 cache
        const auto current = model.solution(0);
        this->updateSolution(vec, 0);
        for (unsigned i = 0; i < model.numGridDof(); ++i) {
            EqVector storage;
            LocalResidual::template computeStorage<Scalar>(storage,
                                                           model.intensiveQuantities(i, 0));
            model.updateCachedStorage(i, /*timeIdx=*/1, storage);
        }
        model.solution(0) = current;
        model.invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
    }

    //! \brief Restart the Newton iteration count at iteration 1.
    //! \details Iteration 0 would recycle the storage of the iterate as
    //! old storage, see setOldState(). Being below the minimum number of
    //! iterations forces a Newton update even if the residual is within
    //! tolerances, which CoSTA corrections often are.
    void restartIterations()
    {
        auto& problem = this->simulator_.problem();
        problem.resetIterationForNewTimestep();
        problem.markTimestepInitialized();
        problem.advanceIteration();
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
        this->solver_->step(timer, nullptr);

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
        this->restartIterations();
        Opm::SimulatorReportSingle dummy;
        this->solver_->model().initialLinearization(dummy, 1, 2, timer);
    }

    //! \brief Solve the corrected equations of the current step.
    //! \return False if the Newton iterations failed
    bool solveCorrected(Opm::SimulatorTimer& timer)
    {
        this->restartIterations();
        const auto& problem = this->simulator_.problem();
        const int maxIter = this->solver_->model().param().newton_max_iter_;
        try {
            while (!this->solver_->model().nonlinearIteration(timer, *this->solver_).converged) {
                if (problem.iterationContext().iteration() > maxIter) {
                    return false;
                }
            }
        }
        catch (const std::exception&) {
            return false;
        }
        catch (const Dune::Exception&) {
            return false;
        }
        return true;
    }

    //! \brief Reset to the start of the step after failed Newton iterations.
    //! \details As for a time step chop, this restores the solution and wells.
    void rollback()
    {
        this->simulator_.problem().updateFailed();
    }

    //! \brief Finish the current step.
    void finishStep(Opm::SimulatorTimer& timer)
    {
        this->simulator_.problem().endTimeStep();
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
    //! \param args Additional command line arguments for the simulator.
    //!             These override the defaults set by the module.
    CoSTAModule(const std::string& infile, bool verbose,
                const std::vector<std::string>& args = {})
    {
        const std::string verb = verbose ? "--enable-terminal-output=true"
                                         : "--enable-terminal-output=false";
        const std::string output_mode = verbose ? "--output-mode=all"
                                                : "--output-mode=none";
        const std::vector<std::string> defaults = {verb, output_mode,
                                                   "--matrix-add-well-contributions=true",
                                                   "--linear-solver=ilu0",
                                                   "--enable-adaptive-time-stepping=0"};
        // Options given in args override the defaults
        const auto key = [](const std::string& option)
                         { return option.substr(0, option.find('=')); };
        std::vector<const char*> arg = {"opm_CoSTA", infile.c_str()};
        for (const auto& option : defaults) {
            if (std::none_of(args.begin(), args.end(),
                             [&](const std::string& a) { return key(a) == key(option); }))
            {
                arg.push_back(option.c_str());
            }
        }
        for (const auto& a : args) {
            arg.push_back(a.c_str());
        }
        int argc = arg.size();
        arg.push_back(nullptr);
        char** argv = const_cast<char**>(arg.data());

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
#ifdef RESERVOIR_COUPLING_ENABLED
        sim->init(timer, 0, nullptr);
#else
        sim->init(timer);
#endif
        sim->create();

        ndof = sim->model().simulator().model().numGridDof() * numEq;
        initial_ = sim->state();
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
            // The step starts from uprev, and uses it as initial guess
            sim->updateSolution(uprev, 1);
            sim->updateSolution(uprev, 0);
        }
        sim->predictStep(timer);

        return sim->state();
    }

    //! \brief Perform a residual calculation step in a CoSTA loop.
    //! \param mu Model parameters
    //! \param uprev State to make a time-step from
    //! \param unext State to calculate residual for
    std::vector<double>
    residual(const ParameterMap&,
             const std::vector<double>& uprev, const std::vector<double>& unext)
    {
        sim->setOldState(uprev);
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
        return this->correctImpl(uprev, sigma);
    }

    //! \brief Returns the active interior faces used by the face-based CoSTA API.
    //! \details Each face appears once, oriented from its \a minus cell to
    //! its \a plus cell (minus < plus in active cell numbering). A positive
    //! flux means flow from minus to plus. The ordering of the returned
    //! arrays defines the ordering of all face vectors (fluxes, tau).
    //! Boundary faces are not included (no-flow boundaries). Non-neighbor
    //! connections are included. The direction is the FaceDir::DirEnum value
    //! of the face seen from the minus cell (0 if unknown).
    pybind11::dict faces()
    {
        this->ensureFaces();
        std::vector<std::array<int,2>> cells, cartesian;
        std::vector<double> trans, area;
        std::vector<int> dir;
        cells.reserve(faces_.size());
        cartesian.reserve(faces_.size());
        trans.reserve(faces_.size());
        area.reserve(faces_.size());
        dir.reserve(faces_.size());
        const auto& nbInfo = sim->model().simulator().model().linearizer().getNeighborInfo();
        for (const auto& face : faces_) {
            const auto& info = nbInfo[face.minus][face.loc].res_nbinfo;
            cells.push_back({face.minus, face.plus});
            cartesian.push_back({static_cast<int>(sim->vanguard().cartesianIndex(face.minus)),
                                 static_cast<int>(sim->vanguard().cartesianIndex(face.plus))});
            trans.push_back(info.trans);
            area.push_back(info.faceArea);
            dir.push_back(static_cast<int>(info.faceDir));
        }

        pybind11::dict result;
        result["cells"] = cells;
        result["cartesian"] = cartesian;
        result["trans"] = trans;
        result["area"] = area;
        result["direction"] = dir;
        return result;
    }

    //! \brief Evaluate integrated face fluxes for a state.
    //! \param mu Model parameters
    //! \param u State to evaluate fluxes for
    //! \return Fluxes ordered as faces(), numEq values per face, in the units
    //! of the residual (surface volume rates by default) and oriented minus to plus
    std::vector<double>
    fluxes(const ParameterMap&, const std::vector<double>& u)
    {
        sim->updateSolution(u, 0);

        return this->computeFluxes();
    }

    //! \brief Component amounts in each cell for a state.
    //! \param mu Model parameters
    //! \param u State
    //! \return numEq values per cell, summable over cells
    std::vector<double> storage(const ParameterMap&, const std::vector<double>& u)
    {
        return sim->storage(u);
    }

    //! \brief Find the state holding given component amounts.
    //! \param mu Model parameters
    //! \param amounts numEq values per cell, as from storage()
    //! \param guess Initial guess for the state
    //! \details Used for a mass consistent transfer of states between
    //! grids, by summing amounts over the cells being merged.
    std::vector<double> stateFromStorage(const ParameterMap&,
                                         const std::vector<double>& amounts,
                                         const std::vector<double>& guess)
    {
        return sim->stateFromStorage(amounts, guess);
    }

    //! \brief Apply the face-cell incidence operator D to a face vector.
    //! \param q Face vector, ordered as faces(), numEq values per face
    //! \return Per-cell balance, numEq values per cell (+q on minus, -q on plus)
    std::vector<double> divergence(const std::vector<double>& q)
    {
        this->ensureFaces();
        if (q.size() != faces_.size() * numEq) {
            throw std::runtime_error("Face vector has wrong size");
        }
        std::vector<double> result(ndof, 0.0);
        for (std::size_t f = 0; f < faces_.size(); ++f) {
            for (unsigned eq = 0; eq < numEq; ++eq) {
                result[faces_[f].minus * numEq + eq] += q[f * numEq + eq];
                result[faces_[f].plus * numEq + eq] -= q[f * numEq + eq];
            }
        }
        return result;
    }

    //! \brief Perform a face correction step in a Face CoSTA loop.
    //! \param mu Model parameters
    //! \param uprev State to make a time-step from
    //! \param tau Face flux correction, ordered as faces()
    //! \param sigma Optional residual correction, numEq values per cell,
    //!              as for correct(). Empty for none.
    //! \return Tuple of corrected state and corrected fluxes, q = Q(u) + tau
    //! \details The correction tau is added to every face flux, which is
    //! equivalent to the cell source -D*tau. The returned fluxes therefore
    //! satisfy the cell balance, with sigma added to the sources, for any
    //! tau. Sigma is meant for what face fluxes cannot represent, such as
    //! differing well models.
    std::tuple<std::vector<double>, std::vector<double>>
    correctFace(const ParameterMap&,
                const std::vector<double>& uprev, const std::vector<double>& tau,
                const std::vector<double>& sigma)
    {
        std::vector<double> total = this->divergence(tau);
        if (!sigma.empty() && sigma.size() != total.size()) {
            throw std::runtime_error("Residual correction has wrong size");
        }
        for (std::size_t i = 0; i < total.size(); ++i) {
            total[i] = (sigma.empty() ? 0.0 : sigma[i]) - total[i];
        }
        auto u = this->correctImpl(uprev, total);

        // Fluxes are evaluated at the corrected state.
        sim->model().simulator().model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
        auto q = this->computeFluxes();
        for (std::size_t i = 0; i < q.size(); ++i) {
            q[i] += tau[i];
        }
        return {u, q};
    }

    //! \brief Get the IDs of all Dirichlet DoFs.
    std::vector<int> dirichletDofs()
    {
        return {};
    }

    //! \brief Returns initial condition for solution.
    std::vector<double> initialCondition(const ParameterMap&)
    {
        return initial_;
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

    //! Number of equations in simulator, i.e. the size of residuals and
    //! residual corrections. States have nstate = ncells * state_vars entries.
    size_t ndof;

    //! \brief Static helper to export to python.
    //! \param m Module to export to
    //! \param name Name of python class
    static void pyExport(pybind11::module& m, const char* name)
    {
      pybind11::class_<CoSTAModule<TypeTag>>(m, name)
          .def(pybind11::init<const std::string&, bool, const std::vector<std::string>&>(),
               pybind11::arg("infile"),
               pybind11::arg("verbose") = true,
               pybind11::arg("args") = std::vector<std::string>{})
          .def("correct", &CoSTAModule<TypeTag>::correct)
          .def("predict", &CoSTAModule<TypeTag>::predict)
          .def("residual", &CoSTAModule<TypeTag>::residual)
          .def("dirichlet_dofs", &CoSTAModule<TypeTag>::dirichletDofs)
          .def("initial_condition", &CoSTAModule<TypeTag>::initialCondition)
          .def("anasol", &CoSTAModule<TypeTag>::anaSol)
          .def("qi", &CoSTAModule<TypeTag>::getQI)
          .def("faces", &CoSTAModule<TypeTag>::faces)
          .def("fluxes", &CoSTAModule<TypeTag>::fluxes)
          .def("divergence", &CoSTAModule<TypeTag>::divergence)
          .def("storage", &CoSTAModule<TypeTag>::storage)
          .def("state_from_storage", &CoSTAModule<TypeTag>::stateFromStorage)
          .def("correct_face", &CoSTAModule<TypeTag>::correctFace,
               pybind11::arg("mu"), pybind11::arg("uprev"), pybind11::arg("tau"),
               pybind11::arg("sigma") = std::vector<double>{})
          .def_property_readonly("nfaces", [](CoSTAModule<TypeTag>& self)
                                 { self.ensureFaces(); return self.faces_.size(); })
          .def_readonly("ndof", &CoSTAModule<TypeTag>::ndof)
          .def_property_readonly("ncells", [](CoSTAModule<TypeTag>& self)
                                 { return self.ndof / numEq; })
          .def_property_readonly_static("state_vars", [](pybind11::object)
                                        { return CoSTASimulator<TypeTag>::stateVars(); })
          .def_property_readonly("nstate", [](CoSTAModule<TypeTag>& self)
                                 { return self.ndof / numEq *
                                          CoSTASimulator<TypeTag>::stateVars().size(); });
    }

protected:
    //! \brief An interior face, oriented from minus to plus cell.
    struct Face
    {
        int minus; //!< Active index of minus cell
        int plus; //!< Active index of plus cell
        unsigned loc; //!< Position of plus cell in neighbor info of minus cell
    };

    //! \brief Set up the face list from the TPFA connectivity.
    void ensureFaces()
    {
        if (!faces_.empty()) {
            return;
        }
        auto& linearizer = msim->model().linearizer();
        if (linearizer.getNeighborInfo().empty()) {
            // Connectivity is set up on first linearization
            linearizer.linearizeDomain();
        }
        const auto& nbInfo = linearizer.getNeighborInfo();
        for (int cell = 0; cell < static_cast<int>(nbInfo.size()); ++cell) {
            unsigned loc = 0;
            for (const auto& nb : nbInfo[cell]) {
                if (static_cast<int>(nb.neighbor) > cell) {
                    faces_.push_back({cell, static_cast<int>(nb.neighbor), loc});
                }
                ++loc;
            }
        }
    }

    //! \brief Compute face fluxes from the current intensive quantities.
    std::vector<double> computeFluxes()
    {
        using LocalResidual = Opm::GetPropType<TypeTag, Opm::Properties::LocalResidual>;
        using RateVector = Opm::GetPropType<TypeTag, Opm::Properties::RateVector>;

        this->ensureFaces();
        const auto& model = sim->model().simulator().model();
        const auto& problem = sim->model().simulator().problem();
        const auto& nbInfo = model.linearizer().getNeighborInfo();
        std::vector<double> result(faces_.size() * numEq);
        RateVector flux, darcy;
        for (std::size_t f = 0; f < faces_.size(); ++f) {
            const auto& face = faces_[f];
            const auto& info = nbInfo[face.minus][face.loc].res_nbinfo;
            LocalResidual::computeFlux(flux, darcy, face.minus, face.plus,
                                       model.intensiveQuantities(face.minus, /*timeIdx=*/0),
                                       model.intensiveQuantities(face.plus, /*timeIdx=*/0),
                                       info, problem.moduleParams());
            for (unsigned eq = 0; eq < numEq; ++eq) {
                result[f * numEq + eq] = Opm::getValue(flux[eq]) * info.faceArea;
            }
        }
        return result;
    }

    //! \brief Convert a residual correction to deck source terms.
    //! \param sigma Correction in residual units, numEq values per cell
    //! \details The correction is added on top of any SOURCE in the deck.
    //! Residuals are in surface volume rates if the model conserves surface
    //! volumes, while source terms are mass rates.
    Opm::Source residualToSource(const std::vector<double>& sigma) const
    {
        using Indices = Opm::GetPropType<TypeTag, Opm::Properties::Indices>;
        constexpr bool surfaceVolume =
            Opm::getPropValue<TypeTag, Opm::Properties::BlackoilConserveSurfaceVolume>();
        static constexpr std::array phases = {
            std::pair{FluidSystem::waterPhaseIdx, Opm::SourceComponent::WATER},
            std::pair{FluidSystem::oilPhaseIdx, Opm::SourceComponent::OIL},
            std::pair{FluidSystem::gasPhaseIdx, Opm::SourceComponent::GAS},
        };

        if (sigma.size() != ndof) {
            throw std::runtime_error("Residual correction has wrong size");
        }

        Opm::Source source = sim->vanguard().schedule()[timer.currentStepNum()].source();
        const auto& problem = sim->model().simulator().problem();
        for (std::size_t i = 0; i < ndof / numEq; ++i) {
            std::array<int,3> ijk;
            sim->vanguard().cartesianCoordinate(i, ijk);
            for (const auto& [phaseIdx, component] : phases) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }
                const unsigned compIdx = FluidSystem::solventComponentIndex(phaseIdx);
                const unsigned eqIdx = Indices::conti0EqIdx +
                                       FluidSystem::canonicalToActiveCompIdx(compIdx);
                double rate = sigma[i * numEq + eqIdx];
                if constexpr (surfaceVolume) {
                    rate *= FluidSystem::referenceDensity(phaseIdx, problem.pvtRegionIndex(i));
                }
                Opm::Source::SourceCell cell;
                cell.component = component;
                cell.rate = source.rate(ijk, component) + rate;
                source.addSourceCell(ijk, cell);
            }
        }
        return source;
    }

    //! \brief Perform a correction step with a residual correction.
    //! \param uprev State to make a time-step from
    //! \param sigma Residual correction, numEq values per cell
    std::vector<double>
    correctImpl(const std::vector<double>& uprev, const std::vector<double>& sigma)
    {
        sim->updateSolution(uprev, 0);
        Opm::Source source = this->residualToSource(sigma);
        auto& origSource = const_cast<Opm::Source&>(sim->vanguard().schedule()[timer.currentStepNum()].source());
        std::swap(source, origSource);
        const bool converged = sim->solveCorrected(timer);
        std::swap(source, origSource);
        if (!converged) {
            // The step is not finished, a new correction can be attempted
            sim->rollback();
            throw std::runtime_error("CoSTA correction did not converge");
        }
        sim->finishStep(timer);

        return sim->state();
    }

    //! \brief Convert a residual to a vector with numEq values per cell.
    std::vector<double> eqVecToVec(const GlobalEqVec& vec)
    {
        std::vector<double> result;
        result.reserve(vec.size() * numEq);
        for (std::size_t i = 0; i < vec.size(); ++i) {
            for (unsigned eq = 0; eq < numEq; ++eq) {
                result.push_back(vec[i][eq]);
            }
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
    std::vector<Face> faces_; //!< Interior faces in face vector ordering
    std::vector<double> initial_; //!< Initial state
};

#endif
