#include <MatrixExponential.hpp>
#include <ProductOfExponentials.hpp>
#include <ScrewTheoryIkProblem.hpp>
#include <ScrewTheoryIkSubproblems.hpp>
#include <ConfigurationSelector.hpp>

#include "PySTS.h"

namespace py = pybind11;
using namespace roboticslab;

void init_solvers(pybind11::module &m)
{
    py::module::import("PyKDL");
    py::class_<MatrixExponential> exp(m, "MatrixExponential");
    py::enum_<MatrixExponential::motion> motion_type(exp, "motion");
    motion_type.value("ROTATION", MatrixExponential::motion::ROTATION);
    motion_type.value("TRANSLATION", MatrixExponential::motion::TRANSLATION);

    exp.def(py::init<MatrixExponential::motion, const KDL::Vector &, const KDL::Vector &>(),
            py::arg("motionType"), py::arg("axis"), py::arg_v("origin", KDL::Vector::Zero(), "Vector.Zero"));

    exp.def("asFrame", &MatrixExponential::asFrame, py::arg("theta"));
    exp.def("getMotionType", &MatrixExponential::getMotionType);
    exp.def("getAxis", &MatrixExponential::getAxis);
    exp.def("getOrigin", &MatrixExponential::getOrigin);
    exp.def("changeBase", &MatrixExponential::changeBase, py::arg("H_new_old"));
    exp.def("cloneWithBase", &MatrixExponential::cloneWithBase, py::arg("H_new_old"));

    py::class_<PoeExpression> poe(m, "PoeExpression");
    poe.def(py::init<const KDL::Frame &>(), py::arg_v("H_S_T", KDL::Frame::Identity(), "Frame.Identity"));
    poe.def("append", static_cast<void (PoeExpression::*)(const MatrixExponential &, const KDL::Frame &)>(&PoeExpression::append),
            py::arg("exp"), py::arg_v("H_new_old", KDL::Frame::Identity(), "Frame.Identity"));
    poe.def("append", static_cast<void (PoeExpression::*)(const PoeExpression &, const KDL::Frame &)>(&PoeExpression::append),
            py::arg("poe"), py::arg_v("H_new_old", KDL::Frame::Identity(), "Frame.Identity"));
    poe.def("getTransform", &PoeExpression::getTransform);
    poe.def("size", &PoeExpression::size);
    poe.def("exponentialAtJoint", &PoeExpression::exponentialAtJoint, py::arg("i"));
    poe.def("changeBaseFrame", &PoeExpression::changeBaseFrame, py::arg("H_new_old"));
    poe.def("changeToolFrame", &PoeExpression::changeToolFrame, py::arg("H_new_old"));
    poe.def("evaluate", &PoeExpression::evaluate, py::arg("q"), py::arg("H"));
    poe.def("reverseSelf", &PoeExpression::reverseSelf);
    poe.def("makeReverse", &PoeExpression::makeReverse);
    poe.def("toChain", &PoeExpression::toChain);
    poe.def_static("fromChain", &PoeExpression::fromChain, py::arg("chain"));

    py::class_<ScrewTheoryIkSubproblem> ik_subproblem(m, "ScrewTheoryIkSubproblem");
    ik_subproblem.def("solve", static_cast<bool (ScrewTheoryIkSubproblem::*)(const KDL::Frame &, const KDL::Frame &, const ScrewTheoryIkSubproblem::JointConfig &, ScrewTheoryIkSubproblem::Solutions &) const>(&ScrewTheoryIkSubproblem::solve),
                      py::arg("rhs"), py::arg("pointTransform"), py::arg("reference"), py::arg("solutions"));
    ik_subproblem.def("solve", static_cast<bool (ScrewTheoryIkSubproblem::*)(const KDL::Frame &, const KDL::Frame &, ScrewTheoryIkSubproblem::Solutions &) const>(&ScrewTheoryIkSubproblem::solve),
                      py::arg("rhs"), py::arg("pointTransform"), py::arg("solutions"));
    ik_subproblem.def("solve", static_cast<bool (ScrewTheoryIkSubproblem::*)(const KDL::Frame &, const KDL::Frame &, const ScrewTheoryIkSubproblem::JointConfig &, ScrewTheoryIkSubproblem::Solutions &, const KDL::Frame &, const KDL::JntArray &) const>(&ScrewTheoryIkSubproblem::solve),
                      py::arg("rhs"), py::arg("pointTransform"), py::arg("reference"), py::arg("solutions"), py::arg("H_S_T_0"), py::arg("c_solutions"));
    ik_subproblem.def("solutions", &ScrewTheoryIkSubproblem::solutions);
    ik_subproblem.def("describe", &ScrewTheoryIkSubproblem::describe);

    py::class_<ScrewTheoryIkProblem, std::shared_ptr<ScrewTheoryIkProblem>> ik_problem(m, "ScrewTheoryIkProblem");
    // ik_problem.def("solve", static_cast<std::vector<bool> (ScrewTheoryIkProblem::*)(const KDL::Frame &, const KDL::JntArray &, ScrewTheoryIkProblem::Solutions &)>(&ScrewTheoryIkProblem::solve),
    //                py::arg("H_S_T"), py::arg("reference"), py::arg("solutions"));
    ik_problem.def("solve", [](ScrewTheoryIkProblem& self, const KDL::Frame & H_S_T)
    {
        ScrewTheoryIkProblem::Solutions solutions;
        auto reachability = self.solve(H_S_T, solutions);
        return std::make_tuple(solutions, reachability);
    });
    ik_problem.def("solutions", &ScrewTheoryIkProblem::solutions);
    ik_problem.def("getSteps", &ScrewTheoryIkProblem::getSteps);
    ik_problem.def("isReversed", &ScrewTheoryIkProblem::isReversed);
    ik_problem.def_static("create", &ScrewTheoryIkProblem::create, py::arg("poe"), py::arg("steps"), py::arg_v("reversed", false, "False"));

    py::class_<ScrewTheoryIkProblemBuilder> ik_problem_builder(m, "ScrewTheoryIkProblemBuilder");
    ik_problem_builder.def(py::init<const PoeExpression &>(), py::arg("poe"));
    ik_problem_builder.def("build", &ScrewTheoryIkProblemBuilder::build);

    py::class_<PadenKahanOne> pk1(m, "PadenKahanOne");
    pk1.def(py::init<const MatrixExponential &, const KDL::Vector &>(), py::arg("exp"), py::arg("p"));
    pk1.def("solve", static_cast<bool (PadenKahanOne::*)(const KDL::Frame &, const KDL::Frame &, const ScrewTheoryIkSubproblem::JointConfig &, ScrewTheoryIkSubproblem::Solutions &) const>(&PadenKahanOne::solve),
            py::arg("rhs"), py::arg("pointTransform"), py::arg("reference"), py::arg("solutions"));
    pk1.def("solutions", &PadenKahanOne::solutions);
    pk1.def("describe", &PadenKahanOne::describe);

    py::class_<PadenKahanTwo> pk2(m, "PadenKahanTwo");
    pk2.def(py::init<const MatrixExponential &, const MatrixExponential &, const KDL::Vector &, const KDL::Vector &>(),
            py::arg("exp1"), py::arg("exp2"), py::arg("p"), py::arg("r"));
    pk2.def("solve", static_cast<bool (PadenKahanTwo::*)(const KDL::Frame &, const KDL::Frame &, const ScrewTheoryIkSubproblem::JointConfig &, ScrewTheoryIkSubproblem::Solutions &) const>(&PadenKahanTwo::solve),
            py::arg("rhs"), py::arg("pointTransform"), py::arg("reference"), py::arg("solutions"));
    pk2.def("solutions", &PadenKahanTwo::solutions);
    pk2.def("describe", &PadenKahanTwo::describe);

    py::class_<PadenKahanThree> pk3(m, "PadenKahanThree");
    pk3.def(py::init<const MatrixExponential &, const KDL::Vector &, const KDL::Vector &>(),
            py::arg("exp"), py::arg("p"), py::arg("k"));
    pk3.def("solve", static_cast<bool (PadenKahanThree::*)(const KDL::Frame &, const KDL::Frame &, const ScrewTheoryIkSubproblem::JointConfig &, ScrewTheoryIkSubproblem::Solutions &) const>(&PadenKahanThree::solve),
            py::arg("rhs"), py::arg("pointTransform"), py::arg("reference"), py::arg("solutions"));
    pk3.def("solutions", &PadenKahanThree::solutions);
    pk3.def("describe", &PadenKahanThree::describe);

    py::class_<PardosGotorOne> pg1(m, "PardosGotorOne");
    pg1.def(py::init<const MatrixExponential &, const KDL::Vector &>(),
            py::arg("exp"), py::arg("p"));
    pg1.def("solve", static_cast<bool (PardosGotorOne::*)(const KDL::Frame &, const KDL::Frame &, const ScrewTheoryIkSubproblem::JointConfig &, ScrewTheoryIkSubproblem::Solutions &) const>(&PardosGotorOne::solve),
            py::arg("rhs"), py::arg("pointTransform"), py::arg("reference"), py::arg("solutions"));
    pg1.def("solutions", &PardosGotorOne::solutions);
    pg1.def("describe", &PardosGotorOne::describe);

    py::class_<PardosGotorTwo> pg2(m, "PardosGotorTwo");
    pg2.def(py::init<const MatrixExponential &, const MatrixExponential &, const KDL::Vector &>(),
            py::arg("exp1"), py::arg("exp2"), py::arg("p"));
    pg2.def("solve", static_cast<bool (PardosGotorTwo::*)(const KDL::Frame &, const KDL::Frame &, const ScrewTheoryIkSubproblem::JointConfig &, ScrewTheoryIkSubproblem::Solutions &) const>(&PardosGotorTwo::solve),
            py::arg("rhs"), py::arg("pointTransform"), py::arg("reference"), py::arg("solutions"));
    pg2.def("solutions", &PardosGotorTwo::solutions);
    pg2.def("describe", &PardosGotorTwo::describe);

    py::class_<PardosGotorThree> pg3(m, "PardosGotorThree");
    pg3.def(py::init<const MatrixExponential &, const KDL::Vector &, const KDL::Vector &>(),
            py::arg("exp1"), py::arg("p"), py::arg("k"));
    pg3.def("solve", static_cast<bool (PardosGotorThree::*)(const KDL::Frame &, const KDL::Frame &, const ScrewTheoryIkSubproblem::JointConfig &, ScrewTheoryIkSubproblem::Solutions &) const>(&PardosGotorThree::solve),
            py::arg("rhs"), py::arg("pointTransform"), py::arg("reference"), py::arg("solutions"));
    pg3.def("solutions", &PardosGotorThree::solutions);
    pg3.def("describe", &PardosGotorThree::describe);

    py::class_<PardosGotorThreePadenKahanOne> pg3_2(m, "PardosGotorThreePadenKahanOne");
    pg3_2.def(py::init<const MatrixExponential &, const MatrixExponential &, const KDL::Vector &, const KDL::Vector &>(),
              py::arg("exp"), py::arg("exp_pk1"), py::arg("p"), py::arg("k"));
    pg3_2.def("solve", static_cast<bool (PardosGotorThreePadenKahanOne::*)(const KDL::Frame &, const KDL::Frame &, const ScrewTheoryIkSubproblem::JointConfig &, ScrewTheoryIkSubproblem::Solutions &) const>(&PardosGotorThreePadenKahanOne::solve),
              py::arg("rhs"), py::arg("pointTransform"), py::arg("reference"), py::arg("solutions"));
    pg3_2.def("solutions", &PardosGotorThreePadenKahanOne::solutions);
    pg3_2.def("describe", &PardosGotorThreePadenKahanOne::describe);

    py::class_<PardosGotorFour> pg4(m, "PardosGotorFour");
    pg4.def(py::init<const MatrixExponential &, const MatrixExponential &, const KDL::Vector &>(),
            py::arg("exp1"), py::arg("exp2"), py::arg("p"));
    pg4.def("solve", static_cast<bool (PardosGotorFour::*)(const KDL::Frame &, const KDL::Frame &, const ScrewTheoryIkSubproblem::JointConfig &, ScrewTheoryIkSubproblem::Solutions &) const>(&PardosGotorFour::solve),
            py::arg("rhs"), py::arg("pointTransform"), py::arg("reference"), py::arg("solutions"));
    pg4.def("solutions", &PardosGotorFour::solutions);
    pg4.def("describe", &PardosGotorFour::describe);

    py::class_<PardosGotorFive> pg5(m, "PardosGotorFive");
    pg5.def(py::init<const MatrixExponential &, const MatrixExponential &, const KDL::Vector &>(),
            py::arg("exp"), py::arg("exp_next"), py::arg("p"));
    pg5.def("solve", static_cast<bool (PardosGotorFive::*)(const KDL::Frame &, const KDL::Frame &, const ScrewTheoryIkSubproblem::JointConfig &, ScrewTheoryIkSubproblem::Solutions &) const>(&PardosGotorFive::solve),
            py::arg("rhs"), py::arg("pointTransform"), py::arg("reference"), py::arg("solutions"));
    pg5.def("solutions", &PardosGotorFive::solutions);
    pg5.def("describe", &PardosGotorFive::describe);

    py::class_<PardosGotorSix> pg6(m, "PardosGotorSix");
    pg6.def(py::init<const MatrixExponential &, const MatrixExponential &, const KDL::Vector &>(),
            py::arg("exp1"), py::arg("exp2"), py::arg("p"));
    pg6.def("solve", static_cast<bool (PardosGotorSix::*)(const KDL::Frame &, const KDL::Frame &, const ScrewTheoryIkSubproblem::JointConfig &, ScrewTheoryIkSubproblem::Solutions &) const>(&PardosGotorSix::solve),
            py::arg("rhs"), py::arg("pointTransform"), py::arg("reference"), py::arg("solutions"));
    pg6.def("solutions", &PardosGotorSix::solutions);
    pg6.def("describe", &PardosGotorSix::describe);

    py::class_<PardosGotorSeven> pg7(m, "PardosGotorSeven");
    pg7.def(py::init<const MatrixExponential &, const MatrixExponential &, const MatrixExponential &, const KDL::Vector &>(),
            py::arg("exp1"), py::arg("exp2"), py::arg("exp3"), py::arg("p"));
    pg7.def("solve", static_cast<bool (PardosGotorSeven::*)(const KDL::Frame &, const KDL::Frame &, const ScrewTheoryIkSubproblem::JointConfig &, ScrewTheoryIkSubproblem::Solutions &) const>(&PardosGotorSeven::solve),
            py::arg("rhs"), py::arg("pointTransform"), py::arg("reference"), py::arg("solutions"));
    pg7.def("solutions", &PardosGotorSeven::solutions);
    pg7.def("describe", &PardosGotorSeven::describe);

    py::class_<Algebraic_UR> algebraic_ur(m, "Algebraic_UR");
    algebraic_ur.def(py::init<const int &, const int &>(), py::arg("j1"), py::arg("j2"));
    algebraic_ur.def("solve", static_cast<bool (Algebraic_UR::*)(const KDL::Frame &, const KDL::Frame &, const ScrewTheoryIkSubproblem::JointConfig &, ScrewTheoryIkSubproblem::Solutions &) const>(&Algebraic_UR::solve),
                     py::arg("rhs"), py::arg("pointTransform"), py::arg("reference"), py::arg("solutions"));
    algebraic_ur.def("solutions", &Algebraic_UR::solutions);
    algebraic_ur.def("describe", &Algebraic_UR::describe);

    py::class_<ConfigurationSelector> config_selector(m, "ConfigurationSelector");
    config_selector.def("configure", &ConfigurationSelector::configure, py::arg("solutions"), py::arg("reachability"));
    config_selector.def("findOptimalConfiguration", &ConfigurationSelector::findOptimalConfiguration, py::arg("qGuess"));
    config_selector.def("retrievePose", &ConfigurationSelector::retrievePose, py::arg("q"));
    config_selector.def("getValidSolutionIndex", &ConfigurationSelector::getValidSolutionIndex);
    // config_selector.def_static("INVALID_CONFIG", &ConfigurationSelector::INVALID_CONFIG);

    py::class_<ConfigurationSelectorLeastOverallAngularDisplacement> config_selector_load(m, "ConfigurationSelectorLeastOverallAngularDisplacement");
    config_selector_load.def(py::init<const KDL::JntArray &, const KDL::JntArray &>(), py::arg("qMin"), py::arg("qMax"));
    config_selector_load.def("configure", &ConfigurationSelectorLeastOverallAngularDisplacement::configure, py::arg("solutions"), py::arg("reachability"));
    config_selector_load.def("findOptimalConfiguration", &ConfigurationSelectorLeastOverallAngularDisplacement::findOptimalConfiguration, py::arg("qGuess"));
    // config_selector_load.def("retrievePose", &ConfigurationSelectorLeastOverallAngularDisplacement::retrievePose, py::arg("q"));
    config_selector_load.def("retrievePose", [](ConfigurationSelectorLeastOverallAngularDisplacement& self)
    {
        KDL::JntArray q;
        self.retrievePose(q);
        return q;
    });
    config_selector_load.def("getValidSolutionIndex", &ConfigurationSelectorLeastOverallAngularDisplacement::getValidSolutionIndex);
}
