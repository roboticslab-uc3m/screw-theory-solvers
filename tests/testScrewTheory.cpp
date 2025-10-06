#include "gtest/gtest.h"

#include <algorithm> // std::all_of

#include <kdl/chain.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/frames.hpp>
#include <kdl/jntarray.hpp>
#include <kdl/joint.hpp>
#include <kdl/utilities/utility.h>

#include "ConfigurationSelector.hpp"
#include "MatrixExponential.hpp"
#include "ProductOfExponentials.hpp"
#include "ScrewTheoryIkProblem.hpp"
#include "ScrewTheoryIkSubproblems.hpp"

namespace roboticslab::test
{

/**
 * @ingroup kinematics-dynamics-tests
 * @brief Tests classes related to @ref ScrewTheoryLib.
 */
class ScrewTheoryTest : public testing::Test
{
public:
    void SetUp() override
    {}

    void TearDown() override
    {}

    static KDL::JntArray fillJointValues(int size, double value)
    {
        KDL::JntArray q(size);

        for (int i = 0; i < size; i++)
        {
            q(i) = value;
        }

        return q;
    }

    static KDL::Chain makeTeoRightArmKinematicsFromDH()
    {
        const KDL::Joint rotZ(KDL::Joint::RotZ);
        KDL::Chain chain;

        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(    0, -KDL::PI_2,        0,          0)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(    0, -KDL::PI_2,        0, -KDL::PI_2)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(    0, -KDL::PI_2, -0.32901, -KDL::PI_2)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(    0,  KDL::PI_2,        0,          0)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(    0, -KDL::PI_2,   -0.215,          0)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(-0.09,          0,        0, -KDL::PI_2)));

        return chain;
    }

    static PoeExpression makeTeoRightArmKinematicsFromPoE()
    {
        KDL::Frame H_S_T({-0.63401, 0, 0});
        PoeExpression poe(H_S_T);

        poe.append(MatrixExponential(MatrixExponential::ROTATION, {0, 0, 1}, KDL::Vector::Zero()));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {0, 1, 0}, KDL::Vector::Zero()));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {1, 0, 0}, KDL::Vector::Zero()));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {0, 0, 1}, {-0.32901, 0, 0}));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {1, 0, 0}, {-0.32901, 0, 0}));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {0, 0, 1}, {-0.54401, 0, 0}));

        return poe;
    }

    static KDL::Chain makeTeoRightLegKinematicsFromDH()
    {
        const KDL::Joint rotZ(KDL::Joint::RotZ);
        KDL::Chain chain;

        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(        0,  KDL::PI_2,      0,         0)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(        0,  KDL::PI_2,      0, KDL::PI_2)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(    -0.33,          0,      0,         0)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(     -0.3,          0, 0.0175,         0)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(        0, -KDL::PI_2,      0,         0)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(-0.123005,          0,      0,         0)));

        return chain;
    }

    static PoeExpression makeTeoRightLegKinematicsFromPoE()
    {
        KDL::Frame H_S_T(KDL::Rotation::RotY(-KDL::PI_2) * KDL::Rotation::RotX(KDL::PI_2), {0.0175, 0, -0.753005});
        PoeExpression poe(H_S_T);

        poe.append(MatrixExponential(MatrixExponential::ROTATION, {0,  0, 1}, KDL::Vector::Zero()));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {0, -1, 0}, KDL::Vector::Zero()));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {1,  0, 0}, KDL::Vector::Zero()));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {1,  0, 0}, {     0, 0, -0.33}));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {1,  0, 0}, {0.0175, 0, -0.63}));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {0, -1, 0}, {0.0175, 0, -0.63}));

        return poe;
    }

    static KDL::Chain makeAbbIrb120KinematicsFromDH()
    {
        const KDL::Joint rotZ(KDL::Joint::RotZ);
        KDL::Chain chain;

        chain.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::None), KDL::Frame(KDL::Rotation::RotX(-KDL::PI_2))));

        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(   0,  KDL::PI_2,  0.29,          0)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(0.27,          0,     0,  KDL::PI_2)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(0.07,  KDL::PI_2,     0,          0)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(   0, -KDL::PI_2, 0.302,          0)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(   0,  KDL::PI_2,     0, -KDL::PI_2)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(   0,          0,  0.16,          0)));

        chain.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::None), KDL::Frame(KDL::Rotation::RotZ(KDL::PI_2))));

        return chain;
    }

    static PoeExpression makeAbbIrb120KinematicsFromPoE()
    {
        KDL::Frame H_S_T(KDL::Rotation::RotX(KDL::PI_2) * KDL::Rotation::RotZ(KDL::PI_2), {0.302, 0.47, 0});
        PoeExpression poe(H_S_T);

        poe.append(MatrixExponential(MatrixExponential::ROTATION, {0,  1, 0}, KDL::Vector::Zero()));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {0,  0, 1}, {    0, 0.29, 0}));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {0,  0, 1}, {    0, 0.56, 0}));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {1,  0, 0}, {0.302, 0.63, 0}));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {0,  0, 1}, {0.302, 0.63, 0}));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {0, -1, 0}, {0.302, 0.63, 0}));

        return poe;
    }

    static KDL::Chain makePumaKinematicsFromDH()
    {
        const KDL::Joint rotZ(KDL::Joint::RotZ);
        KDL::Chain chain;

        chain.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::None), KDL::Frame(KDL::Rotation::RotX(-KDL::PI_2))));

        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(0,  KDL::PI_2, 2,  KDL::PI_2)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(1,          0, 0,  KDL::PI_2)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(0, -KDL::PI_2, 0, -KDL::PI_2)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(0,  KDL::PI_2, 2,          0)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(0, -KDL::PI_2, 0,  KDL::PI_2)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(0,          0, 1,          0)));

        chain.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::None), KDL::Frame(KDL::Rotation::RotZ(-KDL::PI_2))));

        return chain;
    }

    static PoeExpression makePumaKinematicsFromPoE()
    {
        KDL::Frame H_S_T({0, 5, 1});
        PoeExpression poe(H_S_T);

        poe.append(MatrixExponential(MatrixExponential::ROTATION, {0, 1, 0}, KDL::Vector::Zero()));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {1, 0, 0}, {0, 2, 0}));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {1, 0, 0}, {0, 3, 0}));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {0, 1, 0}, {0, 3, 0}));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {1, 0, 0}, {0, 5, 0}));
        poe.append(MatrixExponential(MatrixExponential::ROTATION, {0, 0, 1}, {0, 5, 0}));

        return poe;
    }

    static KDL::Chain makeStanfordKinematicsFromDH()
    {
        const KDL::Joint rotZ(KDL::Joint::RotZ);
        const KDL::Joint translZ(KDL::Joint::TransZ);

        KDL::Chain chain;

        chain.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::None), KDL::Frame(KDL::Rotation::RotX(-KDL::PI_2))));

        chain.addSegment(KDL::Segment(   rotZ, KDL::Frame::DH(0,  KDL::PI_2, 2,  KDL::PI_2)));
        chain.addSegment(KDL::Segment(   rotZ, KDL::Frame::DH(0, -KDL::PI_2, 0,          0)));
        chain.addSegment(KDL::Segment(translZ, KDL::Frame::DH(0,          0, 0,          0)));
        chain.addSegment(KDL::Segment(   rotZ, KDL::Frame::DH(0,  KDL::PI_2, 3,          0)));
        chain.addSegment(KDL::Segment(   rotZ, KDL::Frame::DH(0, -KDL::PI_2, 0,  KDL::PI_2)));
        chain.addSegment(KDL::Segment(   rotZ, KDL::Frame::DH(0,          0, 1,          0)));

        chain.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::None), KDL::Frame(KDL::Rotation::RotZ(-KDL::PI_2))));

        return chain;
    }

    static PoeExpression makeStanfordKinematicsFromPoE()
    {
        KDL::Frame H_S_T({0, 5, 1});
        PoeExpression poe(H_S_T);

        poe.append(MatrixExponential(   MatrixExponential::ROTATION, {0, 1, 0}, KDL::Vector::Zero()));
        poe.append(MatrixExponential(   MatrixExponential::ROTATION, {1, 0, 0}, {0, 2, 0}));
        poe.append(MatrixExponential(MatrixExponential::TRANSLATION, {0, 1, 0}));
        poe.append(MatrixExponential(   MatrixExponential::ROTATION, {0, 1, 0}, {0, 2, 0}));
        poe.append(MatrixExponential(   MatrixExponential::ROTATION, {1, 0, 0}, {0, 5, 0}));
        poe.append(MatrixExponential(   MatrixExponential::ROTATION, {0, 0, 1}, {0, 5, 0}));

        return poe;
    }

    static KDL::Chain makeAbbIrb910scKinematicsFromDH()
    {
        const KDL::Joint rotZ(KDL::Joint::RotZ);
        const KDL::Joint translZ(KDL::Joint::TransZ);

        KDL::Chain chain;

        chain.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::None), KDL::Frame(KDL::Rotation::RotX(-KDL::PI_2))));

        chain.addSegment(KDL::Segment(   rotZ, KDL::Frame::DH( 0.4,       0,      0, 0)));
        chain.addSegment(KDL::Segment(   rotZ, KDL::Frame::DH(0.25, KDL::PI,      0, 0)));
        chain.addSegment(KDL::Segment(translZ, KDL::Frame::DH(   0,       0,      0, 0)));
        chain.addSegment(KDL::Segment(   rotZ, KDL::Frame::DH(   0,       0, -0.125, 0)));

        return chain;
    }

    static PoeExpression makeAbbIrb910scKinematicsFromPoE()
    {
        KDL::Frame H_S_T(KDL::Rotation::RotX(KDL::PI_2), {0.65, 0.125, 0});
        PoeExpression poe(H_S_T);

        poe.append(MatrixExponential(   MatrixExponential::ROTATION, {0,  1, 0}, KDL::Vector::Zero()));
        poe.append(MatrixExponential(   MatrixExponential::ROTATION, {0,  1, 0}, {0.4, 0, 0}));
        poe.append(MatrixExponential(MatrixExponential::TRANSLATION, {0, -1, 0}));
        poe.append(MatrixExponential(   MatrixExponential::ROTATION, {0, -1, 0}, {0.65, 0, 0}));

        return poe;
    }

    static KDL::Chain makeAbbIrb6620lxFromDh()
    {
        const KDL::Joint rotZ(KDL::Joint::RotZ);
        const KDL::Joint translZ(KDL::Joint::TransZ);

        KDL::Chain chain;

        chain.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::None), KDL::Frame({1.088, 2.5, 0})));

        chain.addSegment(KDL::Segment(translZ, KDL::Frame::DH( 0.38,          0,     0,         0)));
        chain.addSegment(KDL::Segment(   rotZ, KDL::Frame::DH(0.975,          0,     0,         0)));
        chain.addSegment(KDL::Segment(   rotZ, KDL::Frame::DH(  0.2,  KDL::PI_2,     0,         0)));
        chain.addSegment(KDL::Segment(   rotZ, KDL::Frame::DH(    0, -KDL::PI_2, 0.887,         0)));
        chain.addSegment(KDL::Segment(   rotZ, KDL::Frame::DH(    0,  KDL::PI_2,     0, KDL::PI_2)));
        chain.addSegment(KDL::Segment(   rotZ, KDL::Frame::DH(    0,          0, 0.357,         0)));

        chain.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::None), KDL::Frame(KDL::Rotation::RotZ(-KDL::PI_2))));

        return chain;
    }

    static PoeExpression makeAbbIrb6620lxFromPoE()
    {
        KDL::Frame H_S_T(KDL::Rotation::RotY(KDL::PI_2), {3, 1.613, 0});
        PoeExpression poe(H_S_T);

        poe.append(MatrixExponential(MatrixExponential::TRANSLATION, {0,  0, 1}));
        poe.append(MatrixExponential(   MatrixExponential::ROTATION, {0,  0, 1}, {1.468,   2.5, 0}));
        poe.append(MatrixExponential(   MatrixExponential::ROTATION, {0,  0, 1}, {2.443,   2.5, 0}));
        poe.append(MatrixExponential(   MatrixExponential::ROTATION, {0, -1, 0}, {2.643, 1.613, 0}));
        poe.append(MatrixExponential(   MatrixExponential::ROTATION, {0,  0, 1}, {2.643, 1.613, 0}));
        poe.append(MatrixExponential(   MatrixExponential::ROTATION, {1,  0, 0}, {2.643, 1.613, 0}));

        return poe;
    }

    static KDL::Chain makeUR16eFromDh()
    {
        const KDL::Joint rotZ(KDL::Joint::RotZ);
        KDL::Chain chain;

        chain.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::None), KDL::Frame(KDL::Rotation::RotZ(0))));

        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(0.0, -KDL::PI_2, 0.181, 0)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(0.478, 0, 0.0, 0)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(0.36, 0, 0.0, 0)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(0.0, -KDL::PI_2, 0.174, 0)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(0.0, KDL::PI_2, 0.12, 0)));
        chain.addSegment(KDL::Segment(rotZ, KDL::Frame::DH(0.0, 0, 0.19, KDL::PI)));

        chain.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::None), KDL::Frame(KDL::Rotation::RotZ(0))));

        return chain;
    }

    static PoeExpression makeUR16eFromPoE()
    {
        KDL::Chain chain = makeUR16eFromDh();
        KDL::ChainFkSolverPos_recursive fkSolver(chain);
        KDL::JntArray q(chain.getNrOfJoints()); 

        KDL::Frame H_DH_0;
        fkSolver.JntToCart(q, H_DH_0);

        KDL::Frame H_S_T = H_DH_0;  
        PoeExpression poe(H_S_T);

        poe.append(MatrixExponential(   MatrixExponential::ROTATION, {0,  0, 1}, {    0,     0, 0.181}));
        poe.append(MatrixExponential(   MatrixExponential::ROTATION, {0,  1, 0}, {    0,     0, 0.181}));
        poe.append(MatrixExponential(   MatrixExponential::ROTATION, {0,  1, 0}, {0.478,     0, 0.181}));
        poe.append(MatrixExponential(   MatrixExponential::ROTATION, {0,  1, 0}, {0.838, 0.174, 0.181}));
        poe.append(MatrixExponential(   MatrixExponential::ROTATION, {0,  0,-1}, {0.838, 0.174, 0.061}));
        poe.append(MatrixExponential(   MatrixExponential::ROTATION, {0,  1, 0}, {0.838, 0.174, 0.061}));

        return poe;
    }

    static void checkSolutions(const ScrewTheoryIkSubproblem::Solutions & actual, const ScrewTheoryIkSubproblem::Solutions & expected)
    {
        ASSERT_EQ(actual.size(), expected.size());

        for (int i = 0; i < actual.size(); i++)
        {
            ASSERT_EQ(actual[i].size(), expected[i].size());

            for (int j = 0; j < actual[i].size(); j++)
            {
                ASSERT_NEAR(actual[i][j], expected[i][j], KDL::epsilon);
            }
        }
    }

    static void checkRobotKinematicsInternal(const KDL::Chain & chain, const PoeExpression & poe, const KDL::JntArray & q, int soln)
    {
        ASSERT_EQ(poe.size(), chain.getNrOfJoints());

        KDL::ChainFkSolverPos_recursive fkSolver(chain);
        KDL::Frame H_S_T_q_DH, H_S_T_q_ST;
        ASSERT_EQ(fkSolver.JntToCart(q, H_S_T_q_DH), KDL::SolverI::E_NOERROR);
        ASSERT_TRUE(poe.evaluate(q, H_S_T_q_ST));
        ASSERT_EQ(H_S_T_q_ST, H_S_T_q_DH);

        ScrewTheoryIkProblemBuilder builder(poe);
        ScrewTheoryIkProblem * ikProblem = builder.build();

        ASSERT_TRUE(ikProblem);
        ASSERT_EQ(ikProblem->solutions(), soln);
        ScrewTheoryIkProblem::Solutions solutions;
        auto reachability = ikProblem->solve(H_S_T_q_ST, q, solutions);

        ASSERT_TRUE(std::all_of(reachability.begin(), reachability.end(), [](bool r) { return r; }));
        delete ikProblem;

        bool match = false;

        for (auto j = 0; j < solutions.size(); j++)
        {
            const auto & solution = solutions[j];

            KDL::Frame H_S_T_q_ST_validate;
            ASSERT_TRUE(poe.evaluate(solution, H_S_T_q_ST_validate));
            ASSERT_EQ(H_S_T_q_ST_validate, H_S_T_q_ST);

            if (solution == q)
            {
                match = true;
            }
        }

        ASSERT_TRUE(match);
    }

    static void checkRobotKinematics(const KDL::Chain & chain, const PoeExpression & poe, int soln)
    {
        const int numJoints = chain.getNrOfJoints();

        checkRobotKinematicsInternal(chain, poe, fillJointValues(numJoints, 0.0), soln);
        checkRobotKinematicsInternal(chain, poe, fillJointValues(numJoints, 0.1), soln);
        checkRobotKinematicsInternal(chain, poe, fillJointValues(numJoints, KDL::PI_2), soln);
        checkRobotKinematicsInternal(chain, poe, fillJointValues(numJoints, KDL::PI), soln);
    }

    static int findTargetConfiguration(const ScrewTheoryIkProblem::Solutions & solutions, const KDL::JntArray & target)
    {
        for (int i = 0; i < solutions.size(); i++)
        {
            bool equal = true;

            for (int j = 0; j < solutions[i].rows(); j++)
            {
                if (!KDL::Equal(solutions[i](j), target(j)))
                {
                    equal = false;
                    break;
                }
            }

            if (equal)
            {
                return i;
            }
        }

        return -1;
    }
};

TEST_F(ScrewTheoryTest, MatrixExponentialInit)
{
    MatrixExponential::motion motionType = MatrixExponential::ROTATION;
    KDL::Vector axis(1, 0, 0);
    KDL::Vector origin(1, 1, 1);

    MatrixExponential exp(MatrixExponential::ROTATION, axis, origin);

    ASSERT_EQ(exp.getMotionType(), motionType);
    ASSERT_EQ(exp.getAxis(), axis);
    ASSERT_EQ(exp.getOrigin(), origin);

    KDL::Frame H_new_old(KDL::Rotation::RotZ(KDL::PI_2), {-1, 0, 0});
    MatrixExponential newExp = exp.cloneWithBase(H_new_old);

    ASSERT_EQ(newExp.getMotionType(), exp.getMotionType());
    ASSERT_EQ(newExp.getAxis(), KDL::Vector(0, 1, 0));
    ASSERT_EQ(newExp.getOrigin(), KDL::Vector(-2, 1, 1));

    exp.changeBase(H_new_old);

    ASSERT_EQ(exp.getMotionType(), newExp.getMotionType());
    ASSERT_EQ(exp.getAxis(), newExp.getAxis());
    ASSERT_EQ(exp.getOrigin(), newExp.getOrigin());
}

TEST_F(ScrewTheoryTest, MatrixExponentialRotation)
{
    double theta = KDL::PI_2;
    MatrixExponential exp(MatrixExponential::ROTATION, {0, 0, 1}, {1, 0, 0});

    KDL::Frame actual = exp.asFrame(theta);
    KDL::Frame expected(KDL::Rotation::RotZ(theta), {1, -1, 0});
    ASSERT_EQ(actual, expected);
}

TEST_F(ScrewTheoryTest, MatrixExponentialTranslation)
{
    KDL::Vector axis(0, 0, 1);
    double theta = KDL::PI;
    MatrixExponential exp(MatrixExponential::TRANSLATION, axis);

    KDL::Frame actual = exp.asFrame(theta);
    KDL::Frame expected(KDL::Rotation::Identity(), theta * axis);

    ASSERT_EQ(actual, expected);
}

TEST_F(ScrewTheoryTest, ProductOfExponentialsInit)
{
    PoeExpression poe;

    ASSERT_EQ(poe.size(), 0);
    ASSERT_EQ(poe.getTransform(), KDL::Frame::Identity());

    MatrixExponential exp(MatrixExponential::ROTATION, {1, 0, 0}, {1, 1, 1});
    poe.append(exp);

    ASSERT_EQ(poe.size(), 1);

    poe.changeBaseFrame(KDL::Frame(KDL::Rotation::RotZ(KDL::PI_2), {-1, 0, 0}));
    poe.changeToolFrame(KDL::Frame({0, 0, 1}));

    ASSERT_EQ(poe.getTransform(), KDL::Frame(KDL::Rotation::RotZ(KDL::PI_2), {-1, 0, 1}));

    KDL::Frame H;
    KDL::JntArray q(1);
    q(0) = KDL::PI_2;

    ASSERT_TRUE(poe.evaluate(q, H));
    ASSERT_EQ(H, KDL::Frame(KDL::Rotation(KDL::Rotation::RotZ(KDL::PI_2) * KDL::Rotation::RotX(KDL::PI_2)), {-2, 0, 0}));

    PoeExpression poe2(KDL::Frame(KDL::Rotation::RotX(KDL::PI_2), {0, 1, 0}));
    MatrixExponential exp2(MatrixExponential::ROTATION, {0, 1, 0}, {1, 0, 0});

    poe2.append(exp2);
    poe.append(poe2, KDL::Frame(KDL::Rotation::RotZ(KDL::PI_2), {0, 0, 1}));

    ASSERT_EQ(poe.size(), 2);
    ASSERT_EQ(poe.getTransform(), KDL::Frame(KDL::Rotation::RotZ(KDL::PI_2) * KDL::Rotation::RotX(KDL::PI_2), {-1, 0, 1}));

    q.resize(2);
    q(1) = -KDL::PI_2;

    ASSERT_TRUE(poe.evaluate(q, H));
    ASSERT_EQ(H, KDL::Frame(KDL::Rotation::RotX(KDL::PI), {-3, 1, 0}));
}

TEST_F(ScrewTheoryTest, ProductOfExponentialsFromChain)
{
    KDL::Chain chain = makeTeoRightArmKinematicsFromDH();
    PoeExpression poe = PoeExpression::fromChain(chain);
    ASSERT_EQ(poe.size(), chain.getNrOfJoints());

    KDL::ChainFkSolverPos_recursive fkSolver(chain);
    KDL::Frame H_S_T_0;

    ASSERT_EQ(fkSolver.JntToCart(KDL::JntArray(chain.getNrOfJoints()), H_S_T_0), KDL::SolverI::E_NOERROR);
    ASSERT_EQ(poe.getTransform(), H_S_T_0);

    KDL::JntArray q = fillJointValues(chain.getNrOfJoints(), KDL::PI_2);;
    KDL::Frame H_S_T_q_DH, H_S_T_q_ST;

    ASSERT_EQ(fkSolver.JntToCart(q, H_S_T_q_DH), KDL::SolverI::E_NOERROR);
    ASSERT_TRUE(poe.evaluate(q, H_S_T_q_ST));
    ASSERT_EQ(H_S_T_q_ST, H_S_T_q_DH);
}

TEST_F(ScrewTheoryTest, ProductOfExponentialsToChain)
{
    PoeExpression poe = makeTeoRightArmKinematicsFromPoE();

    KDL::Chain chainDH = makeTeoRightArmKinematicsFromDH();
    KDL::Chain chainST = poe.toChain();

    KDL::ChainFkSolverPos_recursive fkSolverDH(chainDH);
    KDL::ChainFkSolverPos_recursive fkSolverST(chainST);

    KDL::JntArray q = fillJointValues(poe.size(), KDL::PI_2);;
    KDL::Frame H_S_T_q_DH, H_S_T_q_ST;

    ASSERT_EQ(fkSolverDH.JntToCart(q, H_S_T_q_DH), KDL::SolverI::E_NOERROR);
    ASSERT_EQ(fkSolverST.JntToCart(q, H_S_T_q_ST), KDL::SolverI::E_NOERROR);
    ASSERT_EQ(H_S_T_q_ST, H_S_T_q_DH);
}

TEST_F(ScrewTheoryTest, ProductOfExponentialsIntegrity)
{
    PoeExpression poe = makeTeoRightArmKinematicsFromPoE();

    KDL::Chain chain = poe.toChain();
    KDL::ChainFkSolverPos_recursive fkSolver(chain);

    KDL::JntArray q = fillJointValues(poe.size(), KDL::PI_2);;
    KDL::Frame H_S_T_q_DH, H_S_T_q_ST;

    ASSERT_EQ(fkSolver.JntToCart(q, H_S_T_q_DH), KDL::SolverI::E_NOERROR);

    PoeExpression poeFromChain = PoeExpression::fromChain(chain);

    ASSERT_TRUE(poeFromChain.evaluate(q, H_S_T_q_ST));
    ASSERT_EQ(H_S_T_q_ST, H_S_T_q_DH);
}

TEST_F(ScrewTheoryTest, ProductOfExponentialsReverse)
{
    PoeExpression poe = makeTeoRightArmKinematicsFromPoE();
    PoeExpression poeReversed = poe.makeReverse();

    ASSERT_EQ(poeReversed.size(), poe.size());
    ASSERT_EQ(poeReversed.getTransform(), poe.getTransform().Inverse());

    KDL::Frame H_S_T_0, H_S_T_0_reversed;

    ASSERT_TRUE(poe.evaluate(KDL::JntArray(poe.size()), H_S_T_0));
    ASSERT_TRUE(poeReversed.evaluate(KDL::JntArray(poeReversed.size()), H_S_T_0_reversed));
    ASSERT_EQ(H_S_T_0_reversed, H_S_T_0.Inverse());

    KDL::JntArray q = fillJointValues(poe.size(), KDL::PI_2);;
    KDL::Frame H_S_T_q;

    ASSERT_TRUE(poe.evaluate(q, H_S_T_q));

    q = fillJointValues(poe.size(), -KDL::PI_2);

    KDL::Frame H_S_T_q_reversed;

    ASSERT_TRUE(poeReversed.evaluate(q, H_S_T_q_reversed));
    ASSERT_EQ(H_S_T_q_reversed, H_S_T_q.Inverse());

    poe.reverseSelf();

    ASSERT_TRUE(poe.evaluate(q, H_S_T_q_reversed));
    ASSERT_EQ(H_S_T_q_reversed, H_S_T_q.Inverse());
}

TEST_F(ScrewTheoryTest, PadenKahanOne)
{
    KDL::Vector p(0, 1, 0);
    KDL::Vector k(1, 1, 1);

    MatrixExponential exp(MatrixExponential::ROTATION, {0, 1, 0}, {1, 0, 0});
    PadenKahanOne pk1(exp, p);

    ASSERT_EQ(pk1.solutions(), 1);

    KDL::Frame rhs(k - p);
    ScrewTheoryIkSubproblem::Solutions actual;
    ASSERT_TRUE(pk1.solve(rhs, KDL::Frame::Identity(), actual));

    ASSERT_EQ(actual.size(), 1);
    ASSERT_EQ(actual[0].size(), 1);

    ScrewTheoryIkSubproblem::Solutions expected = {{KDL::PI_2}};

    checkSolutions(actual, expected);

    KDL::Vector k2 = k + KDL::Vector(0, 0, 1);
    KDL::Frame rhs2(k2 - p);
    ASSERT_FALSE(pk1.solve(rhs2, KDL::Frame::Identity(), actual));

    checkSolutions(actual, expected);

    KDL::Vector k3 = k + KDL::Vector(0, 1, 0);
    KDL::Frame rhs3(k3 - p);
    ASSERT_FALSE(pk1.solve(rhs3, KDL::Frame::Identity(), actual));

    checkSolutions(actual, expected);

    KDL::Vector k4 = k + KDL::Vector(0, 0, -1);
    KDL::Frame rhs4(k4 - p);
    ASSERT_FALSE(pk1.solve(rhs4, KDL::Frame::Identity(), actual));

    expected = {{0}};
    checkSolutions(actual, expected);

    expected = {{KDL::PI}};
    ASSERT_FALSE(pk1.solve(rhs4, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Vector p2 = p + KDL::Vector(1, 0, 0);
    PadenKahanOne pk1b(exp, p2);
    KDL::Frame rhs5(k - p2);
    ASSERT_FALSE(pk1b.solve(rhs5, KDL::Frame::Identity(), actual));

    expected = {{0}};
    checkSolutions(actual, expected);

    expected = {{KDL::PI}};
    ASSERT_FALSE(pk1b.solve(rhs5, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Frame rhs6(k4 - p2);
    ASSERT_TRUE(pk1b.solve(rhs6, KDL::Frame::Identity(), actual));

    expected = {{0}};
    checkSolutions(actual, expected);

    expected = {{KDL::PI}};
    ASSERT_TRUE(pk1b.solve(rhs6, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Vector k7 = k4 + KDL::Vector(0, 1, 0);
    KDL::Frame rhs7(k7 - p2);
    ASSERT_FALSE(pk1b.solve(rhs7, KDL::Frame::Identity(), actual));

    expected = {{0}};
    checkSolutions(actual, expected);

    expected = {{KDL::PI}};
    ASSERT_FALSE(pk1b.solve(rhs7, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);
}

TEST_F(ScrewTheoryTest, PadenKahanTwo)
{
    KDL::Vector p(0, 1, 0);
    KDL::Vector k(1, -1, 1);
    KDL::Vector r(1, 0, 0);

    MatrixExponential exp1(MatrixExponential::ROTATION, {1, 0, 0}, r);
    MatrixExponential exp2(MatrixExponential::ROTATION, {0, 1, 0}, r);
    PadenKahanTwo pk2(exp1, exp2, p, r);

    ASSERT_EQ(pk2.solutions(), 2);

    KDL::Frame rhs(k - p);
    ScrewTheoryIkSubproblem::Solutions actual;
    ASSERT_TRUE(pk2.solve(rhs, KDL::Frame::Identity(), actual));

    ASSERT_EQ(actual.size(), 2);
    ASSERT_EQ(actual[0].size(), 2);
    ASSERT_EQ(actual[1].size(), 2);

    ScrewTheoryIkSubproblem::Solutions expected = {
        {KDL::PI, -KDL::PI_2},
        {KDL::PI_2, KDL::PI_2}
    };

    checkSolutions(actual, expected);

    KDL::Vector k2 = k + KDL::Vector(1, 0, 0);
    KDL::Frame rhs2(k2 - p);
    ASSERT_FALSE(pk2.solve(rhs2, KDL::Frame::Identity(), actual));

    expected = {
        {3 * KDL::PI_4, KDL::PI},
        {3 * KDL::PI_4, KDL::PI}
    };

    checkSolutions(actual, expected);

    KDL::Vector p3 = p + KDL::Vector(0.5, 0, 0);
    PadenKahanTwo pk2c(exp1, exp2, p3, r);
    KDL::Frame rhs3(k2 - p3);
    ASSERT_FALSE(pk2c.solve(rhs3, KDL::Frame::Identity(), actual));

    checkSolutions(actual, expected);

    KDL::Vector p4 = p + KDL::Vector(1, std::sqrt(2) - 1, 0);
    PadenKahanTwo pk2d(exp1, exp2, p4, r);
    KDL::Frame rhs4(k - p4);
    ASSERT_TRUE(pk2d.solve(rhs4, KDL::Frame::Identity(), actual));

    expected = {
        {3 * KDL::PI_4, 0},
        {3 * KDL::PI_4, 0}
    };

    checkSolutions(actual, expected);

    expected = {
        {3 * KDL::PI_4, KDL::PI},
        {3 * KDL::PI_4, KDL::PI}
    };

    ASSERT_TRUE(pk2d.solve(rhs4, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Vector p5 = p + KDL::Vector(1, 0, 0);
    PadenKahanTwo pk2e(exp1, exp2, p5, r);
    KDL::Frame rhs5(k - p5);
    ASSERT_FALSE(pk2e.solve(rhs5, KDL::Frame::Identity(), actual));

    expected = {
        {3 * KDL::PI_4, 0},
        {3 * KDL::PI_4, 0}
    };

    checkSolutions(actual, expected);

    expected = {
        {3 * KDL::PI_4, KDL::PI},
        {3 * KDL::PI_4, KDL::PI}
    };

    ASSERT_FALSE(pk2e.solve(rhs5, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Vector p6 = p + KDL::Vector(0, -1, 0);
    KDL::Vector k6 = k + KDL::Vector(1, 1, -1);
    PadenKahanTwo pk2f(exp1, exp2, p6, r);
    KDL::Frame rhs6(k6 - p6);
    ASSERT_TRUE(pk2f.solve(rhs6, KDL::Frame::Identity(), actual));

    expected = {
        {0, KDL::PI},
        {0, KDL::PI}
    };

    checkSolutions(actual, expected);

    expected = {
        {KDL::PI, KDL::PI},
        {KDL::PI, KDL::PI}
    };

    ASSERT_TRUE(pk2f.solve(rhs6, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Vector k7 = k + KDL::Vector(2, 1, -1);
    KDL::Frame rhs7(k7 - p6);
    ASSERT_FALSE(pk2f.solve(rhs7, KDL::Frame::Identity(), actual));

    expected = {
        {0, KDL::PI},
        {0, KDL::PI}
    };

    checkSolutions(actual, expected);

    expected = {
        {KDL::PI, KDL::PI},
        {KDL::PI, KDL::PI}
    };

    ASSERT_FALSE(pk2f.solve(rhs7, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Vector p8 = r;
    KDL::Vector k8 = r;
    PadenKahanTwo pk2g(exp1, exp2, p8, r);
    KDL::Frame rhs8(k8 - p8);
    ASSERT_TRUE(pk2g.solve(rhs8, KDL::Frame::Identity(), actual));

    expected = {
        {0, 0},
        {0, 0}
    };

    checkSolutions(actual, expected);

    expected = {
        {KDL::PI, KDL::PI},
        {KDL::PI, KDL::PI}
    };

    ASSERT_TRUE(pk2g.solve(rhs8, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Vector p9 = p + KDL::Vector(1, 0, 0);
    KDL::Vector k9 = k + KDL::Vector(1, 1, -1);
    PadenKahanTwo pk2h(exp1, exp2, p9, r);
    KDL::Frame rhs9(k9 - p9);
    ASSERT_FALSE(pk2h.solve(rhs9, KDL::Frame::Identity(), actual));

    expected = {
        {0, 0},
        {0, 0}
    };

    checkSolutions(actual, expected);

    expected = {
        {KDL::PI, KDL::PI},
        {KDL::PI, KDL::PI}
    };

    ASSERT_FALSE(pk2h.solve(rhs9, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);
}

TEST_F(ScrewTheoryTest, PadenKahanThree)
{
    KDL::Vector p(0, 1, 0);
    KDL::Vector k(2, 1, 1);
    KDL::Vector delta(-1, 0, 0);

    MatrixExponential exp(MatrixExponential::ROTATION, {0, 1, 0}, {1, 0, 0});
    PadenKahanThree pk3(exp, p, k);

    ASSERT_EQ(pk3.solutions(), 2);

    KDL::Frame rhs(delta - (p - k));
    ScrewTheoryIkSubproblem::Solutions actual;
    ASSERT_TRUE(pk3.solve(rhs, KDL::Frame::Identity(), actual));

    ASSERT_EQ(actual.size(), 2);
    ASSERT_EQ(actual[0].size(), 1);
    ASSERT_EQ(actual[1].size(), 1);

    ScrewTheoryIkSubproblem::Solutions expected = {
        {KDL::PI},
        {KDL::PI_2}
    };

    checkSolutions(actual, expected);

    KDL::Vector k2 = k + KDL::Vector(1, 0, 1);
    PadenKahanThree pk3b(exp, p, k2);
    KDL::Frame rhs2(delta - (p - k2));

    ASSERT_FALSE(pk3b.solve(rhs2, KDL::Frame::Identity(), actual));

    expected = {
        {3 * KDL::PI_4},
        {3 * KDL::PI_4}
    };

    checkSolutions(actual, expected);

    KDL::Vector k3 = k + KDL::Vector(-1, 0, -1);
    PadenKahanThree pk3c(exp, p, k3);
    KDL::Frame rhs3(delta - (p - k3));

    ASSERT_TRUE(pk3c.solve(rhs3, KDL::Frame::Identity(), actual));

    expected = {{0}, {0}};
    checkSolutions(actual, expected);

    expected = {{KDL::PI}, {KDL::PI}};
    ASSERT_TRUE(pk3c.solve(rhs3, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Vector delta4(-1, -1, 0);
    KDL::Frame rhs4(delta4 - (p - k3));

    ASSERT_FALSE(pk3c.solve(rhs4, KDL::Frame::Identity(), actual));

    expected = {{0}, {0}};
    checkSolutions(actual, expected);

    expected = {{KDL::PI}, {KDL::PI}};
    ASSERT_FALSE(pk3c.solve(rhs4, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Vector p5 = p + KDL::Vector(1, 0, 0);
    PadenKahanThree pk3d(exp, p5, k);
    KDL::Frame rhs5(delta4 - (p5 - k));

    ASSERT_TRUE(pk3d.solve(rhs5, KDL::Frame::Identity(), actual));

    expected = {{0}, {0}};
    checkSolutions(actual, expected);

    expected = {{KDL::PI}, {KDL::PI}};
    ASSERT_TRUE(pk3d.solve(rhs5, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Frame rhs6(delta - (p5 - k));
    ASSERT_FALSE(pk3d.solve(rhs6, KDL::Frame::Identity(), actual));

    expected = {{0}, {0}};
    checkSolutions(actual, expected);

    expected = {{KDL::PI}, {KDL::PI}};
    ASSERT_FALSE(pk3d.solve(rhs6, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    PadenKahanThree pk3e(exp, p5, k3);
    KDL::Frame rhs7(KDL::Vector::Zero() - (p5 - k3));

    ASSERT_TRUE(pk3e.solve(rhs7, KDL::Frame::Identity(), actual));

    expected = {{0}, {0}};
    checkSolutions(actual, expected);

    expected = {{KDL::PI}, {KDL::PI}};
    ASSERT_TRUE(pk3e.solve(rhs7, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Frame rhs8(delta - (p5 - k3));
    ASSERT_FALSE(pk3e.solve(rhs8, KDL::Frame::Identity(), actual));

    expected = {{0}, {0}};
    checkSolutions(actual, expected);

    expected = {{KDL::PI}, {KDL::PI}};
    ASSERT_FALSE(pk3e.solve(rhs8, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Vector k9 = k3 + KDL::Vector(0, 1, 0);
    KDL::Vector delta9(0, -1, 0);
    KDL::Frame rhs9(delta9 - (p5 - k9));

    ASSERT_TRUE(pk3e.solve(rhs9, KDL::Frame::Identity(), actual));

    expected = {{0}, {0}};
    checkSolutions(actual, expected);

    expected = {{KDL::PI}, {KDL::PI}};
    ASSERT_TRUE(pk3e.solve(rhs9, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Frame rhs10(delta - (p5 - k9));

    ASSERT_FALSE(pk3e.solve(rhs10, KDL::Frame::Identity(), actual));

    expected = {{0}, {0}};
    checkSolutions(actual, expected);

    expected = {{KDL::PI}, {KDL::PI}};
    ASSERT_FALSE(pk3e.solve(rhs10, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Vector p11 = p + KDL::Vector(-1, 0, 0);
    KDL::Vector k11 = k + KDL::Vector(0, 0, -1);
    KDL::Vector delta11(3, 0, 0);
    PadenKahanThree pk3f(exp, p11, k11);
    KDL::Frame rhs11(delta11 - (p11 - k11));

    ASSERT_TRUE(pk3f.solve(rhs11, KDL::Frame::Identity(), actual));

    expected = {{0}, {0}};

    checkSolutions(actual, expected);

    KDL::Vector k12 = k + KDL::Vector(1, 0, -1);
    PadenKahanThree pk3g(exp, p, k12);
    KDL::Frame rhs12(delta - (p - k12));

    ASSERT_TRUE(pk3g.solve(rhs12, KDL::Frame::Identity(), actual));

    expected = {{KDL::PI}, {KDL::PI}};

    checkSolutions(actual, expected);
}

TEST_F(ScrewTheoryTest, PardosGotorOne)
{
    KDL::Vector p(1, 0, 0);
    KDL::Vector k(1, 1, 0);

    MatrixExponential exp(MatrixExponential::TRANSLATION, {0, 1, 0});
    PardosGotorOne pg1(exp, p);

    ASSERT_EQ(pg1.solutions(), 1);

    KDL::Frame rhs(k - p);
    ScrewTheoryIkSubproblem::Solutions actual;
    ASSERT_TRUE(pg1.solve(rhs, KDL::Frame::Identity(), actual));

    ASSERT_EQ(actual.size(), 1);
    ASSERT_EQ(actual[0].size(), 1);

    ScrewTheoryIkSubproblem::Solutions expected = {{1.0}};

    checkSolutions(actual, expected);
}

TEST_F(ScrewTheoryTest, PardosGotorTwo)
{
    KDL::Vector p(1, 1, 0);
    KDL::Vector k(2, 3, 0);

    MatrixExponential exp1(MatrixExponential::TRANSLATION, {0, 1, 0});
    MatrixExponential exp2(MatrixExponential::TRANSLATION, {1, 0, 0});
    PardosGotorTwo pg2(exp1, exp2, p);

    ASSERT_EQ(pg2.solutions(), 1);

    KDL::Frame rhs(k - p);
    ScrewTheoryIkSubproblem::Solutions actual;
    ASSERT_TRUE(pg2.solve(rhs, KDL::Frame::Identity(), actual));

    ASSERT_EQ(actual.size(), 1);
    ASSERT_EQ(actual[0].size(), 2);

    ScrewTheoryIkSubproblem::Solutions expected = {{k.y() - p.y(), k.x() - p.x()}};

    checkSolutions(actual, expected);
}

TEST_F(ScrewTheoryTest, PardosGotorThree)
{
    KDL::Vector p(1, 0, 0);
    KDL::Vector k(1, 2, 0);
    KDL::Vector delta(0, 1, 0);

    MatrixExponential exp(MatrixExponential::TRANSLATION, {0, 1, 0});
    PardosGotorThree pg3(exp, p, k);

    ASSERT_EQ(pg3.solutions(), 2);

    KDL::Frame rhs(delta - (p - k));
    ScrewTheoryIkSubproblem::Solutions actual;
    ASSERT_TRUE(pg3.solve(rhs, KDL::Frame::Identity(), actual));

    ASSERT_EQ(actual.size(), 2);
    ASSERT_EQ(actual[0].size(), 1);
    ASSERT_EQ(actual[1].size(), 1);

    ScrewTheoryIkSubproblem::Solutions expected = {
        {k.y() - p.y() + delta.y()},
        {k.y() - p.y() - delta.y()}
    };

    checkSolutions(actual, expected);

    KDL::Vector k2 = k + KDL::Vector(2, 0, 0);
    PardosGotorThree pg3b(exp, p, k2);
    KDL::Frame rhs2(delta - (p - k2));

    ASSERT_FALSE(pg3b.solve(rhs2, KDL::Frame::Identity(), actual));

    expected = {{2.0}, {2.0}};

    checkSolutions(actual, expected);
}

TEST_F(ScrewTheoryTest, PardosGotorFour)
{
    KDL::Vector p(0, 1, 0);
    KDL::Vector k(3, 1, 1);

    MatrixExponential exp1(MatrixExponential::ROTATION, {0, 1, 0}, {2, 0, 0});
    MatrixExponential exp2(MatrixExponential::ROTATION, {0, 1, 0}, {1, 0, 0});
    PardosGotorFour pg4(exp1, exp2, p);

    ASSERT_EQ(pg4.solutions(), 2);

    KDL::Frame rhs(k - p);
    ScrewTheoryIkSubproblem::Solutions actual;
    ASSERT_TRUE(pg4.solve(rhs, KDL::Frame::Identity(), actual));

    ASSERT_EQ(actual.size(), 2);
    ASSERT_EQ(actual[0].size(), 2);
    ASSERT_EQ(actual[1].size(), 2);

    ScrewTheoryIkSubproblem::Solutions expected = {
        {KDL::PI_2, KDL::PI_2},
        {KDL::PI, -KDL::PI_2}
    };

    checkSolutions(actual, expected);

    KDL::Vector k2 = k + KDL::Vector(0, 1, 0);
    KDL::Frame rhs2(k2 - p);
    ASSERT_FALSE(pg4.solve(rhs2, KDL::Frame::Identity(), actual));

    checkSolutions(actual, expected);

    KDL::Vector p3 = p + KDL::Vector(0.75, 0, 0);
    KDL::Vector k3 = k + KDL::Vector(-0.75, 0, -0.75);
    PardosGotorFour pg4c(exp1, exp2, p3);
    KDL::Frame rhs3(k3 - p3);
    ASSERT_FALSE(pg4c.solve(rhs3, KDL::Frame::Identity(), actual));

    expected = {
        {3 * KDL::PI_4, KDL::PI},
        {3 * KDL::PI_4, KDL::PI}
    };

    checkSolutions(actual, expected);

    KDL::Vector k4 = k + KDL::Vector(-1, 0, -1);
    KDL::Frame rhs4(k4 - p);
    ASSERT_TRUE(pg4.solve(rhs4, KDL::Frame::Identity(), actual));

    expected = {
        {0, KDL::PI},
        {0, KDL::PI}
    };

    checkSolutions(actual, expected);

    expected = {
        {KDL::PI, KDL::PI},
        {KDL::PI, KDL::PI}
    };

    ASSERT_TRUE(pg4.solve(rhs4, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Vector k5 = k4 + KDL::Vector(0, 1, 0);
    KDL::Frame rhs5(k5 - p);
    ASSERT_FALSE(pg4.solve(rhs5, KDL::Frame::Identity(), actual));

    expected = {
        {0, KDL::PI},
        {0, KDL::PI}
    };

    checkSolutions(actual, expected);

    expected = {
        {KDL::PI, KDL::PI},
        {KDL::PI, KDL::PI}
    };

    ASSERT_FALSE(pg4.solve(rhs5, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Vector p6 = p + KDL::Vector(1, 0, 0);
    KDL::Vector k6 = k + KDL::Vector(-1, 0, 0);
    PardosGotorFour pg4d(exp1, exp2, p6);
    KDL::Frame rhs6(k6 - p6);
    ASSERT_TRUE(pg4d.solve(rhs6, KDL::Frame::Identity(), actual));

    expected = {
        {KDL::PI_2, 0},
        {KDL::PI_2, 0}
    };

    checkSolutions(actual, expected);

    expected = {
        {KDL::PI_2, KDL::PI},
        {KDL::PI_2, KDL::PI}
    };

    ASSERT_TRUE(pg4d.solve(rhs6, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Frame rhs7(k - p6);
    ASSERT_FALSE(pg4d.solve(rhs7, KDL::Frame::Identity(), actual));

    expected = {
        {3 * KDL::PI_4, 0},
        {3 * KDL::PI_4, 0}
    };

    checkSolutions(actual, expected);

    expected = {
        {3 * KDL::PI_4, KDL::PI},
        {3 * KDL::PI_4, KDL::PI}
    };

    ASSERT_FALSE(pg4d.solve(rhs7, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Frame rhs8(k4 - p6);
    ASSERT_FALSE(pg4d.solve(rhs8, KDL::Frame::Identity(), actual));

    expected = {
        {0, 0},
        {0, 0}
    };

    checkSolutions(actual, expected);

    expected = {
        {KDL::PI, KDL::PI},
        {KDL::PI, KDL::PI}
    };

    ASSERT_FALSE(pg4d.solve(rhs8, KDL::Frame::Identity(), expected[0], actual));

    checkSolutions(actual, expected);

    KDL::Vector p7(0, 1, 0); 
    KDL::Vector k7(2, 1, 2); 
    PardosGotorFour pg4e(exp1, exp2, p7);
    KDL::Frame rhs9(k7 - p7);

    ASSERT_TRUE(pg4e.solve(rhs9, KDL::Frame::Identity(), actual));

    expected = {
        {KDL::PI_2, 0},
        {KDL::PI_2, 0}
    };

    checkSolutions(actual, expected);

}

TEST_F(ScrewTheoryTest, PardosGotorFive)
{
    KDL::Vector p(0, 1, 1);
    KDL::Vector k(2, 0, 0);

    MatrixExponential exp(MatrixExponential::ROTATION, {0, 0, 1}, {0, 0, 0});
    MatrixExponential exp2(MatrixExponential::ROTATION, {1, 0, 0}, {0, 0, 0});
    PardosGotorFive pg5(exp, exp2, p);

    ASSERT_EQ(pg5.solutions(), 2);  

    KDL::Frame rhs(k - p);
    ScrewTheoryIkSubproblem::Solutions actual;
    ASSERT_TRUE(pg5.solve(rhs, KDL::Frame::Identity(), actual));

    ASSERT_EQ(actual.size(), 2);
    ASSERT_EQ(actual[0].size(), 1);
    ASSERT_EQ(actual[1].size(), 1);

    ScrewTheoryIkSubproblem::Solutions expected = {
        {-KDL::PI_2},
        {KDL::PI_2},
    };

    checkSolutions(actual, expected);

    //caso del ajuste

    KDL::Vector p2(0.3, 1, 1);

    PardosGotorFive pg5b(exp, exp2, p2);

    KDL::Frame rhs2(k - p2);
    ASSERT_TRUE(pg5b.solve(rhs2, KDL::Frame::Identity(), actual));

    expected = {
        {-1.138451},
        {2.304278},
    };

    checkSolutions(actual, expected);
}

TEST_F(ScrewTheoryTest, PardosGotorSix)
{
    //Dos ejes que se cruzan normales
    KDL::Vector p(-1, 0, 0);
    KDL::Vector k(2, 1, 2);

    MatrixExponential exp1(MatrixExponential::ROTATION, KDL::Vector(0, 1, 0), KDL::Vector(2, 0, 0));
    MatrixExponential exp2(MatrixExponential::ROTATION, KDL::Vector(0, 0, 1), KDL::Vector(0, 0, 0));
    PardosGotorSix pg6(exp1, exp2, p);

    ASSERT_EQ(pg6.solutions(), 1);

    KDL::Frame rhs(k - p);
    ScrewTheoryIkSubproblem::Solutions actual;
    ASSERT_TRUE(pg6.solve(rhs, KDL::Frame::Identity(), actual));

    ASSERT_EQ(actual.size(), 1);
    ASSERT_EQ(actual[0].size(), 2);

    ScrewTheoryIkSubproblem::Solutions expected = {{KDL::PI_2, -KDL::PI_2}};

    checkSolutions(actual, expected);

    //Desplaza el punto k del caso anterior a una posición que no es alcanzable
    KDL::Vector k2(2, 2, 2);
    KDL::Frame rhs2(k2 - p);
    ASSERT_FALSE(pg6.solve(rhs2, KDL::Frame::Identity(), actual));

    checkSolutions(actual, expected);

    //Desplaza el punto k del caso anterior a una posición que no es alcanzable, aunque los resultados aproximados cambiarán
    KDL::Vector k3(2, 0, 0.25);
    KDL::Frame rhs3(k3 - p);
    ASSERT_FALSE(pg6.solve(rhs3, KDL::Frame::Identity(), actual));

    expected = {{KDL::PI_2, KDL::PI}};

    checkSolutions(actual, expected);

    //similar al caso anterior
    KDL::Vector p4(-1, 0, 2);
    KDL::Vector k4(2, 1, -2);
    KDL::Frame rhs4(k4 - p4);
    PardosGotorSix pg6d(exp1, exp2, p4);
    ASSERT_FALSE(pg6d.solve(rhs4, KDL::Frame::Identity(), actual));

    expected = {{-3 * KDL::PI_4, -KDL::PI_2}};

    checkSolutions(actual, expected);

    //si el punto p se encuentra en el eje 2
    KDL::Vector p5(0, 0, 0);
    KDL::Vector k5(2, 0, -2);
    KDL::Frame rhs5(k5 - p5);
    PardosGotorSix pg6e(exp1, exp2, p5);
    ASSERT_TRUE(pg6e.solve(rhs5, KDL::Frame::Identity(), actual));

    expected = {{-KDL::PI_2, 0}};

    checkSolutions(actual, expected);

    //si el punto p se encuentra en el eje 2, pero no es alcanzable
    KDL::Vector p6(0, 0, 0);
    KDL::Vector k6(2, 0, -1.5);
    KDL::Frame rhs6(k6 - p6);
    PardosGotorSix pg6f(exp1, exp2, p6);
    ASSERT_FALSE(pg6f.solve(rhs6, KDL::Frame::Identity(), actual));

    expected = {{-KDL::PI_2, 0}};

    checkSolutions(actual, expected);

    //si el punto p se encuentra en el eje 2, no es alcanzable y el angulo cambia
    KDL::Vector p7(0, 0, 0);
    KDL::Vector k7(5, 0, 0);
    KDL::Frame rhs7(k7 - p7);
    PardosGotorSix pg6g(exp1, exp2, p7);
    ASSERT_FALSE(pg6g.solve(rhs7, KDL::Frame::Identity(), actual));

    expected = {{KDL::PI, 0}};

    checkSolutions(actual, expected);

    //si el punto p se encuentra en el eje 1
    KDL::Vector p8(0, 2, 0);
    KDL::Vector k8(2, 0, 0);
    KDL::Frame rhs8(k8 - p8);
    PardosGotorSix pg6h(exp1, exp2, p8);
    ASSERT_TRUE(pg6h.solve(rhs8, KDL::Frame::Identity(), actual));

    expected = {{0, -KDL::PI_2}};

    checkSolutions(actual, expected);

    //si el punto p se encuentra en el eje 2, pero no es alcanzable
    KDL::Vector p9(0, 2, 0);
    KDL::Vector k9(3, 0, 0);
    KDL::Frame rhs9(k9 - p9);
    PardosGotorSix pg6i(exp1, exp2, p9);
    ASSERT_FALSE(pg6i.solve(rhs9, KDL::Frame::Identity(), actual));

    expected = {{0, -KDL::PI_2}};

    checkSolutions(actual, expected);

    //si el punto p se encuentra en el eje 2, no es alcanzable y el ángulo cambia
    KDL::Vector p10(-1, 0, 0);
    KDL::Vector k10(2, 0, 0);
    KDL::Frame rhs10(k10 - p10);
    PardosGotorSix pg6j(exp1, exp2, p10);
    ASSERT_FALSE(pg6j.solve(rhs10, KDL::Frame::Identity(), actual));

    expected = {{0, KDL::PI}};

    checkSolutions(actual, expected);

    //si p y k están en el mismo punto, que pertenece al eje 2
    KDL::Vector p11(0, 0, 0);
    KDL::Vector k11(0, 0, 0);
    KDL::Frame rhs11(k11 - p11);
    PardosGotorSix pg6k(exp1, exp2, p11);
    ASSERT_TRUE(pg6k.solve(rhs11, KDL::Frame::Identity(), actual));

    expected = {{0, 0}};

    checkSolutions(actual, expected);

    //si p y k están en el mismo punto, que pertenece al eje 1
    KDL::Vector p12(2, 0, 0);
    KDL::Vector k12(2, 0, 0);
    KDL::Frame rhs12(k12 - p12);
    PardosGotorSix pg6l(exp1, exp2, p12);
    ASSERT_TRUE(pg6l.solve(rhs12, KDL::Frame::Identity(), actual));

    expected = {{0, 0}};

    checkSolutions(actual, expected);

    //si p y k están en el mismo punto, que no pertenece a ningún eje. DABA ERROR, he añadido una condición en el else de pg6.solve
    //para que no devuelva falso si los puntos p y k son iguales
    KDL::Vector p13(1, 0, 0);
    KDL::Vector k13(1, 0, 0);
    KDL::Frame rhs13(k13 - p13);
    PardosGotorSix pg6m(exp1, exp2, p13);
    ASSERT_TRUE(pg6m.solve(rhs13, KDL::Frame::Identity(), actual));

    expected = {{0, 0}};

    checkSolutions(actual, expected);

    //si p y k se encuentran en los ejes 1 y 2
    KDL::Vector p14(0, 0, 0);
    KDL::Vector k14(2, 0, 0);
    KDL::Frame rhs14(k14 - p14);
    PardosGotorSix pg6n(exp1, exp2, p14);
    ASSERT_FALSE(pg6n.solve(rhs14, KDL::Frame::Identity(), actual));

    expected = {{0, 0}};

    checkSolutions(actual, expected);
}

TEST_F(ScrewTheoryTest, PardosGotorSeven)
{
    //dos paralelas y una que se cruza, caso con una única solución, ya que los puntos intermedios de los ejes paralelos son iguales
    KDL::Vector p(-3, 0, 0);
    KDL::Vector k(2, 1, 2);

    MatrixExponential exp1(MatrixExponential::ROTATION, {0, 1, 0}, {2, 0, 0});
    MatrixExponential exp2(MatrixExponential::ROTATION, {0, 0, 1}, {0, 0, 0});
    MatrixExponential exp3(MatrixExponential::ROTATION, {0, 0, 1}, {-2, 0, 0});

    PardosGotorSeven pg7(exp1, exp2, exp3, p);
    ASSERT_EQ(pg7.solutions(), 4);

    KDL::Frame rhs(k - p);
    ScrewTheoryIkSubproblem::Solutions actual;
    ASSERT_TRUE(pg7.solve(rhs, KDL::Frame::Identity(), actual));

    ASSERT_EQ(actual.size(), 4);
    ASSERT_EQ(actual[0].size(), 3);
    ASSERT_EQ(actual[1].size(), 3);
    ASSERT_EQ(actual[2].size(), 3);
    ASSERT_EQ(actual[3].size(), 3);

    ScrewTheoryIkSubproblem::Solutions expected = {
        {KDL::PI_2, -KDL::PI_2, KDL::PI},
        {KDL::PI_2, -KDL::PI_2, KDL::PI},
        {KDL::PI_2, -KDL::PI_2, KDL::PI},
        {KDL::PI_2, -KDL::PI_2, KDL::PI}
    };

    checkSolutions(actual, expected);

    //caso anterior pero con posición no alcanzable
    KDL::Vector p2(-3, 0, 0);
    KDL::Vector k2(2, 0, 2);

    PardosGotorSeven pg7b(exp1, exp2, exp3, p2);

    KDL::Frame rhs2(k2 - p2);

    ASSERT_FALSE(pg7b.solve(rhs2, KDL::Frame::Identity(), actual));

    expected = {
        {KDL::PI_2, -KDL::PI_2, KDL::PI},
        {KDL::PI_2, -KDL::PI_2, KDL::PI},
        {KDL::PI_2, -KDL::PI_2, KDL::PI},
        {KDL::PI_2, -KDL::PI_2, KDL::PI}
    };

    checkSolutions(actual, expected);

    //dos paralelas y una que se cruza, caso con dos soluciones, ya que los puntos intermedios de los ejes paralelos son diferentes
    KDL::Vector p3(0, 1, 0);
    KDL::Vector k3(4, 0, 1);

    MatrixExponential exp1_b(MatrixExponential::ROTATION, {0, 0, 1}, {3, 0, 0});
    MatrixExponential exp2_b(MatrixExponential::ROTATION, {0, 1, 0}, {2, 0, 0});
    MatrixExponential exp3_b(MatrixExponential::ROTATION, {0, 1, 0}, {1, 0, 0});
    PardosGotorSeven pg7c(exp1_b, exp2_b, exp3_b, p3);

    KDL::Frame rhs3(k3 - p3);

    ASSERT_TRUE(pg7c.solve(rhs3, KDL::Frame::Identity(), actual));

    expected = {
        {-KDL::PI_2, KDL::PI_2, KDL::PI_2},
        {-KDL::PI_2, KDL::PI, -KDL::PI_2},
        {-KDL::PI_2, KDL::PI_2, KDL::PI_2},
        {-KDL::PI_2, KDL::PI, -KDL::PI_2}
    };

    checkSolutions(actual, expected);

    //caso anterior con posición no alcanzable
    KDL::Vector p4(0, 1, 0);
    KDL::Vector k4(4, 0, 2);

    PardosGotorSeven pg7d(exp1_b, exp2_b, exp3_b, p3);

    KDL::Frame rhs4(k4 - p4);

    ASSERT_FALSE(pg7d.solve(rhs4, KDL::Frame::Identity(), actual));

    expected = {
        {-KDL::PI_2, KDL::PI_2, KDL::PI_2},
        {-KDL::PI_2, KDL::PI, -KDL::PI_2},
        {-KDL::PI_2, KDL::PI_2, KDL::PI_2},
        {-KDL::PI_2, KDL::PI, -KDL::PI_2}
    };

    checkSolutions(actual, expected);

    //dos paralelas y una que se corta, caso con cuatro soluciones, ya que los puntos intermedios de los ejes paralelos son diferentes
    KDL::Vector p5(0, 1, 0);
    KDL::Vector k5(2, 0, 1);

    MatrixExponential exp1_c(MatrixExponential::ROTATION, {0, 0, 1}, {2, 1, 0});
    MatrixExponential exp2_c(MatrixExponential::ROTATION, {0, 1, 0}, {2, 0, 0});
    MatrixExponential exp3_c(MatrixExponential::ROTATION, {0, 1, 0}, {1, 0, 0});
    PardosGotorSeven pg7e(exp1_c, exp2_c, exp3_c, p5);

    KDL::Frame rhs5(k5 - p5);

    ASSERT_TRUE(pg7e.solve(rhs5, KDL::Frame::Identity(), actual));

    expected = {
        {KDL::PI_2, 0, KDL::PI_2},
        {KDL::PI_2, KDL::PI_2, -KDL::PI_2},
        {-KDL::PI_2, KDL::PI_2, KDL::PI_2},
        {-KDL::PI_2, KDL::PI, -KDL::PI_2}
    };

    checkSolutions(actual, expected);

    //caso anterior con posición no alcanzable    
    KDL::Vector p6(0, 1, 0);
    KDL::Vector k6(2, 0, 2);

    PardosGotorSeven pg7f(exp1_c, exp2_c, exp3_c, p6);

    KDL::Frame rhs6(k6 - p6);

    ASSERT_FALSE(pg7f.solve(rhs6, KDL::Frame::Identity(), actual));

    expected = {
        {KDL::PI_2, 0, KDL::PI_2},
        {KDL::PI_2, KDL::PI_2, -KDL::PI_2},
        {-KDL::PI_2, KDL::PI_2, KDL::PI_2},
        {-KDL::PI_2, KDL::PI, -KDL::PI_2}
    };

    checkSolutions(actual, expected);

    //dos paralelas y una que se corta, caso con dos soluciones, ya que los puntos intermedios de los ejes paralelos son iguales
    KDL::Vector p7(-3, 0, 0);
    KDL::Vector k7(0, 1, 0);

    MatrixExponential exp1_d(MatrixExponential::ROTATION, {0, 1, 0}, {0, 0, 0});
    MatrixExponential exp2_d(MatrixExponential::ROTATION, {0, 0, 1}, {0, 0, 0});
    MatrixExponential exp3_d(MatrixExponential::ROTATION, {0, 0, 1}, {-2, 0, 0});
    PardosGotorSeven pg7g(exp1_d, exp2_d, exp3_d, p7);

    KDL::Frame rhs7(k7 - p7);

    ASSERT_TRUE(pg7g.solve(rhs7, KDL::Frame::Identity(), actual));

    expected = {
        {0, -KDL::PI_2, KDL::PI},
        {0, -KDL::PI_2, KDL::PI},
        {0, -KDL::PI_2, KDL::PI},
        {0, -KDL::PI_2, KDL::PI}
    };

    checkSolutions(actual, expected);

    //caso anterior con posición no alcanzable  
    KDL::Vector p8(-3, 0, 0);
    KDL::Vector k8(0, 0, 0);

    PardosGotorSeven pg7h(exp1_d, exp2_d, exp3_d, p8);

    KDL::Frame rhs8(k8 - p8);

    ASSERT_FALSE(pg7h.solve(rhs8, KDL::Frame::Identity(), actual));

    expected = {
        {0, -KDL::PI_2, KDL::PI},
        {0, -KDL::PI_2, KDL::PI},
        {0, -KDL::PI_2, KDL::PI},
        {0, -KDL::PI_2, KDL::PI}
    };

    checkSolutions(actual, expected);

    //p pertenece al eje 3  
    KDL::Vector p9(-2, 0, 0);
    KDL::Vector k9(2, 2, 2);

    PardosGotorSeven pg7i(exp1, exp2, exp3, p9);

    KDL::Frame rhs9(k9 - p9);

    ASSERT_TRUE(pg7i.solve(rhs9, KDL::Frame::Identity(), actual));

    expected = {
        {KDL::PI_2, -KDL::PI_2, 0},
        {KDL::PI_2, -KDL::PI_2, 0},
        {KDL::PI_2, -KDL::PI_2, 0},
        {KDL::PI_2, -KDL::PI_2, 0}
    };

    checkSolutions(actual, expected);

    //k pertenece al punto final de la rotación entre los paralelos
    KDL::Vector p10(-3, 0, 0);
    KDL::Vector k10(0, 1, 0);

    PardosGotorSeven pg7j(exp1, exp2, exp3, p10);

    KDL::Frame rhs10(k10 - p10);

    ASSERT_TRUE(pg7j.solve(rhs10, KDL::Frame::Identity(), actual));

    expected = {
        {0, -KDL::PI_2, KDL::PI},
        {0, -KDL::PI_2, KDL::PI},
        {0, -KDL::PI_2, KDL::PI},
        {0, -KDL::PI_2, KDL::PI}
    };

    checkSolutions(actual, expected);

    //k pertenece al eje 1 
    KDL::Vector p11(0, 1, 0);
    KDL::Vector k11(3, 1, 1);

    MatrixExponential exp1_e(MatrixExponential::ROTATION, {0, 0, 1}, {3, 1, 0});
    MatrixExponential exp2_e(MatrixExponential::ROTATION, {0, 1, 0}, {3, 0, 0});
    MatrixExponential exp3_e(MatrixExponential::ROTATION, {0, 1, 0}, {1, 0, 0});
    PardosGotorSeven pg7k(exp1_e, exp2_e, exp3_e, p11);

    KDL::Frame rhs11(k11 - p11);

    ASSERT_TRUE(pg7k.solve(rhs11, KDL::Frame::Identity(), actual));

    expected = {
        {0, KDL::PI_2, KDL::PI},
        {0, KDL::PI_2, KDL::PI},
        {0, KDL::PI_2, KDL::PI},
        {0, KDL::PI_2, KDL::PI}
    };

    checkSolutions(actual, expected);

    //p y k son iguales
    KDL::Vector p12(-3, 0, 0);
    KDL::Vector k12(-3, 0, 0);

    PardosGotorSeven pg7l(exp1, exp2, exp3, p12);

    KDL::Frame rhs12(k12 - p12);

    ASSERT_TRUE(pg7l.solve(rhs12, KDL::Frame::Identity(), actual));

    expected = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
    };

    checkSolutions(actual, expected);
}

TEST_F(ScrewTheoryTest, PardosGotorEight)
{
    KDL::Vector p(-1, 1, 0);
    KDL::Vector k(4, 1, 1);

    PoeExpression poe = makeUR16eFromPoE();
    KDL::Chain chain = makeUR16eFromDh();

    KDL::ChainFkSolverPos_recursive fkSolver(chain);
    KDL::Frame H_S_T_q_DH, H_S_T_q_ST;
    KDL::JntArray q = fillJointValues(6, 0);
    ASSERT_TRUE(poe.evaluate(q, H_S_T_q_ST));

    MatrixExponential exp1(MatrixExponential::ROTATION, {0, 1, 0}, {0, 0, 0.181});
    MatrixExponential exp2(MatrixExponential::ROTATION, {0, 1, 0}, {0.478, 0, 0.181});
    MatrixExponential exp3(MatrixExponential::ROTATION, {0, 1, 0}, {0.838, 0.174, 0.181}); //la coordenada x era -1
    PardosGotorEight pg8(exp1, exp2, exp3, p, 1, 3, poe);

    ASSERT_EQ(pg8.solutions(), 2);

    KDL::JntArray sol(6);  // 6 es el número de articulaciones
    sol(0) = 0; sol(1) =  0; sol(2) =  0; sol(3) =  0; sol(4) = 0; sol(5) =  0;

    ScrewTheoryIkSubproblem::Solutions expected = {
        {0, 0, 0},
        {0, 0, 0}
    };

    KDL::Frame rhs(k - p);
    ScrewTheoryIkSubproblem::Solutions actual;
    ASSERT_TRUE(pg8.solve(rhs, KDL::Frame::Identity(), expected[0], actual, H_S_T_q_ST, sol, poe.getTransform()));

    ASSERT_EQ(actual.size(), 2);
    ASSERT_EQ(actual[0].size(), 3);
    ASSERT_EQ(actual[1].size(), 3);

    checkSolutions(actual, expected);

    q = fillJointValues(6, KDL::PI_2);
    ASSERT_TRUE(poe.evaluate(q, H_S_T_q_ST));
    PardosGotorEight pg8b(exp1, exp2, exp3, p, 1, 3, poe);

    sol(0) = -2.825412711346422; sol(1) =  0; sol(2) =  0; sol(3) =  0; sol(4) = -1.570796326794897; sol(5) =  2.825412711346421;

    expected = {
        {0.279783, 1.570796, -0.279783},
        {1.570796, -1.570796, 1.570796}
    };

    ASSERT_TRUE(pg8.solve(rhs, KDL::Frame::Identity(), expected[0], actual, H_S_T_q_ST, sol, poe.getTransform()));

    checkSolutions(actual, expected);
}

TEST_F(ScrewTheoryTest, AbbIrb120Kinematics)
{
    KDL::Chain chain = makeAbbIrb120KinematicsFromDH();
    PoeExpression poe = makeAbbIrb120KinematicsFromPoE();

    checkRobotKinematics(chain, poe, 8);
}

TEST_F(ScrewTheoryTest, PumaKinematics)
{
    KDL::Chain chain = makePumaKinematicsFromDH();
    PoeExpression poe = makePumaKinematicsFromPoE();

    checkRobotKinematics(chain, poe, 8);
}

TEST_F(ScrewTheoryTest, StanfordKinematics)
{
    KDL::Chain chain = makeStanfordKinematicsFromDH();
    PoeExpression poe = makeStanfordKinematicsFromPoE();

    checkRobotKinematics(chain, poe, 8);
}

TEST_F(ScrewTheoryTest, AbbIrb910scKinematics)
{
    KDL::Chain chain = makeAbbIrb910scKinematicsFromDH();
    PoeExpression poe = makeAbbIrb910scKinematicsFromPoE();

    checkRobotKinematics(chain, poe, 2);
}

TEST_F(ScrewTheoryTest, AbbIrb6620lxKinematics)
{
    KDL::Chain chain = makeAbbIrb6620lxFromDh();
    PoeExpression poe = makeAbbIrb6620lxFromPoE();

    checkRobotKinematics(chain, poe, 4);
}

TEST_F(ScrewTheoryTest, UR16eKinematics)
{
    KDL::Chain chain = makeUR16eFromDh();
    PoeExpression poe = makeUR16eFromPoE();

    checkRobotKinematics(chain, poe, 8);
}

TEST_F(ScrewTheoryTest, TeoRightArmKinematics)
{
    KDL::Chain chain = makeTeoRightArmKinematicsFromDH();
    PoeExpression poe = makeTeoRightArmKinematicsFromPoE();

    checkRobotKinematics(chain, poe, 8);
}

TEST_F(ScrewTheoryTest, TeoRightLegKinematics)
{
    KDL::Chain chain = makeTeoRightLegKinematicsFromDH();
    PoeExpression poe = makeTeoRightLegKinematicsFromPoE();

    checkRobotKinematics(chain, poe, 8);
}

TEST_F(ScrewTheoryTest, ConfigurationSelector)
{
    PoeExpression poe = makeTeoRightArmKinematicsFromPoE();

    KDL::JntArray q(poe.size());
    q(3) = KDL::PI_2; // elbow

    KDL::Frame H;
    ASSERT_TRUE(poe.evaluate(q, H));

    ScrewTheoryIkProblemBuilder builder(poe);
    ScrewTheoryIkProblem * ikProblem = builder.build();

    ASSERT_TRUE(ikProblem);
    ASSERT_EQ(ikProblem->solutions(), 8);

    ScrewTheoryIkProblem::Solutions solutions;
    auto reachability = ikProblem->solve(H, q, solutions);
    ASSERT_TRUE(std::all_of(reachability.begin(), reachability.end(), [](bool b) { return b; }));

    KDL::JntArray qMin = fillJointValues(poe.size(), -KDL::PI);
    KDL::JntArray qMax = fillJointValues(poe.size(), KDL::PI);

    ConfigurationSelectorLeastOverallAngularDisplacementFactory confFactory(qMin, qMax);
    ConfigurationSelector * config = confFactory.create();

    ASSERT_TRUE(config);
    ASSERT_TRUE(config->configure(solutions, reachability));
    ASSERT_TRUE(config->findOptimalConfiguration(q));

    KDL::JntArray qSolved;
    config->retrievePose(qSolved);
    int n1 = findTargetConfiguration(solutions, q);

    ASSERT_NE(n1, -1);

    H.p += {0.01, 0, 0}; // add a tiny displacement

    reachability = ikProblem->solve(H, q, solutions);
    ASSERT_TRUE(std::all_of(reachability.begin(), reachability.end(), [](bool b) { return b; }));
    delete ikProblem;

    ASSERT_TRUE(config->configure(solutions, reachability));
    ASSERT_TRUE(config->findOptimalConfiguration(q));

    config->retrievePose(qSolved);
    int n2 = findTargetConfiguration(solutions, qSolved);

    ASSERT_EQ(n2, n1);
    delete config;
}

TEST_F(ScrewTheoryTest, ConfigurationSelectorGait)
{
    PoeExpression poe = makeTeoRightLegKinematicsFromPoE();

    KDL::JntArray q(poe.size());
    q(2) = -0.3; // approx. 20 degrees
    q(3) = 0.6;
    q(4) = -0.3;

    KDL::Frame H;
    ASSERT_TRUE(poe.evaluate(q, H));

    ScrewTheoryIkProblemBuilder builder(poe);
    ScrewTheoryIkProblem * ikProblem = builder.build();

    ASSERT_TRUE(ikProblem);
    ASSERT_EQ(ikProblem->solutions(), 8);

    ScrewTheoryIkProblem::Solutions solutions;
    auto reachability = ikProblem->solve(H, q, solutions);
    ASSERT_TRUE(std::all_of(reachability.begin(), reachability.end(), [](bool b) { return b; }));

    KDL::JntArray qMin = fillJointValues(poe.size(), -KDL::PI);
    KDL::JntArray qMax = fillJointValues(poe.size(), KDL::PI);

    ConfigurationSelectorHumanoidGaitFactory confFactory(qMin, qMax);
    ConfigurationSelector * config = confFactory.create();

    KDL::JntArray qInitial(poe.size());

    ASSERT_TRUE(config);
    ASSERT_TRUE(config->configure(solutions, reachability));
    ASSERT_TRUE(config->findOptimalConfiguration(qInitial));

    KDL::JntArray qSolved;
    config->retrievePose(qSolved);
    int n1 = findTargetConfiguration(solutions, q);

    ASSERT_NE(n1, -1);
    delete config;
}

} // namespace roboticslab::test
