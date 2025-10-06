// -*- mode:C++; tab-width:4; c-basic-offset:4; indent-tabs-mode:nil -*-

#ifndef __SCREW_THEORY_IK_SUBPROBLEMS_HPP__
#define __SCREW_THEORY_IK_SUBPROBLEMS_HPP__

#include <kdl/frames.hpp>

#include "ScrewTheoryIkProblem.hpp"
#include "MatrixExponential.hpp"

namespace roboticslab
{

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief First Paden-Kahan subproblem
 *
 * Single solution, single revolute joint geometric IK subproblem given by
 * @f$ e\,^{\hat{\xi}\,{\theta}} \cdot p = k @f$
 * (rotation screw applied to a point).
 */
class PadenKahanOne : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param exp POE term.
     * @param p Characteristic point.
     */
    PadenKahanOne(const MatrixExponential & exp, const KDL::Vector & p);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override;

    int solutions() const override
    { return 1; }

    const char * describe() const override
    { return "PK1"; }

private:
    const MatrixExponential exp;
    const KDL::Vector p;
    const KDL::Rotation axisPow;
};

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief Second Paden-Kahan subproblem
 *
 * Dual solution, double revolute joint geometric IK subproblem given by
 * @f$ e\,^{\hat{\xi_1}\,{\theta_1}} \cdot e\,^{\hat{\xi_2}\,{\theta_2}} \cdot p = k @f$
 * (consecutive crossing rotation screws to a point).
 */
class PadenKahanTwo : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param exp1 First POE term.
     * @param exp2 Second POE term.
     * @param p Characteristic point.
     * @param r Point of intersection between both screw axes.
     */
    PadenKahanTwo(const MatrixExponential & exp1, const MatrixExponential & exp2, const KDL::Vector & p, const KDL::Vector & r);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override;

    int solutions() const override
    { return 2; }

    const char * describe() const override
    { return "PK2"; }

private:
    const MatrixExponential exp1, exp2;
    const KDL::Vector p, r, axesCross;
    const KDL::Rotation axisPow1, axisPow2;
    const double axesDot;
};

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief Third Paden-Kahan subproblem
 *
 * Dual solution, single revolute joint geometric IK subproblem given by
 * @f$ \left \| e\,^{\hat{\xi}\,{\theta}} \cdot p - k \right \| = \delta @f$
 * (rotation screw for moving @f$ p @f$ to a distance @f$ \delta @f$ from @f$ k @f$).
 */
class PadenKahanThree : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param exp POE term.
     * @param p First characteristic point.
     * @param k Second characteristic point.
     */
    PadenKahanThree(const MatrixExponential & exp, const KDL::Vector & p, const KDL::Vector & k);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override;

    int solutions() const override
    { return 2; }

    const char * describe() const override
    { return "PK3"; }

private:
    const MatrixExponential exp;
    const KDL::Vector p, k;
    const KDL::Rotation axisPow;
};

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief First Pardos-Gotor subproblem
 *
 * Single solution, single prismatic joint geometric IK subproblem given by
 * @f$ e\,^{\hat{\xi}\,{\theta}} \cdot p = k @f$
 * (translation screw applied to a point, see @cite pardosgotor2018str
 * @cite pardosgotor2022str).
 */
class PardosGotorOne : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param exp POE term.
     * @param p Characteristic point.
     */
    PardosGotorOne(const MatrixExponential & exp, const KDL::Vector & p);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override;

    int solutions() const override
    { return 1; }

    const char * describe() const override
    { return "PG1"; }

private:
    const MatrixExponential exp;
    const KDL::Vector p;
};

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief Second Pardos-Gotor subproblem
 *
 * Single solution, double prismatic joint geometric IK subproblem given by
 * @f$ e\,^{\hat{\xi_1}\,{\theta_1}} \cdot e\,^{\hat{\xi_2}\,{\theta_2}} \cdot p = k @f$
 * (consecutive translation screws to a point, see @cite pardosgotor2018str
 * @cite pardosgotor2022str).
 */
class PardosGotorTwo : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param exp1 First POE term.
     * @param exp2 Second POE term.
     * @param p Characteristic point.
     */
    PardosGotorTwo(const MatrixExponential & exp1, const MatrixExponential & exp2, const KDL::Vector & p);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override;

    int solutions() const override
    { return 1; }

    const char * describe() const override
    { return "PG2"; }

private:
    const MatrixExponential exp1, exp2;
    const KDL::Vector p, crossPr2;
    const double crossPr2Norm;
};

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief Third Pardos-Gotor subproblem
 *
 * Dual solution, single prismatic joint geometric IK subproblem given by
 * @f$ \left \| e\,^{\hat{\xi}\,{\theta}} \cdot p - k \right \| = \delta @f$
 * (translation screw for moving @f$ p @f$ to a distance @f$ \delta @f$ from @f$ k @f$,
 * see @cite pardosgotor2018str @cite pardosgotor2022str).
 */
class PardosGotorThree : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param exp POE term.
     * @param p First characteristic point.
     * @param k Second characteristic point.
     */
    PardosGotorThree(const MatrixExponential & exp, const KDL::Vector & p, const KDL::Vector & k);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override;

    int solutions() const override
    { return 2; }

    const char * describe() const override
    { return "PG3"; }

private:
    const MatrixExponential exp;
    const KDL::Vector p, k;
};

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief Fourth Pardos-Gotor subproblem
 *
 * Dual solution, double revolute joint geometric IK subproblem given by
 * @f$ e\,^{\hat{\xi_1}\,{\theta_1}} \cdot e\,^{\hat{\xi_2}\,{\theta_2}} \cdot p = k @f$
 * (consecutive parallel rotation screws applied to a point,
 * see @cite pardosgotor2018str @cite pardosgotor2022str).
 */
class PardosGotorFour : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param exp1 First POE term.
     * @param exp2 Second POE term.
     * @param p Characteristic point.
     */
    PardosGotorFour(const MatrixExponential & exp1, const MatrixExponential & exp2, const KDL::Vector & p);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override;

    int solutions() const override
    { return 2; }

    const char * describe() const override
    { return "PG4"; }

private:
    const MatrixExponential exp1, exp2;
    const KDL::Vector p, n;
    const KDL::Rotation axisPow;
};

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief Fifth Pardos-Gotor subproblem
 * 
 * Dual solution, single revolute joint geometric IK subproblem given by
 * @f$ e\,^{\hat{\xi}\,{\theta}} \cdot \omega_p = \omega_k \ ; e\,^{\hat{\xi}\,{\theta}} \cdot \Pi_p = \Pi_k @f$
 * (rotation around one single axis applied to a perpendicular line or plane,
 * see @cite pardosgotor2022str).
 */
class PardosGotorFive : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param exp POE term.
     * @param exp_next POE term.
     * @param p Characteristic point.
     */
    PardosGotorFive(const MatrixExponential & exp, const MatrixExponential & exp_next, const KDL::Vector & p);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override;

    int solutions() const override
    { return 2; }

    const char * describe() const override
    { return "PG5"; }

private:
    const MatrixExponential exp;
    const MatrixExponential exp_next;
    const KDL::Vector p;
    const KDL::Rotation axisPow;
};

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief Sixth Pardos-Gotor subproblem
 *
 * Single solution, double revolute joint geometric IK subproblem given by
 * @f$ e\,^{\hat{\xi_1}\,{\theta_1}} \cdot e\,^{\hat{\xi_2}\,{\theta_2}} \cdot p = k @f$
 * (consecutive skew rotation screws applied to a point,
 * see @cite pardosgotor2022str).
 */
class PardosGotorSix : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;
    /**
     * @brief Constructor
     *
     * @param exp1 First POE term.
     * @param exp2 Second POE term.
     * @param p Characteristic point.
     */
    PardosGotorSix(const MatrixExponential & exp1, const MatrixExponential & exp2, const KDL::Vector & p);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override;

    int solutions() const override
    { return 1; }

    const char * describe() const override
    { return "PG6"; }

private:
    const MatrixExponential exp1, exp2;
    const KDL::Vector p, axesCross;
    const KDL::Rotation axisPow1, axisPow2;
    const double axesDot;
};

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief Seventh Pardos-Gotor subproblem
 * 
 * Triple solution, triple revolute joint geometric IK subproblem given by
 * @f$ e\,^{\hat{\xi_1}\,{\theta_1}} \cdot e\,^{\hat{\xi_2}\,{\theta_2}} \cdot e\,^{\hat{\xi_3}\,{\theta_3}} \cdot p = k @f$
 * (rotation around three subsequent axes (one skew + two parallel) applied to a point,
 * see @cite pardosgotor2022str).
 */
class PardosGotorSeven : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param exp1 First POE term.
     * @param exp2 Second POE term.
     * @param exp3 Third POE term.
     * @param p Characteristic point.
     */
    PardosGotorSeven(const MatrixExponential & exp1, const MatrixExponential & exp2, const MatrixExponential & exp3, const KDL::Vector & p);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override;

    int solutions() const override
    { return 4; }

    const char * describe() const override
    { return "PG7"; }

private:
    const MatrixExponential exp1, exp2, exp3;
    const KDL::Vector p, axesCross, axesCross_inverted;
    const KDL::Rotation axisPow1, axisPow2;
    const double axesDot;
};

/**
 * @ingroup ScrewTheoryLibd
 *
 * @brief Eighth Pardos-Gotor subproblem
 *
 * Triple solution, triple revolute joint geometric IK subproblem given by
 * @f$ e\,^{\hat{\xi_1}\,{\theta_1}} \cdot e\,^{\hat{\xi_2}\,{\theta_2}} \cdot e\,^{\hat{\xi_3}\,{\theta_3}} \cdot Hp = Hk @f$
 * (rotation around three subsequent parallel axes applied to a pose or coordinate system,
 * see @cite pardosgotor2022str).
 */
class PardosGotorEight : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param exp1 First POE term.
     * @param exp2 Second POE term.
     * @param exp3 Third POE term.  
     * @param p Characteristic point.
     * @param firstID ID of first axis
     * @param lastID ID of third axis
     * @param poe POE terms
     */
    PardosGotorEight(const MatrixExponential & exp1, const MatrixExponential & exp2, const MatrixExponential & exp3, const KDL::Vector & p, const int firstID, const int lastId, const PoeExpression poe);


    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions, const KDL::Frame & H_S_T, const KDL::JntArray & c_solutions, const KDL::Frame & H_S_T_0) const override;

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override
    {
        return false;
    }
    int solutions() const override
    { return 2; }//two sets of triple solutions

    const char * describe() const override
    { return "PG8"; }

private:
    const MatrixExponential exp1, exp2, exp3;
    const KDL::Vector p, n;
    const KDL::Rotation axisPow;
    const int firstID, lastID;
    const PoeExpression poe;
};

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief Third Pardos-Gotor and First Paden Kahan subproblems
 *
 * Dual solution, single prismatic joint geometric IK subproblem given by
 * @f$ e\,^{\hat{\xi}\,{\theta}} \cdot p = k_2 @f$
 * (combination of subproblems pg3 and pk1 used to solve the UR robot family,
 * see @cite pardosgotor2022str).
 */
class PardosGotorThreePadenKahanOne : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param exp_pk1 axis used for the pk1 subproblem.
     * @param exp_pg3 axis used for the pg3 subproblem.
     * @param p First characteristic point.
     * @param k Second characteristic point.
     */
    PardosGotorThreePadenKahanOne(const MatrixExponential & exp, const MatrixExponential & exp_pk1, const KDL::Vector & p, const KDL::Vector & k);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override;

    int solutions() const override
    { return 2; }

    const char * describe() const override
    { return "PG3+PK1"; }

private:
    const MatrixExponential exp_pg3, exp_pk1;
    const KDL::Vector p, k;
    const KDL::Rotation axisPow;
};

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief Algebraic approach
 *
 * Single solution, algebraic approach given by
 * @f$ H_{ST}(\theta) = e\,^{\hat{\xi_1}\,{\theta_1}} \cdot e\,^{\hat{\xi_2}\,{\theta_2}} \cdot e\,^{\hat{\xi_3}\,{\theta_3}} \cdot e\,^{\hat{\xi_4}\,{\theta_4}} \cdot e\,^{\hat{\xi_5}\,{\theta_5}} \cdot e\,^{\hat{\xi_6}\,{\theta_6}} \cdot H_{ST}(0) @f$
 * (algebraic approach used to solve the UR robot family,
 * see @cite pardosgotor2022str).
 */
class Algebraic_UR : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param q1 Index for the first joint.
     * @param q2  Index for the second joint.
     *
     */
    Algebraic_UR(int j1, int j2);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions, const KDL::Frame & H_S_T, const KDL::JntArray & c_solutions) const override;

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions, const KDL::Frame & H_S_T, const KDL::JntArray & c_solutions, const KDL::Frame & H_S_T_0) const override
    {
        return solve(rhs, pointTransform, reference, solutions, H_S_T, c_solutions);
    }

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override
    { 
        return false; 
    }

    int solutions() const override
    { return 1; }

    const char * describe() const override
    { return "Algebraic"; }

private:
    const int j1, j2;
};

} // namespace roboticslab

#endif // __SCREW_THEORY_IK_SUBPROBLEMS_HPP__
