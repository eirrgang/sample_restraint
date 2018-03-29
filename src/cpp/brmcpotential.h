//
// Created by Jennifer Hays on 3/28/2018
//

#ifndef GROMACS_BRMCPOTENTIAL_H
#define GROMACS_BRMCPOTENTIAL_H

#include <iostream>

#include "gmxapi/gromacsfwd.h"
#include "gmxapi/md/mdmodule.h"

#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/real.h"

namespace plugin
{

class BRMC
{
    public:
        BRMC(real equilibrium, real couplingconstant) :
            target{equilibrium},
            alpha{couplingconstant}
        {};

        BRMC() :
            BRMC{0.0, 0.0}
        {};

        // Allow easier automatic generation of bindings.
        struct input_param_type {
            float whateverIwant;
        };

        struct output_type
        {};

        /*!
         * \brief Calculate BRMC force on particle at position v in reference to position v0.
         *
         * \param v position at which to evaluate force
         * \param v0 position of BRMC bond reference
         * \return F = -k ((v - v0)/|v - v0| - R0);
         *
         * R0 == 1.0 is the equilibrium distance in the BRMC potential.
         * k == 1.0 is the spring constant.
         *
         * In the case of a pair of BRMCally bonded particles, the force on particle i is evaluated with particle j as
         * the reference point with
         * \code
         * auto force = calculateForce(r_i, r_j);
         * \endcode
         *
         * The force on particle j is the opposite as the force vector for particle i. E.g.
         * \code
         * assert(-1 * force, calculateForce(r_j, r_i));
         * \endcode
         */
        gmx::PotentialPointData calculate(gmx::Vector v,
                                          gmx::Vector v0,
                                          gmx_unused double t);

        // Cache of historical distance data. Not thread safe
//        std::vector<float> history{};

        // The class will either be inherited as a mix-in or inherit a CRTP base class. Either way, it probably needs proper virtual destructor management.
        virtual ~BRMC() {
//            for (auto&& distance: history)
//            {
//                std::cout << distance << "\n";
//            }
//            std::cout << std::endl;
        }

    private:
        // set equilibrium separation distance
        // TODO: be clearer about units
        real target;
        // set spring constant
        // TODO: be clearer about units
        real alpha;
};

// implement IRestraintPotential in terms of BRMC
// To be templated and moved.
class BRMCRestraint : public ::gmx::IRestraintPotential, private BRMC
{
    public:
        BRMCRestraint(unsigned long int site1,
                          unsigned long int site2,
                          real target,
                          real alpha) :
            BRMC{target, alpha},
            site1_{site1},
            site2_{site2}
        {};

        std::vector<unsigned long int> sites() const override;

        // \todo provide this facility automatically
        gmx::PotentialPointData evaluate(gmx::Vector r1,
                                         gmx::Vector r2,
                                         double t) override;

    private:
        unsigned long int site1_{0};
        unsigned long int site2_{0};
};

class BRMCModule : public gmxapi::MDModule
{
    public:
        using param_t = BRMC::input_param_type;

        BRMCModule(unsigned long int site1,
                       unsigned long int site2,
                       real target,
                       real alpha)
        {
            site1_ = site1;
            site2_ = site2;
            target_ = target;
            alpha_ = alpha;
        }


        const char *name() override
        {
            return "BRMCModule";
        }

        /*!
         * \brief implement gmxapi::MDModule::getRestraint()
         *
         * \return Handle to configured library object.
         */
        std::shared_ptr<gmx::IRestraintPotential> getRestraint() override
        {
            auto restraint = std::make_shared<BRMCRestraint>(site1_, site2_, target_, alpha_);
            return restraint;
        }

        /*!
         * \brief Set restraint parameters.
         *
         * \todo generalize this
         * \param site1
         * \param site2
         * \param k
         * \param R0
         */
        void setParams(unsigned long int site1,
                        unsigned long int site2,
                        real target,
                        real alpha)
        {
            site1_ = site1;
            site2_ = site2;
            target_ = target;
            alpha_ = alpha;
        }

    private:
        unsigned long int site1_;
        unsigned long int site2_;
        real target_;
        real alpha_;
};

} // end namespace plugin

#endif //GROMACS_BRMCPOTENTIAL_H
