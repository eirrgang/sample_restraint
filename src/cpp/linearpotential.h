//
// Created by Jennifer Hays on 5/3/2018
//

#ifndef GROMACS_LINEARPOTENTIAL_H
#define GROMACS_LINEARPOTENTIAL_H

#include <iostream>

#include "gmxapi/gromacsfwd.h"
#include "gmxapi/md/mdmodule.h"

#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/real.h"

namespace plugin
{

class Linear
{
    public:
        Linear(real equilibrium, real couplingconstant) :
            R0{equilibrium},
            k{couplingconstant}
        {};

        Linear() :
            Linear{0.0, 0.0}
        {};

        // Allow easier automatic generation of bindings.
        struct input_param_type {
            float whateverIwant;
        };

        struct output_type
        {};

        /*!
         * \brief Calculate Linear force on particle at position v in reference to position v0.
         *        U = k |v-v0|
         * \param v position at which to evaluate force
         * \param v0 position of Linear bond reference
         * \return F = -k (v - v0)/|v - v0|;
         *
         * R0 == 1.0 is the equilibrium distance in the Linear potential.
         * k == 1.0 is the spring constant.
         *
         * In the case of a pair of Linearally bonded particles, the force on particle i is evaluated with particle j as
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
        virtual ~Linear() {
//            for (auto&& distance: history)
//            {
//                std::cout << distance << "\n";
//            }
//            std::cout << std::endl;
        }

    private:
        // set equilibrium separation distance
        // TODO: be clearer about units
        real R0;
        // set spring constant
        // TODO: be clearer about units
        real k;
};

// implement IRestraintPotential in terms of Linear
// To be templated and moved.
class LinearRestraint : public ::gmx::IRestraintPotential, private Linear
{
    public:
        LinearRestraint(unsigned long int site1,
                          unsigned long int site2,
                          real R0,
                          real k) :
            Linear{R0, k},
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

class LinearModule : public gmxapi::MDModule
{
    public:
        using param_t = Linear::input_param_type;

        LinearModule(unsigned long int site1,
                       unsigned long int site2,
                       real R0,
                       real k)
        {
            site1_ = site1;
            site2_ = site2;
            R0_ = R0;
            k_ = k;
        }


        const char *name() override
        {
            return "LinearModule";
        }

        /*!
         * \brief implement gmxapi::MDModule::getRestraint()
         *
         * \return Handle to configured library object.
         */
        std::shared_ptr<gmx::IRestraintPotential> getRestraint() override
        {
            auto restraint = std::make_shared<LinearRestraint>(site1_, site2_, R0_, k_);
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
                        real R0,
                        real k)
        {
            site1_ = site1;
            site2_ = site2;
            R0_ = R0;
            k_ = k;
        }

    private:
        unsigned long int site1_;
        unsigned long int site2_;
        real R0_;
        real k_;
};

} // end namespace plugin

#endif //GROMACS_LINEARPOTENTIAL_H
