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

#include "make_unique.h"

namespace plugin
{
// Stop-gap for cross-language data exchange pending SharedData implementation and inclusion of Eigen.
// Adapted from pybind docs.
    template<class T>
    class Matrix {
    public:
        Matrix(size_t rows, size_t cols) :
                rows_(rows),
                cols_(cols),
                data_(rows_*cols_, 0)
        {
        }

        explicit Matrix(std::vector<T>&& captured_data) :
                rows_{1},
                cols_{captured_data.size()},
                data_{std::move(captured_data)}
        {
        }

        std::vector<T> *vector() { return &data_; }
        T* data() { return data_.data(); };
        size_t rows() const { return rows_; }
        size_t cols() const { return cols_; }
    private:
        size_t rows_;
        size_t cols_;
        std::vector<T> data_;
    };

// Defer implicit instantiation to ensemblepotential.cpp
    extern template class Matrix<double>;

/*!
 * \brief An active handle to ensemble resources provided by the Context.
 *
 * The semantics of holding this handle aren't determined yet, but it should be held as briefly as possible since it
 * may involve locking global resources or preventing the simulation from advancing. Basically, though, it allows the
 * Context implementation flexibility in how or where it provides services.
 */
    class EnsembleResourceHandle
    {
    public:
        /*!
         * \brief Ensemble reduce.
         *
         * For first draft, assume an all-to-all sum. Reduce the input into the stored Matrix.
         * // Template later... \tparam T
         * \param data
         */
//        void reduce(const Matrix<double>& input);

        /*!
         * \brief Ensemble reduce.
         * \param send Matrices to be summed across the ensemble using Context resources.
         * \param receive destination of reduced data instead of updating internal Matrix.
         */
        void reduce(const Matrix<double> &send,
                    Matrix<double> *receive) const;

        /*!
         * \brief Apply a function to each input and accumulate the output.
         *
         * \tparam I Iterable.
         * \tparam T Output type.
         * \param iterable iterable object to produce inputs to function
         * \param output structure that should be present and up-to-date on all ranks.
         * \param function map each input in iterable through this function to accumulate output.
         */
//        template<typename I, typename T>
//        void map_reduce(const I& iterable, T* output, void (*function)(double, const PairHist&, PairHist*));

        // to be abstracted and hidden...
        const std::function<void(const Matrix<double>&, Matrix<double>*)>* _reduce;
    };

/*!
 * \brief Reference to workflow-level resources managed by the Context.
 *
 * Provides a connection to the higher-level workflow management with which to access resources and operations. The
 * reference provides no resources directly and we may find that it should not extend the life of a Session or Context.
 * Resources are accessed through Handle objects returned by member functions.
 */
    class EnsembleResources
    {
    public:
        explicit EnsembleResources(std::function<void(const Matrix<double>&, Matrix<double>*)>&& reduce) :
                reduce_(reduce)
        {};

        EnsembleResourceHandle getHandle() const;

    private:
//        std::shared_ptr<Matrix> _matrix;
        std::function<void(const Matrix<double>&, Matrix<double>*)> reduce_;
    };

/*!
 * \brief Template for MDModules from restraints.
 *
 * \tparam R a class implementing the gmx::IRestraintPotential interface.
 */
    template<class R>
    class RestraintModule : public gmxapi::MDModule // consider names
    {
    public:
        using param_t = typename R::input_param_type;

        RestraintModule(std::string name,
                        std::vector<unsigned long int> sites,
                        const typename R::input_param_type& params,
                        std::shared_ptr<EnsembleResources> resources) :
                sites_{std::move(sites)},
                params_{params},
                resources_{std::move(resources)},
                name_{std::move(name)}
        {

        };

        ~RestraintModule() override = default;

        // \todo make member function const
        const char *name() override
        {
            return name_.c_str();
        }

        std::shared_ptr<gmx::IRestraintPotential> getRestraint() override
        {
            auto restraint = std::make_shared<R>(sites_, params_, resources_);
            return restraint;
        }

    private:
        std::vector<unsigned long int> sites_;
        param_t params_;

        // Need to figure out if this is copyable or who owns it.
        std::shared_ptr<EnsembleResources> resources_;

        const std::string name_;
    };

    struct brmc_input_param_type
    {
        /// learned coupling constant
        double alpha{0};
        double alpha_prev{0};

        /// keep track of mean and variance
        double mean{0};
        double variance{0};

        /// parameters for training coupling constant (Adagrad)
        double A{0};
        double tau{0};
        double g{0};
        double gsqrsum{0};
        double eta{0};
        bool converged{0};

        /// target distance
        double target{0};

        /// Number of samples to store during each tau window.
        unsigned int nSamples{0};
        unsigned int currentSample{0};
        double samplePeriod{0};
        double nextUpdateTime{0};

    };

// \todo We should be able to automate a lot of the parameter setting stuff
// by having the developer specify a map of parameter names and the corresponding type, but that could get tricky.
// The statically compiled fast parameter structure would be generated with a recursive variadic template
// the way a tuple is. ref https://eli.thegreenplace.net/2014/variadic-templates-in-c/

    std::unique_ptr<brmc_input_param_type>
    makeBRMCParams(double A,
                   double tau,
                   double target,
                   unsigned int nSamples)
    {
        using gmx::compat::make_unique;
        auto params = make_unique<brmc_input_param_type>();
        params->A = A;
        params->tau = tau;
        params->target = target;
        params->nSamples = nSamples;

        return params;
    };

class BRMC
{
    public:
        using input_param_type = brmc_input_param_type;

//        EnsembleHarmonic();

        explicit BRMC(const input_param_type &params);

        BRMC(double alpha,
             double alpha_prev,
             double mean,
             double variance,
             double A,
             double tau,
             double g,
             double gsqrsum,
             double eta,
             bool converged,
             double target,
             unsigned int nSamples,
             double samplePeriod);

        // If dispatching this virtual function is not fast enough, the compiler may be able to better optimize a free
        // function that receives the current restraint as an argument.
        gmx::PotentialPointData calculate(gmx::Vector v,
                                          gmx::Vector v0,
                                          gmx_unused double t);

        // An update function to be called on the simulation master rank/thread periodically by the Restraint framework.
        void callback(gmx::Vector v,
                      gmx::Vector v0,
                      double t,
                      const EnsembleResources &resources);

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
        /// learned coupling constant
        double alpha_;
        double alpha_prev_;

        /// keep track of mean and variance
        double mean_;
        double variance_;

        /// parameters for training coupling constant (Adagrad)
        double A_;
        double tau_;
        double g_;
        double gsqrsum_;
        double eta_;
        bool converged_;

        /// target distance
        double target_;

        /// Number of samples to store during each window.
        unsigned int nSamples_;
        unsigned int currentSample_;
        double samplePeriod_;
        double nextSampleTime_;

        double windowStartTime_;
        double nextUpdateTime_;
};

// implement IRestraintPotential in terms of BRMC
// To be templated and moved.
class BRMCRestraint : public ::gmx::IRestraintPotential, private BRMC
{
    public:
        using BRMC::input_param_type;
        BRMCRestraint(const std::vector<unsigned long> &sites,
                      const input_param_type &params,
                      std::shared_ptr<EnsembleResources> resources):
            BRMC(params),
            sites_{sites},
            resources_{std::move(resources)}
        {}

    std::vector<unsigned long int> sites() const override
    {
        return sites_;
    }

    gmx::PotentialPointData evaluate(gmx::Vector r1,
                                     gmx::Vector r2,
                                     double t) override
    {
        return calculate(r1, r2, t);
    };


    // An update function to be called on the simulation master rank/thread periodically by the Restraint framework.
    void update(gmx::Vector v,
                gmx::Vector v0,
                double t) override
    {
        // Todo: use a callback period to mostly bypass this and avoid excessive mutex locking.
        callback(v,
                 v0,
                 t,
                 *resources_);
    };

    void setResources(std::unique_ptr<EnsembleResources>&& resources)
    {
        resources_ = std::move(resources);
    }

private:
    std::vector<unsigned long int> sites_;
//        double callbackPeriod_;
//        double nextCallback_;
    std::shared_ptr<EnsembleResources> resources_;
};


// Just declare the template instantiation here for client code.
// We will explicitly instantiate a definition in the .cpp file where the input_param_type is defined.
    extern template class RestraintModule<BRMCRestraint>;
} // end namespace plugin

#endif //GROMACS_BRMCPOTENTIAL_H
