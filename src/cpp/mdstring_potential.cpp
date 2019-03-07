/*! \file
 * \brief Code to implement the potential declared in mdstring_potential.h
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 */

#include "mdstring_potential.h"

#include <cassert>
#include <cmath>

#include <memory>
#include <vector>

#include "gmxapi/context.h"
#include "gmxapi/session.h"
#include "gmxapi/md/mdsignals.h"

#include "sessionresources.h"

namespace plugin
{

/*!
 * \brief Discretize a density field on a grid.
 *
 * Apply a Gaussian blur when building a density grid for a list of values.
 * Normalize such that the area under each sample is 1.0/num_samples.
 */
class BlurToGrid
{
    public:
        /*!
         * \brief Construct the blurring functor.
         *
         * \param low The coordinate value of the first grid point.
         * \param gridSpacing Distance between grid points.
         * \param sigma Gaussian parameter for blurring inputs onto the grid.
         */
        BlurToGrid(double low,
                   double gridSpacing,
                   double sigma) :
            low_{low},
            binWidth_{gridSpacing},
            sigma_{sigma}
        {
        };

        /*!
         * \brief Callable for the functor.
         *
         * \param samples A list of values to be blurred onto the grid.
         * \param grid Pointer to the container into which to accumulate a blurred histogram of samples.
         *
         * Example:
         *
         *     # Acquire 3 samples to be discretized with blurring.
         *     std::vector<double> someData = {3.7, 8.1, 4.2};
         *
         *     # Create an empty grid to store magnitudes for points 0.5, 1.0, ..., 10.0.
         *     std::vector<double> histogram(20, 0.);
         *
         *     # Specify the above grid and a Gaussian parameter of 0.8.
         *     auto blur = BlurToGrid(0.5, 0.5, 0.8);
         *
         *     # Collect the density grid for the samples.
         *     blur(someData, &histogram);
         *
         */
        void operator()(const std::vector<double>& samples,
                        std::vector<double>* grid)
        {
            const auto nbins = grid->size();
            const double& dx{binWidth_};
            const auto num_samples = samples.size();

            const double denominator = 1.0 / (2 * sigma_ * sigma_);
            const double normalization = 1.0 / (num_samples * sqrt(2.0 * M_PI * sigma_ * sigma_));
            // We aren't doing any filtering of values too far away to contribute meaningfully, which
            // is admittedly wasteful for large sigma...
            for (size_t i = 0;i < nbins;++i)
            {
                double bin_value{0};
                const double bin_x{low_ + i * dx};
                for (const auto distance : samples)
                {
                    const double relative_distance{bin_x - distance};
                    const auto numerator = -relative_distance * relative_distance;
                    bin_value += normalization * exp(numerator * denominator);
                }
                grid->at(i) = bin_value;
            }
        };

    private:
        /// Minimum value of bin zero
        const double low_;

        /// Size of each bin
        const double binWidth_;

        /// Smoothing factor
        const double sigma_;
};

MDStringPotential::MDStringPotential(size_t nbins,
                                     double binWidth,
                                     double minDist,
                                     double maxDist,
                                     const std::vector<double>& experimental,
                                     unsigned int nSamples,
                                     double samplePeriod,
                                     unsigned int nWindows,
                                     double k,
                                     double sigma) :
    state_{nbins, binWidth, minDist, maxDist, experimental, nSamples, samplePeriod, nWindows, k, sigma}
{
    state_.histogram = std::move(std::vector<double>(nbins, 0.));
    state_.nextSampleTime = state_.samplePeriod;
    state_.nextWindowUpdateTime = state_.nSamples * state_.samplePeriod;
}

MDStringPotential::MDStringPotential(const input_param_type& params) :
    MDStringPotential(params.nBins,
                     params.binWidth,
                     params.minDist,
                     params.maxDist,
                     params.experimental,
                     params.nSamples,
                     params.samplePeriod,
                     params.nWindows,
                     params.k,
                     params.sigma)
{
}

//
//
// HERE is the (optional) function that updates the state of the restraint periodically.
// It is called before calculate() once per timestep per simulation (on the master rank of
// a parallelized simulation).
//
//
void MDStringPotential::callback(gmx::Vector v,
                                gmx::Vector v0,
                                double t,
                                const Resources& resources)
{
    const auto rdiff = v - v0;
    const auto Rsquared = dot(rdiff,
                              rdiff);
    const auto R = sqrt(Rsquared);

    // Store historical data every sample_period steps
    if (t >= state_.nextSampleTime)
    {
        state_.distanceSamples[state_.currentSample++] = R;
        state_.nextSampleTime = (state_.currentSample + 1) * state_.samplePeriod + state_.windowStartTime;
    };

    // Every nsteps:
    if (t >= state_.nextWindowUpdateTime)
    {
        // Reduce sampled data for this restraint in this simulation, applying a Gaussian blur to fill a grid.
        auto blur = BlurToGrid(0.0,
                               state_.binWidth,
                               state_.sigma);
        assert(state_.distanceSamples.size() == state_.nSamples);
        assert(state_.currentSample == state_.nSamples);
//        blur(state_.distanceSamples,
//             new_window.vector());
        // We can just do the blur locally since there aren't many bins. Bundling these operations for
        // all restraints could give us a chance at some parallelism. We should at least use some
        // threading if we can.

        // We request a handle each time before using resources to make error handling easier if there is a failure in
        // one of the mdstring member processes and to give more freedom to how resources are managed from step to step.
        auto ensemble = resources.getHandle();
        // Get global reduction (sum) and checkpoint.
//        ensemble.reduce(send_buffer,
//                        receive_buffer);

        // Note we do not have the integer timestep available here. Therefore, we can't guarantee that updates occur
        // with the same number of MD steps in each interval, and the interval will effectively lose digits as the
        // simulation progresses, so _update_period should be cleanly representable in binary. When we extract this
        // to a facility, we can look for a part of the code with access to the current timestep.
        state_.windowStartTime = t;
        state_.nextWindowUpdateTime = state_.nSamples * state_.samplePeriod + state_.windowStartTime;

        // Reset sample bufering.
        state_.currentSample = 0;
        // Reset sample times.
        state_.nextSampleTime = t + state_.samplePeriod;
    };

}


//
//
// HERE is the function that does the calculation of the restraint force.
//
//
gmx::PotentialPointData MDStringPotential::calculate(gmx::Vector v,
                                                    gmx::Vector v0,
                                                    double t)
{
    // This is not the vector from v to v0. It is the position of a site
    // at v, relative to the origin v0. This is a potentially confusing convention...
    const auto rdiff = v - v0;
    const auto Rsquared = dot(rdiff,
                              rdiff);
    const auto R = sqrt(Rsquared);


    // Compute output
    gmx::PotentialPointData output;
    // Energy not needed right now.
//    output.energy = 0;

    if (R != 0) // Direction of force is ill-defined when v == v0
    {

        double f{0};

        double f_scal{0};

        const size_t numBins = state_.histogram.size();
        double normConst = sqrt(2 * M_PI) * state_.sigma * state_.sigma * state_.sigma;

        for (size_t n = 0;n < numBins;n++)
        {
            const double x{n * state_.binWidth - R};
            const double argExp{-0.5 * x * x / (state_.sigma * state_.sigma)};
            f_scal += state_.histogram.at(n) * exp(argExp) * x / normConst;
        }
        f = -state_.k * f_scal;

        const auto magnitude = f / norm(rdiff);
        output.force = rdiff * static_cast<decltype(rdiff[0])>(magnitude);
    }
    return output;
}

std::unique_ptr<mdstring_data_t>
makeMDStringParams(size_t nbins,
                   double binWidth,
                   double minDist,
                   double maxDist,
                   const std::vector<double>& experimental,
                   unsigned int nSamples,
                   double samplePeriod,
                   unsigned int nWindows,
                   double k,
                   double sigma)
{
    using std::make_unique;
    auto params = make_unique<mdstring_data_t>();
    params->nBins = nbins;
    params->binWidth = binWidth;
    params->minDist = minDist;
    params->maxDist = maxDist;
    params->experimental = experimental;
    params->nSamples = nSamples;
    params->samplePeriod = samplePeriod;
    params->nWindows = nWindows;
    params->k = k;
    params->sigma = sigma;

    return params;
};

// Important: Explicitly instantiate a definition for the templated class declared in mdstringpotential.h.
// Failing to do this will cause a linker error.
template
class ::plugin::RestraintModule<Restraint<MDStringPotential>>;

} // end namespace plugin
