//
// Created by Jennifer Hays on 3/27/2018
//

#include "brmcpotential.h"
//#include <cmath>

#include <array>

namespace plugin
{

    template class ::plugin::Matrix<double>;

    BRMC::BRMC(double alpha,
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
               double samplePeriod,
               unsigned int currentSample,
               double windowStartTime):
    alpha_{alpha},
    alpha_prev_{alpha_prev},
    mean_{mean},
    variance_{variance},
    A_{A},
    tau_{tau},
    g_{g},
    gsqrsum_{gsqrsum},
    eta_{eta},
    converged_{converged},
    target_{target},
    nSamples_{nSamples},
    samplePeriod_{samplePeriod},
    nextSampleTime_{samplePeriod},
    nextUpdateTime_{nSamples*samplePeriod},
    currentSample_{currentSample},
    windowStartTime_{windowStartTime}
{};

    BRMC::BRMC(const input_param_type &params) :
            BRMC(params.alpha,
                 params.alpha_prev,
                 params.mean,
                 params.variance,
                 params.A,
                 params.tau,
                 params.g,
                 params.gsqrsum,
                 params.eta,
                 params.converged,
                 params.target,
                 params.nSamples,
                 params.samplePeriod,
                 params.currentSample,
                 params.windowStartTime)
    {}
    void BRMC::callback(gmx::Vector v, gmx::Vector v0, double t,
                                           const EnsembleResources &resources) {

        if (!converged_){
            auto rdiff = v - v0;
            const auto Rsquared = dot(rdiff,
                                      rdiff);
            const auto R = sqrt(Rsquared);
            if (t == 0){
                mean_ = R;
            }

            if (t >= nextSampleTime_){
                // update mean and variance
                int j = currentSample_+1;
                auto difference = (R - mean_);
                auto diffsqr = difference*difference;
//                printf("Difference, j: %f, %d\n", difference, j);
//                printf("Mean and variance: %f, %f\n", mean_, variance_);
                variance_ = variance_ + (j - 1) * diffsqr / j;
                mean_ = mean_ + difference / j;
//                printf("Difference, j: %f, %d\n", difference, j);
//                printf("Mean and variance: %f, %f\n", mean_, variance_);
                // Update next time to take a sample
                currentSample_++;
                nextSampleTime_ = (currentSample_ + 1)*samplePeriod_ + windowStartTime_;
            }

            if(t >= nextUpdateTime_){
                assert(currentSample_ == nSamples_);
                printf("Alpha: %f", alpha_);
                g_ = (1-mean_/target_)*variance_;
                eta_ = A_/sqrt(gsqrsum_);
                alpha_prev_ = alpha_;
                alpha_ = alpha_prev_ - eta_*g_;

                printf("alpha: %f", alpha_);
                // Reset mean and variance
                mean_ = R;
                variance_ = 0;
                windowStartTime_ = t;
                nextUpdateTime_ = nSamples_*samplePeriod_ + windowStartTime_;

                // Reset sample buffering.
                currentSample_ = 0;
                // Reset sample times.
                nextSampleTime_ = t + samplePeriod_;
            }
        }
    }


gmx::PotentialPointData BRMC::calculate(gmx::Vector v,
                                   gmx::Vector v0,
                                   gmx_unused double t)
{
    // Our convention is to calculate the force that will be applied to v.
    // An equal and opposite force is applied to v0.
    //auto rdiff = v - v0;
    auto rdiff = v0 - v; // Taking v0-v just let's me not apply a negative sign for output.force.
    const auto Rsquared = dot(rdiff, rdiff);
    const auto R = sqrt(Rsquared);
    // TODO: find appropriate math header and namespace

    // In White & Voth, the additional energy is alpha * f(r)/favg

    gmx::PotentialPointData output;

    output.energy = alpha_ * R/target_;
    // Direction of force is ill-defined when v == v0
    if (R != 0)
    {
        // For harmonic: output.force = k * (double(R0)/R - 1.0)*rdiff;
        // For BRMC: outpu.force = - alpha/target * (unit vector in direction v-v0).
        output.force = (alpha_/target_/double(R)) * rdiff; // Why is there a double cast here?
    }

//    history.emplace_back(magnitude - R0);
    return output;
}

    template class ::plugin::RestraintModule<BRMCRestraint>;
} // end namespace plugin
