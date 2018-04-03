//
// Created by Jennifer Hays on 3/27/2018
//

#include "brmcpotential.h"
//#include <cmath>

#include <array>

namespace plugin
{

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

    output.energy = alpha * R/target;
    // Direction of force is ill-defined when v == v0
    if (R != 0)
    {
        // For harmonic: output.force = k * (double(R0)/R - 1.0)*rdiff;
        // For BRMC: outpu.force = - alpha/target * (unit vector in direction v-v0).
        output.force = (alpha/target/double(R)) * rdiff; // Why is there a double cast here?
    }

//    history.emplace_back(magnitude - R0);
    return output;
}

gmx::PotentialPointData BRMCRestraint::evaluate(gmx::Vector r1,
                                                 gmx::Vector r2,
                                                 double t)
{
    return calculate(r1, r2, t);
}

std::vector<unsigned long int> BRMCRestraint::sites() const
{
    return {site1_, site2_};
}

} // end namespace plugin
