//
// Created by Eric Irrgang on 3/24/18.
//

#include "testingconfiguration.h"

#include <iostream>
#include <vector>

#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/arrayref.h"

#include "mdstring_potential.h"
#include "sessionresources.h"

#include <gtest/gtest.h>

using ::gmx::Vector;

namespace {

std::ostream& operator<<(std::ostream& stream, const Vector& vec)
{
    stream << "(" << vec[0] << "," << vec[1] << "," << vec[2] << ")";
    return stream;
}

TEST(PotentialPlugin, ForceCalc)
{
    const Vector zerovec = {0, 0, 0};
    // define some unit vectors
    const Vector e1{real(1), real(0), real(0)};
    const Vector e2{real(0), real(1), real(0)};
    const Vector e3{real(0), real(0), real(1)};

//    const real R0{1.0};
//    const real k{1.0};
//
//    // store temporary values long enough for inspection
//    Vector force{};
//
//    // Define a reference distribution with a triangular peak at the 1.0 bin.
//    const std::vector<double>
//        experimental{{0, 1, 0, 0, 0, 0, 0, 0, 0, 0}};
//
//
//    plugin::EnsembleHarmonic restraint{10, // nbins
//                                       1.0, // binWidth
//                                       0.0, // minDist
//                                       10.0, // maxDist
//                                       experimental, // experimental reference histogram
//                                       1, // nSamples
//                                       0.001, // samplePeriod
//                                       1, // nWindows
//                                       100., // k
//                                       1.0 // sigma
//    };
//
//    auto calculateForce =
//        [&restraint](const Vector& a, const Vector& b, double t)
//        {
//            return restraint.calculate(a,b,t).force;
//        };
//
//    // With the initial histogram (all zeros) the force should be zero no matter where the particles are.
//    ASSERT_EQ(static_cast<real>(0.0), norm(calculateForce(e1, e1, 0.)));
//    ASSERT_EQ(static_cast<real>(0.0), norm(calculateForce(e1, e2, 0.)));
//    ASSERT_EQ(static_cast<real>(0.0), norm(calculateForce(e1, static_cast<real>(-1)*e1, 0.)));
}

} // end anonymous namespace

int main(int argc, char* argv[])
{
::testing::InitGoogleTest(&argc, argv);
return RUN_ALL_TESTS();
}
