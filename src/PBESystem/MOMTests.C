// Copyright 2005, Google Inc.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// A sample program demonstrating using Google C++ testing framework.
//
// Author: wan@google.com (Zhanyong Wan)


// This sample shows how to write a more complex unit test for a class
// that has multiple member functions.
//
// Usually, it's a good idea to have one test for each method in your
// class.  You don't have to do that exactly, but it helps to keep
// your tests organized.  You may also throw in additional tests as
// needed.

#include "gtest/gtest.h"
#include "PBESystems-internal.H"

//This is in order to tunnel the data from argc and argv into OF routines.
int my_argc;
char** my_argv;

namespace Foam
{
namespace PBESystems
{
namespace internal
{

void expect_appr_eq(scalar expected, scalar value)
{
    EXPECT_LT(mag(expected - value) / expected, 0.01);
}

//Testing particular values
TEST(gammaTest, IntegerValues)
{
    //Values taken from Wikipedia (oracle test!)
    expect_appr_eq(1,gamma(1));
    expect_appr_eq(1,gamma(2));
    expect_appr_eq(2,gamma(3));
    expect_appr_eq(6,gamma(4));
    expect_appr_eq(24,gamma(5));
}

//Testing particular values
TEST(gammaTest, RealValues)
{
    //Values taken from Wikipedia (oracle test!)
    expect_appr_eq(sqrt(Foam::constant::mathematical::pi), gamma(0.5));
    expect_appr_eq(0.5 * sqrt(Foam::constant::mathematical::pi), gamma(1.5));
    expect_appr_eq(0.75 * sqrt(Foam::constant::mathematical::pi), gamma(2.5));
    expect_appr_eq(4.0/3.0 * sqrt(Foam::constant::mathematical::pi), gamma(-1.5));
    expect_appr_eq(-2.0 * sqrt(Foam::constant::mathematical::pi), gamma(-0.5));
}

//Testing for \gamma{t+1} = t \gamma{t}
TEST(gammaTest, BasicProperty)
{
    //Some random numbers
    scalar t[5] = {0.21312, -0.53, 1.23, -1.2, 11};
    for (int i = 0; i < 5; i++)
    {
        expect_appr_eq(gamma(t[i] + 1), t[i] * gamma(t[i]));
    }
}

}
}
}
int main(int argc, char** argv)
{    
    ::testing::InitGoogleTest(&argc, argv);
    my_argc = argc;
    my_argv = argv;
    return RUN_ALL_TESTS();
}

