/* Copyright (c) 2019-20 M. Grady Saunders
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 
 *   1. Redistributions of source code must retain the above
 *      copyright notice, this list of conditions and the following
 *      disclaimer.
 * 
 *   2. Redistributions in binary form must reproduce the above
 *      copyright notice, this list of conditions and the following
 *      disclaimer in the documentation and/or other materials
 *      provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
/*+-+*/
#include <preform/medium.hpp>
#include <layered-sqt/medium/henyey_greenstein2.hpp>

namespace ls {

// Henyey-Greenstein phase.
typedef pre::hg_phase_stack<Float, 2> HenyeyGreenstein2Phase;

// Phase function.
Float HenyeyGreenstein2Medium::phase(
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>& wi) const
{
    (void) pcg;
    Vec2<Float> g = {g0, g1};
    Vec2<Float> w = {1 - b, b};
    return HenyeyGreenstein2Phase(g, w).ps(wo, wi);
}

// Phase function sample.
Vec3<Float> HenyeyGreenstein2Medium::phaseSample(
            Pcg32& pcg,
            Float& tau,
            const Vec3<Float>& wo) const
{
    (void) tau; // tau *= 1
    Vec2<Float> g = {g0, g1};
    Vec2<Float> w = {1 - b, b};
    return HenyeyGreenstein2Phase(g, w).ps_sample(generateCanonical2(pcg), wo);
}

// Initialize from argument string.
void HenyeyGreenstein2Medium::init(const std::string& arg)
{
    std::stringstream iss(arg);
    std::string str;

    try {
        // Read parameters.
        // Note std::stod() throws if string is invalid.
        while (iss >> str) {
            if (!str.compare(0, 3, "g0=", 3)) {
                g0 = std::stod(str.substr(3));
            }
            else
            if (!str.compare(0, 3, "g1=", 3)) {
                g1 = std::stod(str.substr(3));
            }
            else
            if (!str.compare(0, 2, "b=", 2)) {
                b = std::stod(str.substr(2));
            }
            else {
                // Trigger catch block.
                throw std::exception();
            }
        }
    }
    catch (...) {
        // Runtime error.
        throw
            std::runtime_error(
            std::string(__PRETTY_FUNCTION__)
                .append(": invalid argument '").append(str).append("'"));
    }

    // Check.
    const char* error_message = nullptr;
    if (!(g0 > -1 && g0 < 1)) {
        error_message = ": g0 is outside (-1, 1)";
    }
    else
    if (!(g1 > -1 && g1 < 1)) {
        error_message = ": g1 is outside (-1, 1)";
    }
    else
    if (!(b >= 0 && b <= 1)) {
        error_message = ": b is outside [0, 1]";
    }
    // Error?
    if (error_message) {
        throw 
            std::runtime_error(
            std::string(__PRETTY_FUNCTION__).append(error_message));
    }
}

} // namespace ls
