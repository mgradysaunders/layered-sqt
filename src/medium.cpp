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
#include <sstream>
#include <layered-sqt/medium.hpp>

namespace ls {

// Transmittance.
Float Medium::transmittance(Float d) const
{
    assert(d >= 0);
    if (d == 0 || mu == 0) {
        return 1;
    }
    else {
        return pr::exp(-mu * 
               pr::min(d, pr::numeric_limits<Float>::max()));
    }
}

// Transmittance sample.
Float Medium::transmittanceSample(
        Pcg32& pcg, 
        Float& tau,
        Float dmax) const 
{
    if (mu != 0) {

        // Sample distance.
        Float u = generateCanonical(pcg);
        Float d = (u < Float(0.5) ? -pr::log1p(-u) : -pr::log(1 - u)) / mu;

        // Scattering event occurred?
        if (d < dmax) {
            // Update throughput.
            tau *= mus / mu;

            // Return distance of scattering event.
            return d;
        }
        else {
            // Indicate no scattering event.
            return dmax;
        }
    }
    else {
        // Indicate no scattering event.
        return dmax;
    }
}

// Phase function.
Float Medium::phase(
        const Vec3<Float>& wo,
        const Vec3<Float>& wi) const
{
    (void) wo;
    (void) wi;
    throw 
        std::runtime_error(
        std::string(__PRETTY_FUNCTION__).append(": invalid call"));
}

// Phase function sample.
Vec3<Float> Medium::phaseSample(
        Pcg32& pcg,
        Float& tau,
        const Vec3<Float>& wo) const
{
    (void) pcg;
    (void) tau;
    (void) wo;
    throw 
        std::runtime_error(
        std::string(__PRETTY_FUNCTION__).append(": invalid call"));
}

// Initialize from argument string.
void Medium::init(const std::string& arg)
{
    (void) arg;
}

} // namespace ls
