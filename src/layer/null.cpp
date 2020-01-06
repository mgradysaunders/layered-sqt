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
#include <layered-sqt/layer/null.hpp>

namespace ls {

// Initialize from argument string.
void NullBsdfLayer::init(const std::string& arg)
{
    std::stringstream iss(arg);
    std::string str;
    if (iss >> str) {
        // Runtime error.
        throw
            std::runtime_error(
            std::string(__PRETTY_FUNCTION__)
                .append(": invalid argument '").append(str).append("'"));
    }
}

// BSDF.
Float NullBsdfLayer::bsdf(
            Pcg32&,
            const Vec3<Float>&,
            const Vec3<Float>&,
            Float* f_pdf) const
{
    // This function should never be called, return zeros.
    if (f_pdf) {
        *f_pdf = 0;
    }
    return 0;
}

// BSDF sample.
Vec3<Float> NullBsdfLayer::bsdfSample(
            Pcg32&, 
            Float&,
            const Vec3<Float>& wo) const
{
    return -wo;
}

// Is transmissive?
bool NullBsdfLayer::isTransmissive() const
{
    return true;
}

// Is null?
bool NullBsdfLayer::isNull() const
{
    return true;
}

} // namespace ls
