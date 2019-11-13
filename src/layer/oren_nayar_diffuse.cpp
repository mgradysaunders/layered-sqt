/* Copyright (c) 2019 M. Grady Saunders
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
#include <preform/microsurface.hpp>
#include <layered-sqt/layer/oren_nayar_diffuse.hpp>

namespace ls {

typedef pr::oren_nayar_diffuse_brdf<Float> OrenNayarDiffuseBrdf;

// BSDF.
Float OrenNayarDiffuseBrdfLayer::bsdf(
            Pcg32&,
            const Vec3<Float>& wo,
            const Vec3<Float>& wi,
            Float* f_pdf) const
{
    if (pr::signbit(wo[2]) != 
        pr::signbit(wi[2])) {
        if (f_pdf) {
            *f_pdf = 0;
        }
        return 0;
    }
    else {
        OrenNayarDiffuseBrdf surf = {sigma};
        if (f_pdf) {
            *f_pdf = 
            surf.fs_pdf(wo, wi); 
        }
        return fR * surf.fs(wo, wi);
    }
}

// BSDF sample.
Vec3<Float> OrenNayarDiffuseBrdfLayer::bsdfSample(
            Pcg32& pcg, 
            Float& tau,
            const Vec3<Float>& wo) const
{
    OrenNayarDiffuseBrdf surf = {sigma};
    Vec3<Float> wi = surf.fs_pdf_sample(generateCanonical2(pcg), wo);
    tau *= fR * surf.fs(wo, wi) / 
                surf.fs_pdf(wo, wi);
    return wi;
}

// Is transmissive?
bool OrenNayarDiffuseBrdfLayer::isTransmissive() const
{
    return false;
}

// Initialize from argument string.
void OrenNayarDiffuseBrdfLayer::init(const std::string& arg)
{
    std::istringstream iss(arg);
    std::string str;

    try {
        // Read arguments.
        // Note std::stod() throws if string is invalid.
        while (iss >> str) {
            if (!str.compare(0, 3, "fR=", 3)) {
                fR = std::stod(str.substr(3));
            }
            else 
            if (!str.compare(0, 6, "sigma=", 6)) {
                sigma = std::stod(str.substr(6));
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
    if (!(fR >= 0 && fR <= 1)) {
        error_message = ": fR is outside [0, 1]";
    }
    else 
    if (!(sigma >= 0 && sigma <= 0.8)) {
        error_message = ": alpha is non-positive";
    }
    // Error?
    if (error_message) {
        throw 
            std::runtime_error(
            std::string(__PRETTY_FUNCTION__).append(error_message));
    }
}

} // namespace ls
