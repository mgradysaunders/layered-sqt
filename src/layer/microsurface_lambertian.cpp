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
#include <preform/microsurface.hpp>
#include <layered-sqt/layer/microsurface_lambertian.hpp>

namespace ls {

// Microsurface from Lambertian BSDF.
typedef 
    pr::microsurface_lambertian_bsdf<
        Float,
        pr::microsurface_trowbridge_reitz_slope,
        pr::microsurface_uniform_height> MicrosurfaceLambertianBsdf;

// BSDF.
Float MicrosurfaceLambertianLayer::bsdf(
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>& wi,
            Float* f_pdf) const
{
    // Surface.
    MicrosurfaceLambertianBsdf surf = {
        fR, fT,
        Vec2<Float>{alpha, alpha}
    };

    return surf.fs(pcg, wo, wi, 0, 
                   use_multiple_scattering ? 0 : 1, f_pdf);
}

// BSDF sample.
Vec3<Float> MicrosurfaceLambertianLayer::bsdfSample(
            Pcg32& pcg, 
            Float& tau,
            const Vec3<Float>& wo) const
{
    // Surface.
    MicrosurfaceLambertianBsdf surf = {
        fR, fT,
        Vec2<Float>{alpha, alpha}
    };

    unsigned k; 
    Vec3<Float> wi = surf.fs_pdf_sample(pcg, wo, k);

    // If not perfectly energy-conserving, then BSDF
    // is not identical to BSDF-PDF.
    if (!(fR + fT == 1)) {

        // Update throughput.
        Float f_pdf;
        Float f = surf.fs(pcg, wo, wi, 0, 
                          use_multiple_scattering ? 0 : 1, &f_pdf);
        tau *= f / f_pdf;
    }

    return wi;
}

// Is transmissive?
bool MicrosurfaceLambertianLayer::isTransmissive() const
{
    return false;
}

// Initialize from argument string.
void MicrosurfaceLambertianLayer::init(const std::string& arg)
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
            if (!str.compare(0, 3, "fT=", 3)) {
                fT = std::stod(str.substr(3));
            }
            else 
            if (!str.compare(0, 6, "alpha=", 6)) {
                alpha = std::stod(str.substr(6));
            }
            else
            if (!str.compare(0, 24, "use_multiple_scattering=", 24)) {
                if (str.substr(24) == "true") {
                    use_multiple_scattering = true;
                }
                else
                if (str.substr(24) == "false") {
                    use_multiple_scattering = false;
                }
                else {
                    // Trigger catch block.
                    throw std::exception();
                }
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
    if (!(fT >= 0 && fT <= 1)) {
        error_message = ": fT is outside [0, 1]";
    }
    else 
    if (!(fR + fT <= 1)) {
        error_message = ": fR + fT is greater than 1";
    }
    else 
    if (!(alpha > 0)) {
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
