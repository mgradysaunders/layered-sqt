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
#include <layered-sqt/layer/microsurface_lambertian.hpp>

namespace ls {

// Lambertian BRDF.
typedef 
    pr::microsurface_lambertian_brdf<
        Float,
        pr::microsurface_trowbridge_reitz_slope,
        pr::microsurface_uniform_height> MicrosurfaceLambertianBrdf;

// BSDF.
Float MicrosurfaceLambertianBrdfLayer::bsdf(
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>& wi,
            Float* f_pdf) const
{
    // Surface.
    MicrosurfaceLambertianBrdf surf = {
        fR, // fT,
        Vec2<Float>{alpha, alpha}
    };

    // Use multiple scattering?
    if (use_multiple_scattering) {

        // Function to generate canonical random numbers.
        auto uk = [&pcg]() {
            return generateCanonical(pcg);
        };

        // Multiple-scattering version.
        Float f = 0;
        if (f_pdf) {
            *f_pdf = 0;
        }
        for (int iter = 0; 
                 iter < iter_count; iter++) {
            Float tmp_f_pdf;
            Float tmp_f = surf.fm(uk, wo, wi, 0, &tmp_f_pdf);
            f = 
            f + (tmp_f - f) / (iter + 1);
            if (f_pdf) {
                *f_pdf = 
                *f_pdf + (tmp_f_pdf - *f_pdf) / (iter + 1);
            }
        }
        return f;
    }
    else {

        // Single-scattering version.
        Float f = 0;
        for (int iter = 0;
                 iter < iter_count; iter++) {
            f = 
            f + (surf.fs(generateCanonical2(pcg), wo, wi) - f) / (iter + 1);
        }
        if (f_pdf) {
            *f_pdf = surf.fs_pdf(wo, wi);
        }
        return f;
    }
}

// BSDF sample.
Vec3<Float> MicrosurfaceLambertianBrdfLayer::bsdfSample(
            Pcg32& pcg, 
            Float& tau,
            const Vec3<Float>& wo) const
{
    // Surface.
    MicrosurfaceLambertianBrdf surf = {
        fR, // fT,
        Vec2<Float>{alpha, alpha}
    };

    // Use multiple scattering?
    if (use_multiple_scattering) {

        // Function to generate canonical random numbers.
        auto uk = [&pcg]() {
            return generateCanonical(pcg);
        };

        // Multiple-scattering version.
        int k; 
        Vec3<Float> wi = surf.fm_pdf_sample(uk, wo, k);

        // If not perfectly energy-conserving, then BSDF
        // is not identical to BSDF-PDF.
        if (!(fR == 1)) {

            // Update throughput.
            Float f_pdf;
            Float f = surf.fm(uk, wo, wi, 0, &f_pdf);
            tau *= f / f_pdf;
        }

        return wi;
    }
    else {

        // Single-scattering version.
        Vec3<Float> wi = 
        surf.fs_pdf_sample(
                generateCanonical2(pcg), wo);

        // Update throughput.
        Float f = 0;
        for (int iter = 0;
                 iter < iter_count; iter++) {
            f = 
            f + (surf.fs(generateCanonical2(pcg), wo, wi) - f) / (iter + 1);
        }
        tau *= f / surf.fs_pdf(wo, wi);

        return wi;
    }
}

// Is transmissive?
bool MicrosurfaceLambertianBrdfLayer::isTransmissive() const
{
    return false;
}

// Initialize from argument string.
void MicrosurfaceLambertianBrdfLayer::init(const std::string& arg)
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
            else
            if (!str.compare(0, 11, "iter_count=", 11)) {
                iter_count = std::stoi(str.substr(11));
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
    if (!(alpha > 0)) {
        error_message = ": alpha is non-positive";
    }
    else 
    if (!(iter_count > 0)) {
        error_message = ": iter_count is non-positive";
    }
    // Error?
    if (error_message) {
        throw 
            std::runtime_error(
            std::string(__PRETTY_FUNCTION__).append(error_message));
    }

    // Prevent outlandish computation time?
    if (iter_count > 128) {
        iter_count = 128;
        std::cerr << "from " << __PRETTY_FUNCTION__ << ": ";
        std::cerr << "Warning, clamping iter_count to maximum of 128\n";
    }
}

} // namespace ls
