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
#include <layered-sqt/layer/dielectric_microsurface.hpp>

namespace ls {

// Dielectric microsurface.
typedef 
    pr::microsurface_dielectric_bsdf<
        Float,
        pr::microsurface_trowbridge_reitz_slope,
        pr::microsurface_uniform_height> DielectricMicrosurfaceBsdf;

// BSDF.
Float DielectricMicrosurfaceBsdfLayer::bsdf(
        Pcg32& pcg,
        const Vec3<Float>& wo,
        const Vec3<Float>& wi) const
{
    // Surface.
    DielectricMicrosurfaceBsdf surf = {
        kR, kT, 
        this->medium_above->eta /
        this->medium_below->eta,
        Vec2<Float>{alpha, alpha}
    };

    if (fm_iters <= 0) {

        // Single-scattering version.
        return surf.fs(wo, wi);
    }
    else {

        // Function to generate canonical random numbers.
        auto uk = [&pcg]() {
            return generateCanonical(pcg);
        };

        // Multiple-scattering version.
        return surf.fm(uk, wo, wi, 0, fm_iters);
    }
}

// BSDF-PDF.
Float DielectricMicrosurfaceBsdfLayer::bsdfPdf(
        Pcg32& pcg,
        const Vec3<Float>& wo,
        const Vec3<Float>& wi) const
{
    // Surface.
    DielectricMicrosurfaceBsdf surf = {
        kR, kT, 
        this->medium_above->eta /
        this->medium_below->eta,
        Vec2<Float>{alpha, alpha}
    };

    if (fm_iters <= 0) {

        // Single-scattering version.
        return surf.fs_pdf(wo, wi);
    }
    else {

        // Function to generate canonical random numbers.
        auto uk = [&pcg]() {
            return generateCanonical(pcg);
        };

        // Multiple-scattering version.
        return surf.fm_pdf(uk, wo, wi, 0, fm_iters);
    }
}

// BSDF-PDF sample.
Vec3<Float> DielectricMicrosurfaceBsdfLayer::bsdfPdfSample(
        Pcg32& pcg, 
        Float& tau,
        const Vec3<Float>& wo) const
{
    // Surface.
    DielectricMicrosurfaceBsdf surf = {
        kR, kT, 
        this->medium_above->eta /
        this->medium_below->eta,
        Vec2<Float>{alpha, alpha}
    };

    if (fm_iters <= 0) {

        // Single-scattering version.
        Vec3<Float> wi = 
        surf.fs_pdf_sample(
                generateCanonical(pcg),
                generateCanonical2(pcg), wo);

        // Update throughput.
        tau *= 
            surf.fs(wo, wi) / 
            surf.fs_pdf(wo, wi);

        return wi;
    }
    else {
        // Function to generate canonical random numbers.
        auto uk = [&pcg]() {
            return generateCanonical(pcg);
        };

        // Multiple-scattering version.
        int k; 
        Vec3<Float> wi = surf.fm_pdf_sample(uk, wo, k);

        // If not perfectly energy-conserving, then BSDF
        // is not identical to BSDF-PDF.
        if (!(kR == 1 &&
              kT == 1)) {

            // Update throughput.
            Float fm_pdf;
            Float fm = surf.fm(uk, wo, wi, 0, fm_iters, &fm_pdf);
            tau *= fm / fm_pdf;
        }

        return wi;
    }
}

// Is transmissive?
bool DielectricMicrosurfaceBsdfLayer::isTransmissive() const
{
    return kT > 0;
}

// Load.
void DielectricMicrosurfaceBsdfLayer::load(const std::string& strln)
{
    std::stringstream iss(strln);
    std::string str;

    try {
        // Read parameters.
        // Note std::stod() throws if string is invalid.
        while (iss >> str) {
            if (!str.compare(0, 3, "kR=", 3)) {
                kR = std::stod(str.substr(3));
            }
            else 
            if (!str.compare(0, 3, "kT=", 3)) {
                kT = std::stod(str.substr(3));
            }
            else 
            if (!str.compare(0, 6, "alpha=", 6)) {
                alpha = std::stod(str.substr(6));
            }
            else
            if (!str.compare(0, 9, "fm_iters=", 9) ) {
                fm_iters = std::stoi(str.substr(9));
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
                .append(": invalid parameter '").append(str).append("'"));
    }

    // Check.
    const char* error_message = nullptr;
    if (!(kR >= 0 && kR <= 1)) {
        error_message = ": kR is outside [0, 1]";
    }
    else 
    if (!(kT >= 0 && kT <= 1)) {
        error_message = ": kT is outside [0, 1]";
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

    // Prevent outlandish computation time?
    if (fm_iters > 128) {
        fm_iters = 128;
        std::cerr << "from " << __PRETTY_FUNCTION__ << ": ";
        std::cerr << "Warning, clamping fm_iters to maximum value of 128\n";
    }
}

} // namespace ls
