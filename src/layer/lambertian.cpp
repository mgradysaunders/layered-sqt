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
#include <layered-sqt/layer/lambertian.hpp>

namespace ls {

// BSDF.
Float LambertianBsdfLayer::bsdf(
        Pcg32&,
        const Vec3<Float>& wo,
        const Vec3<Float>& wi) const
{
    return pr::numeric_constants<Float>::M_1_pi() * pr::abs(wi[2]) * 
          (pr::signbit(wo[2]) == pr::signbit(wi[2]) ? fR : fT);
}

// BSDF-PDF.
Float LambertianBsdfLayer::bsdfPdf(
        Pcg32&,
        const Vec3<Float>& wo,
        const Vec3<Float>& wi) const
{
    return Vec3<Float>::cosine_hemisphere_pdf(pr::abs(wi[2])) *
          (pr::signbit(wo[2]) == pr::signbit(wi[2]) ? fR : fT) / (fR + fT);
}

// BSDF-PDF sample.
Vec3<Float> LambertianBsdfLayer::bsdfPdfSample(
        Pcg32& pcg, 
        Float& tau,
        const Vec3<Float>& wo) const
{
    Vec3<Float> wi = 
    Vec3<Float>::cosine_hemisphere_pdf_sample(generateCanonical2(pcg));
    if (generateCanonical(pcg) < fR / (fR + fT)) {
        wi[2] = pr::copysign(wi[2], +wo[2]);
    }
    else {
        wi[2] = pr::copysign(wi[2], -wo[2]);
    }

    // Update throughput.
    tau *= fR + fT;

    return wi;
}

// Is transmissive?
bool LambertianBsdfLayer::isTransmissive() const
{
    return fT > 0;
}

// Initialize from argument string.
void LambertianBsdfLayer::init(const std::string& arg)
{
    std::stringstream iss(arg);
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
    // Error?
    if (error_message) {
        throw 
            std::runtime_error(
            std::string(__PRETTY_FUNCTION__).append(error_message));
    }
}

} // namespace ls
