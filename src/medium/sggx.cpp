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
#include <preform/medium.hpp>
#include <layered-sqt/medium/sggx.hpp>

namespace ls {

// Microvolume SGGX specular phase.
typedef pr::microvolume_sggx_specular_phase<Float> 
        SggxSpecularPhase;

// Microvolume SGGX diffuse phase.
typedef pr::microvolume_sggx_diffuse_phase<Float> 
        SggxDiffusePhase;

Float SggxPhaseMedium::phase(
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>& wi) const
{
    // Is specular?
    if (type == Type::eSpecular) {

        // Use specular.
        SggxSpecularPhase func(
                Mat3<Float>{
                    {1, 0, 0},
                    {0, 1, 0},
                    {0, 0, 1}
                },
                Vec3<Float>{
                    Apara * Apara,
                    Apara * Apara,
                    Aperp * Aperp
                });
        return func.ps(wo, wi);
    }
    else {

        // Use diffuse.
        SggxDiffusePhase func(
                Mat3<Float>{
                    {1, 0, 0},
                    {0, 1, 0},
                    {0, 0, 1}
                },
                Vec3<Float>{
                    Apara * Apara,
                    Apara * Apara,
                    Aperp * Aperp
                });
        return func.ps(generateCanonical2(pcg), wo, wi);
    }
}

Vec3<Float> SggxPhaseMedium::phaseSample(
            Pcg32& pcg,
            Float& tau,
            const Vec3<Float>& wo) const
{
    (void) tau;

    if (type == Type::eSpecular) {

        // Use specular.
        SggxSpecularPhase func(
                Mat3<Float>{
                    {1, 0, 0},
                    {0, 1, 0},
                    {0, 0, 1}
                },
                Vec3<Float>{
                    Apara * Apara,
                    Apara * Apara,
                    Aperp * Aperp
                });

        return func.ps_sample(generateCanonical2(pcg), wo);
    }
    else {

        // Use diffuse.
        SggxDiffusePhase func(
                Mat3<Float>{
                    {1, 0, 0},
                    {0, 1, 0},
                    {0, 0, 1}
                },
                Vec3<Float>{
                    Apara * Apara,
                    Apara * Apara,
                    Aperp * Aperp
                });

        return func.ps_sample(
                    generateCanonical2(pcg), 
                    generateCanonical2(pcg), wo);
    }
}

// Initialize from argument string.
void SggxPhaseMedium::init(const std::string& arg)
{
    std::stringstream iss(arg);
    std::string str;

    try {
        // Read parameters.
        // Note std::stod() throws if string is invalid.
        while (iss >> str) {
            if (!str.compare(0, 6, "Apara=", 6)) {
                Apara = std::stod(str.substr(6));
            }
            else
            if (!str.compare(0, 6, "Aperp=", 6)) {
                Aperp = std::stod(str.substr(6));
            }
            else
            if (str == "type=Specular") {
                type = Type::eSpecular;
            }
            else
            if (str == "type=Diffuse") {
                type = Type::eDiffuse;
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
    if (!(Apara > 0)) {
        error_message = ": Apara non-positive";
    }
    else
    if (!(Aperp > 0)) {
        error_message = ": Aperp non-positive";
    }
    // Error?
    if (error_message) {
        throw 
            std::runtime_error(
            std::string(__PRETTY_FUNCTION__).append(error_message));
    }

    // Clamp to be safe?
    Apara = pr::fmax(Apara, 1e-2);
    Aperp = pr::fmax(Aperp, 1e-2);
}

} // namespace ls
