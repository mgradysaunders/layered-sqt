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
#pragma once
#ifndef LAYERED_SQT_MEDIUM_HPP
#define LAYERED_SQT_MEDIUM_HPP

#include <layered-sqt/common.hpp>

namespace ls {

/**
 * @defgroup medium_interface Medium interface
 *
 * `<layered-sqt/medium.hpp>`
 */
/**@{*/

/**
 * @brief Medium interface.
 */
class Medium
{
public:

    /**
     * @brief Default constructor.
     */
    Medium() = default;

    /**
     * @brief Destructor.
     */
    virtual ~Medium()
    {
    }

    /**
     * @brief Layer above.
     */
    const Layer* layer_above = nullptr;

    /**
     * @brief Layer below.
     */
    const Layer* layer_below = nullptr;

    /**
     * @brief Refractive index @f$ \eta @f$.
     */
    Float eta = 1;

    /**
     * @brief Volume absorption coefficient @f$ \mu_a @f$.
     */
    Float mua = 0;

    /**
     * @brief Volume scattering coefficient @f$ \mu_s @f$.
     */
    Float mus = 0;

    /**
     * @brief Volume extinction coefficient @f$ \mu = \mu_a + \mu_s @f$.
     */
    Float mu = 0;

public:

    /**
     * @brief Initialize from argument string.
     *
     * @param[in] arg
     * Argument.
     */
    virtual 
    void init(const std::string& arg);

    /**
     * @brief Transmittance.
     *
     * @param[in] d
     * Distance.
     */
    Float transmittance(Float d) const;

    /**
     * @brief Transmittance sample.
     *
     * @param[inout] pcg
     * Generator.
     *
     * @param[inout] tau
     * Path throughput.
     *
     * @param[in] dmax
     * Distance maximum.
     *
     * @returns
     * Distance of scattering event, or `dmax` if no scattering event.
     */
    Float transmittanceSample(
            Pcg32& pcg, 
            Float& tau,
            Float dmax) const;

    /**
     * @brief Phase function.
     *
     * @param[inout] pcg
     * Generator.
     *
     * @param[in] wo
     * Outgoing direction.
     *
     * @param[in] wi
     * Incident direction.
     */
    virtual 
    Float phase(
            Pcg32& pcg, 
            const Vec3<Float>& wo,
            const Vec3<Float>& wi) const;

    /**
     * @brief Phase function sample.
     *
     * @param[inout] pcg
     * Generator.
     *
     * @param[inout] tau
     * Path throughput.
     *
     * @param[in] wo
     * Outgoing direction.
     *
     * @returns
     * Incident direction.
     */
    virtual 
    Vec3<Float> phaseSample(
            Pcg32& pcg,
            Float& tau,
            const Vec3<Float>& wo) const;
};

/**@}*/

} // namespace ls

#endif // #ifndef LAYERED_SQT_MEDIUM_HPP
