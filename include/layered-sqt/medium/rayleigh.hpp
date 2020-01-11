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
#ifndef LAYERED_SQT_MEDIUM_RAYLEIGH_HPP
#define LAYERED_SQT_MEDIUM_RAYLEIGH_HPP

#include <layered-sqt/medium.hpp>

namespace ls {

/**
 * @addtogroup mediums Mediums
 *
 * `<layered-sqt/medium/rayleigh.hpp>`
 */
/**@{*/

/**
 * @brief Rayleigh medium.
 */
class RayleighMedium final : public Medium
{
public:

    /**
     * @brief Default constructor.
     */
    RayleighMedium() = default;

    /**
     * @brief Depolarization factor @f$ \rho \in (-1, +1) @f$.
     */
    Float rho = 0;

public:

    /**
     * @copydoc Medium::init()
     *
     * Accepts arguments
     * - `rho=`(float).
     *
     * @throw std::runtime_error
     * If
     * - invalid argument, or
     * - `rho` is outside `(-1, 1)`.
     */
    void init(const std::string& arg);

    /**
     * @copydoc Medium::phase()
     */
    Float phase(
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>& wi) const;

    /**
     * @copydoc Medium::phaseSample()
     */
    Vec3<Float> phaseSample(
            Pcg32& pcg,
            Float& tau,
            const Vec3<Float>& wo) const;
};

/**@}*/

} // namespace ls

#endif // #ifndef LAYERED_SQT_MEDIUM_RAYLEIGH_HPP
