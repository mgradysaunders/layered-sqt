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
#pragma once
#ifndef LAYERED_SQT_LAYER_HPP
#define LAYERED_SQT_LAYER_HPP

#include <layered-sqt/common.hpp>
#include <layered-sqt/medium.hpp>

namespace ls {

/**
 * @defgroup layer_interface Layer interface
 *
 * `<layered-sqt/layer.hpp>`
 */
/**@{*/

/**
 * @brief Layer interface.
 */
class Layer
{
public:

    /**
     * @brief Destructor.
     */
    virtual ~Layer()
    {
    }

    /**
     * @brief Z-height.
     */
    Float zheight = 0;

    /**
     * @brief Medium above.
     */
    const Medium* medium_above = nullptr;

    /**
     * @brief Medium below.
     */
    const Medium* medium_below = nullptr;

public:

    /**
     * @name Interface
     */
    /**@{*/

    /**
     * @brief Load.
     */
    virtual
    void load(const std::string& strln) = 0;

    /**
     * @brief BSDF.
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
    Float bsdf(
            Pcg32& pcg,
            const Vec3<Float>& wo, 
            const Vec3<Float>& wi) const = 0;

    /**
     * @brief BSDF probability density function.
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
    Float bsdfPdf(
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>& wi) const = 0;

    /**
     * @brief BSDF probability density function sample.
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
    Vec3<Float> bsdfPdfSample(
            Pcg32& pcg, 
            Float& tau,
            const Vec3<Float>& wo) const = 0;

    /**
     * @brief Is transmissive?
     */
    virtual 
    bool isTransmissive() const = 0;

    /**@}*/

public:

    /**
     * @brief Intersect.
     *
     * @param[in] ray
     * Ray.
     *
     * @param[out] hit
     * Hit.
     */
    bool intersect(const Ray& ray, Hit& hit) const;
};

/**@}*/

} // namespace ls

#endif // #ifndef LAYERED_SQT_LAYER_HPP