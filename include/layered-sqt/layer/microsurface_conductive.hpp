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
#ifndef LAYERED_SQT_LAYER_MICROSURFACE_CONDUCTIVE_HPP
#define LAYERED_SQT_LAYER_MICROSURFACE_CONDUCTIVE_HPP

#include <layered-sqt/layer.hpp>

namespace ls {

/**
 * @addtogroup layers Layers
 *
 * `<layered-sqt/layer/microsurface_conductive.hpp>`
 */
/**@{*/

/**
 * @brief Microsurface conductive layer.
 */
class MicrosurfaceConductiveLayer final : public Layer
{
public:

    /**
     * @brief Default constructor.
     */
    MicrosurfaceConductiveLayer() = default;

    /**
     * @brief Roughness @f$ \alpha @f$.
     */
    Float alpha = 0.5;

    /**
     * @brief Use multiple scattering?
     */
    bool use_multiple_scattering = false;

public:

    /**
     * @name Interface
     */
    /**@{*/

    /**
     * @copydoc Layer::init()
     *
     * Accepts arguments
     * - `alpha=`(float),
     * - `use_multiple_scattering=`(bool).
     *
     * @throw std::runtime_error
     * If
     * - invalid argument, or
     * - `alpha` is non-positive.
     */
    void init(const std::string& arg);

    /**
     * @copydoc Layer::validate()
     */
    void validate() const;

    /**
     * @copydoc Layer::bsdf()
     */
    Float bsdf(
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>& wi,
            Float* fs_pdf = nullptr) const;

    /**
     * @copydoc Layer::bsdfSample()
     */
    Vec3<Float> bsdfSample(
            Pcg32& pcg, 
            Float& tau,
            const Vec3<Float>& wo) const;

    /**
     * @copydoc Layer::isTransmissive()
     */
    bool isTransmissive() const;

    /**@}*/
};

/**@}*/

} // namespace ls

#endif // #ifndef LAYERED_SQT_LAYER_MICROSURFACE_CONDUCTIVE_HPP
