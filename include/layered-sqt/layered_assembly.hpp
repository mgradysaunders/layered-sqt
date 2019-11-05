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
#ifndef LAYERED_SQT_LAYERED_ASSEMBLY_HPP
#define LAYERED_SQT_LAYERED_ASSEMBLY_HPP

#include <layered-sqt/layer.hpp>

namespace ls {

/**
 * @defgroup layered_assembly Layered assembly
 *
 * `<layered-sqt/layered_assembly.hpp>`
 */
/**@{*/

/**
 * @brief Layered assembly.
 */
class LayeredAssembly
{
public:

    /**
     * @brief Default constructor.
     */
    LayeredAssembly() = default;

    /**
     * @brief Destructor.
     */
    ~LayeredAssembly()
    {
        clear();
    }

    /**
     * @brief Initialize from input stream.
     *
     * @param[inout] is
     * Input stream.
     */
    void init(std::istream& is);

    /**
     * @brief Clear.
     */
    void clear();

    /**
     * @brief Compute BSDF/BSDF-PDF.
     *
     * @param[inout] pcg
     * Generator.
     *
     * @param[in] wo
     * Outgoing direction.
     *
     * @param[in] wi
     * Incident directions.
     *
     * @param[in] wi_count
     * Incident direction count.
     *
     * @param[out] f
     * _Optional_. BSDFs per incident direction.
     *
     * @param[out] f_pdf
     * _Optional_. BSDF-PDFs per incident direction.
     */
    void compute(
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>* wi, int wi_count,
            Float* f,
            Float* f_pdf) const;

    /**
     * @brief Compute BSDF/BSDF-PDF average.
     *
     * @param[in] path_count_before
     * Path count in average before this call.
     *
     * @param[in] path_count
     * Path count in average after this call.
     *
     * @param[inout] pcg
     * Generator.
     *
     * @param[in] wo
     * Outgoing direction.
     *
     * @param[in] wi
     * Incident directions.
     *
     * @param[in] wi_count
     * Incident direction count.
     *
     * @param[inout] f
     * _Optional_. BSDFs per incident direction.
     *
     * @param[inout] f_pdf
     * _Optional_. BSDF-PDFs per incident direction.
     */
    void computeAverage(
            int path_count_before,
            int path_count,
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>* wi, int wi_count,
            Float* f,
            Float* f_pdf) const;

    /**
     * @brief Random scatter direction.
     *
     * @param[inout] pcg
     * Generator.
     *
     * @param[in] wo
     * Outgoing direction.
     */
    Vec3<Float> randomScatterDirection(
            Pcg32& pcg, 
            const Vec3<Float>& wo) const;

    /**
     * @brief Is transmissive?
     *
     * If every layer is transmissive, then the assembly is transmissive. 
     * Alternatively, if even a single layer is not transmissive, then the
     * assembly is not transmissive.
     */
    bool isTransmissive() const
    {
        for (const Layer* layer : layers_) {
            if (!layer->isTransmissive()) {
                return false;
            }
        }
        return true;
    }

private:

    /**
     * @brief Mediums (media).
     */
    std::vector<Medium*> mediums_;

    /**
     * @brief Layers.
     */
    std::vector<Layer*> layers_;

    /**
     * @brief Non-null BSDF layer at top.
     */
    const Layer* layers_top_ = nullptr;

    /**
     * @brief Non-null BSDF layer at bottom.
     */
    const Layer* layers_bottom_ = nullptr;
};

/**@}*/

} // namespace ls

#endif // #ifndef LAYERED_SQT_LAYERED_ASSEMBLY_HPP
