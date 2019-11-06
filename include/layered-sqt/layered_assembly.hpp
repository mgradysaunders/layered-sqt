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
 *
 * This characterizes a _layered assembly_, consisting of `Layer`s 
 * separated by `Medium`s, and allows for Monte Carlo evaluation
 * of the emergent BSDF/BSDF-PDF.
 *
 * @par Layered assembly structure
 *
 * A layered assembly must contain at least 1 layer and at least 2
 * participating media. In general, a layered assembly may contain @f$ N @f$
 * layers and @f$ N + 1 @f$ participating media. In essence, a layer
 * characterizes a surface, while a medium characterizes the space between
 * surfaces (or, in the boundary cases, the spaces above and below the top 
 * and bottom surfaces respectively).
 *
 * Each `Layer` stores its Z-height in the member variable `zheight` and 
 * implements a BSDF model via virtual functions. The implementation requires 
 * that Z-heights of layers in an assembly be strictly decreasing, so 
 * that layers appear in top-to-bottom order. `NullBsdfLayer` is a special 
 * sub-class representing a layer without a BSDF, useful to separate media
 * _without_ surface scattering at a layer in between. The implementation 
 * allows `NullBsdfLayer`s to appear anywhere in the assembly. Moreover, an
 * assembly may consist entirely of `NullBsdfLayers` to represent scattering
 * in a layered cloud of sorts&mdash;beware the implicit delta function,
 * this is inadvisable for low scattering coefficients!
 *
 * Each `Medium` stores its refractive index @f$ \eta @f$ in the member
 * variable `eta`, its Henyey-Greenstein parameter @f$ g @f$ in the member
 * variable `g`, and its absorption, scattering, and extinction coefficients 
 * @f$ \mu_a @f$, @f$ \mu_s @f$, and @f$ \mu = \mu_a + \mu_s @f$ in the 
 * member variables `mua`, `mus`, and `mu` respectively. The implementation 
 * considers that &ldquo;vacuum&rdquo; is a medium with @f$ \eta = 1 @f$ 
 * and @f$ \mu_a = \mu_s = \mu = 0 @f$, and _not a null pointer_. The 
 * implementation requires that the boundary media be non-absorbing and 
 * non-scattering such that @f$ \mu_a = \mu_s = \mu = 0 @f$, though 
 * @f$ \eta @f$ need not be unity. 
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
     * Initialize from input stream, which is by assumption a 
     * line-by-line plain-text description of the layered assembly
     * from top to bottom. So, the first line describes the top
     * medium, the second line describes the layer beneath the
     * medium, the third line describes the medium beneath the
     * previous layer, and so on until the bottom medium.
     *
     * A line describing a medium must begin with the keyword `Medium`
     * and may be followed by optional keyword arguments with syntax
     * `key=val` (no spaces!) where `key` is either `eta`, `g`, `mus`, or 
     * `mua`. By default, `eta=1`, `g=0`, `mus=0`, and `mua=0`. Hence, to
     * specify a vacuum medium,
     * ~~~~~~
     * Medium
     * ~~~~~~
     * is sufficient. As a slightly more interesting example, 
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Medium eta=1.5 mua=0 mus=0.8 g=-0.2 
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * specifies a medium with the refractive index of glass, no volume
     * absorption, and a bit of volume back-scattering. It may be worth
     * noting that `mua=0` in the above line could be left out since 
     * `mua` is `0` by default.
     *
     * A line describing a layer must begin with the keyword `Layer` and
     * be followed by a required keyword argument `z=val` to specify 
     * Z-height. This must be followed by the keyword for a particular 
     * sub-class (either `NullBsdf`, `LambertianBsdf`,
     * `MicrosurfaceLambertianBrdf`, or `MicrosurfaceDielectricBsdf`), 
     * which in turn may be followed by optional keyword arguments with 
     * syntax `key=val` (no spaces!) associated with the sub-class. 
     * For example,
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Layer z=3.2 LambertianBsdf fR=0.7 fT=0.2
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * specifies a layer at Z-height 3.2 with a Lambertian BSDF that is
     * 70% reflective, 20% transmissive, and (per energy conservation)
     * 10% absorptive.
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
