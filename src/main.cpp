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
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <preform/random.hpp>
#include <preform/float_interval.hpp>
#include <preform/multi.hpp>
#include <preform/multi_math.hpp>
#include <preform/multi_random.hpp>
#include <preform/microsurface.hpp>
#include <preform/option_parser.hpp>

#ifndef DIELECTRIC_FM_ITERATIONS
#define DIELECTRIC_FM_ITERATIONS 8
#endif // #ifndef DIELECTRIC_FM_ITERATIONS

#ifndef MAX_BOUNCES
#define MAX_BOUNCES 32
#endif // #ifndef MAX_BOUNCES

/**
 * @brief Floating point type.
 */
typedef double Float;

/**
 * @brief 2-dimensional vector.
 */
template <typename T>
using Vec2 = pr::vec2<T>;

/**
 * @brief 3-dimensional vector.
 */
template <typename T>
using Vec3 = pr::vec3<T>;

/**
 * @brief 3-dimensional matrix.
 */
template <typename T>
using Mat3 = pr::mat3<T>;

/**
 * @brief Permuted congruential generator.
 */
typedef pr::pcg32 Pcg32;

/**
 * @brief Generate canonical random sample.
 */
Float generateCanonical(Pcg32& pcg)
{
    return pr::generate_canonical<Float>(pcg);
}

/**
 * @brief Generate canonical 2-dimensional random sample.
 */
Vec2<Float> generateCanonical2(Pcg32& pcg)
{
    return pr::generate_canonical<Float, 2>(pcg);
}

/**
 * @brief Generate canonical 3-dimensional random sample.
 */
Vec3<Float> generateCanonical3(Pcg32& pcg)
{
    return pr::generate_canonical<Float, 3>(pcg);
}

#if !DOXYGEN

    // Prototypes.
    class Medium;
    class Layer;

#endif // #if !DOXYGEN

/**
 * @brief Ray.
 */
class Ray
{
public:

    /**
     * @brief Position.
     */
    Vec3<Float> pos;

    /**
     * @brief Direction.
     */
    Vec3<Float> dir;

    /**
     * @brief Associated medium.
     */
    const Medium* medium = nullptr;
};

/**
 * @brief Hit.
 */
class Hit
{
public:

    /**
     * @brief Position.
     */
    Vec3<Float> pos;
};

/**
 * @brief Medium.
 *
 * TODO
 */
class Medium
{
public:

    /**
     * @brief Default constructor.
     */
    Medium() = default;

    /**
     * @brief Constructor.
     *
     * @param[in] g
     * Henyey-Greenstein parameter @f$ g @f$.
     *
     * @param[in] mua
     * Volume absorption coefficient @f$ \mu_a @f$.
     *
     * @param[in] mus
     * Volume scattering coefficient @f$ \mu_s @f$.
     *
     * @param[in] eta
     * Refractive index @f$ \eta @f$.
     *
     * @throw std::invalid_argument
     * If
     * - `g` is outside `(-1, 1)`,
     * - `mua` is less than `0`,
     * - `mus` is less than `0`, or
     * - `eta` is less than `1`.
     */
    Medium(Float g, Float mua, Float mus, Float eta) :
            g(g),
            mua(mua),
            mus(mus),
            mu(mua + mus),
            eta(eta)
    {
        const char* error_message = nullptr;

        if (!(g > -1 && g < 1)) {
            error_message = ": Henyey-Greenstein parameter is outside (-1, 1)";
        }
        else 
        if (!(mua >= 0)) {
            error_message = ": volume absorption coefficient is negative";
        }
        else 
        if (!(mus >= 0)) {
            error_message = ": volume scattering coefficient is negative";
        }
        else
        if (!(eta >= 1)) {
            error_message = ": refractive index is less than 1";
        }

        if (error_message) {
            throw 
                std::invalid_argument(
                std::string(__PRETTY_FUNCTION__).append(error_message));
        }
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
     * @brief Henyey-Greenstein parameter @f$ g \in (-1, 1) @f$.
     */
    Float g = 0;

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

    /**
     * @brief Refractive index @f$ \eta @f$.
     */
    Float eta = 1;

    /**
     * @brief Transmittance.
     *
     * @param[in] d
     * Distance.
     */
    Float transmittance(Float d) const
    {
        assert(d >= 0);
        if (d == 0 || mu == 0) {
            return 1;
        }
        else {
            return pr::exp(-mu * 
                   pr::min(d, pr::numeric_limits<Float>::max()));
        }
    }

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
            Float dmax) const
    {
        if (mu != 0) {

            // Sample distance.
            Float u = generateCanonical(pcg);
            Float d = (u < Float(0.5) ? -pr::log1p(-u) : -pr::log(1 - u)) / mu;

            // Scattering event occurred?
            if (d < dmax) {
                // Update throughput.
                tau *= mus / mu;

                // Return distance of scattering event.
                return d;
            }
            else {
                // Indicate no scattering event.
                return dmax;
            }
        }
        else {
            // Indicate no scattering event.
            return dmax;
        }
    }

    /**
     * @brief Phase function.
     *
     * @param[in] wo
     * Outgoing direction.
     *
     * @param[in] wi
     * Incident direction.
     */
    Float phase(
            const Vec3<Float>& wo,
            const Vec3<Float>& wi) const
    {
        return Vec3<Float>::hg_phase_pdf(g, pr::dot(wo, wi));
    }

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
    Vec3<Float> phaseSample(
            Pcg32& pcg,
            Float& tau,
            const Vec3<Float>& wo) const
    {
        (void) tau; // tau *= 1
        return pr::dot(
                    Mat3<Float>::build_onb(wo),
                    Vec3<Float>::hg_phase_pdf_sample(g, 
                                    generateCanonical2(pcg)));
    }
};

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
     * @brief BSDF sample.
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
    Vec3<Float> bsdfSample(
            Pcg32& pcg, 
            Float& tau,
            const Vec3<Float>& wo) const = 0;

    /**
     * @brief Intersect.
     *
     * @param[in] ray
     * Ray.
     *
     * @param[out] hit
     * Hit.
     */
    bool intersect(const Ray& ray, Hit& hit) const
    {
        Float t = (zheight - ray.pos[2]) / ray.dir[2];
        if (!(t > 0 &&
              t < pr::numeric_limits<Float>::infinity())) {
            return false;
        }

        hit.pos = {
            ray.pos[0] + ray.dir[0] * t,
            ray.pos[1] + ray.dir[1] * t,
            zheight
        };
        return true;
    }
};

/**
 * @brief Lambertian layer.
 */
class LambertianLayer final : public Layer
{
public:

    /**
     * @brief Default constructor.
     */
    LambertianLayer() = default;

    /**
     * @brief Constructor.
     *
     * @param[in] fR
     * BRDF coefficient @f$ f_R @f$.
     *
     * @param[in] fT
     * BTDF coefficient @f$ f_T @f$.
     *
     * @throw std::invalid_argument
     * If
     * - `fR` is outside `[0, 1]`,
     * - `fT` is outside `[0, 1]`, or
     * - `fR` plus `fT` is greater than `1`.
     */
    LambertianLayer(Float fR, Float fT) : 
            fR(fR), 
            fT(fT)
    {
        const char* error_message = nullptr;

        if (!(fR >= 0 && fR <= 1)) {
            error_message = ": BRDF coefficient is outside [0, 1]";
        }
        else 
        if (!(fT >= 0 && fT <= 1)) {
            error_message = ": BTDF coefficient is outside [0, 1]";
        }
        else 
        if (!(fR + fT <= 1)) {
            error_message = ": BRDF/BTDF coefficient sum is greater than 1";
        }

        if (error_message) {
            throw 
                std::invalid_argument(
                std::string(__PRETTY_FUNCTION__).append(error_message));
        }
    }

    /**
     * @brief BRDF coefficient @f$ f_R @f$.
     */
    Float fR = 0.5;

    /**
     * @brief BTDF coefficient @f$ f_T @f$.
     */
    Float fT = 0.5;

    /**
     * @copydoc Layer::bsdf()
     */
    Float bsdf(
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>& wi) const
    {
        (void) pcg;
        return pr::numeric_constants<Float>::M_1_pi() * pr::abs(wi[2]) * 
              (pr::signbit(wo[2]) == pr::signbit(wi[2]) ? fR : fT);
    }

    /**
     * @copydoc Layer::bsdfSample()
     */
    Vec3<Float> bsdfSample(
            Pcg32& pcg, 
            Float& tau,
            const Vec3<Float>& wo) const
    {
        tau *= fR + fT;
        Vec3<Float> wi = 
        Vec3<Float>::cosine_hemisphere_pdf_sample(generateCanonical2(pcg));
        if (generateCanonical(pcg) < fR / (fR + fT)) {
            wi[2] = pr::copysign(wi[2], +wo[2]);
        }
        else {
            wi[2] = pr::copysign(wi[2], -wo[2]);
        }
        return wi;
    }
};

/**
 * @brief Dielectric layer.
 */
class DielectricLayer final : public Layer
{
public:

    /**
     * @brief Constructor.
     */
    DielectricLayer() = default;

    /**
     * @brief Constructor.
     *
     * @param[in] kR
     * BRDF coefficient @f$ k_R @f$.
     *
     * @param[in] kT
     * BTDF coefficient @f$ k_T @f$.
     *
     * @param[in] alpha
     * Roughness @f$ \alpha @f$.
     *
     * @throw std::invalid_argument
     * If
     * - `kR` is outside `[0, 1]`,
     * - `kT` is outside `[0, 1]`, or
     * - `alpha` is non-positive.
     */
    DielectricLayer(Float kR, Float kT, Float alpha) :
            kR(kR), 
            kT(kT),
            alpha(alpha)
    {
        const char* error_message = nullptr;

        if (!(kR >= 0 && kR <= 1)) {
            error_message = ": BRDF coefficient is outside [0, 1]";
        }
        else 
        if (!(kT >= 0 && kT <= 1)) {
            error_message = ": BTDF coefficient is outside [0, 1]";
        }
        else 
        if (!(alpha > 0)) {
            error_message = ": roughness is non-positive";
        }

        if (error_message) {
            throw 
                std::invalid_argument(
                std::string(__PRETTY_FUNCTION__).append(error_message));
        }
    }

    /**
     * @brief BRDF coefficient @f$ k_R @f$.
     */
    Float kR = 1;

    /**
     * @brief BTDF coefficient @f$ k_T @f$.
     */
    Float kT = 1;

    /**
     * @brief Roughness @f$ \alpha @f$.
     */
    Float alpha = 0.5;

    /**
     * @copydoc Layer::bsdf()
     */
    Float bsdf(
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>& wi) const
    {
        pr::microsurface_dielectric_bsdf<
            Float,
            pr::microsurface_trowbridge_reitz_slope,
            pr::microsurface_uniform_height> model(
                    kR, kT, 
                    this->medium_above->eta /
                    this->medium_below->eta,
                    Vec2<Float>{alpha, alpha});

        return 
            model.fm(
            [&pcg]() -> Float {
                return generateCanonical(pcg);
            },
            wo, wi, 0,
            DIELECTRIC_FM_ITERATIONS);
    }

    /**
     * @copydoc Layer::bsdfSample()
     */
    Vec3<Float> bsdfSample(
            Pcg32& pcg, 
            Float& tau,
            const Vec3<Float>& wo) const
    {
        pr::microsurface_dielectric_bsdf<
            Float,
            pr::microsurface_trowbridge_reitz_slope,
            pr::microsurface_uniform_height> model(
                    kR, kT, 
                    this->medium_above->eta /
                    this->medium_below->eta,
                    Vec2<Float>{alpha, alpha});

        int k; 
        Vec3<Float> wi = 
            model.fm_pdf_sample(
            [&pcg]() -> Float {
                return generateCanonical(pcg);
            },
            wo, k);
        if (!(kR == 1 &&
              kT == 1)) {
            Float fm_pdf;
            Float fm = 
                model.fm(
                [&pcg]() -> Float {
                    return generateCanonical(pcg);
                },
                wo, wi, 0, 
                DIELECTRIC_FM_ITERATIONS,
                &fm_pdf);
            tau *= fm / fm_pdf;
        }
        return wi;
    }
};

/**
 * @brief Assembly.
 *
 * TODO
 */
class Assembly
{
public:

    /**
     * @brief Destructor.
     */
    ~Assembly()
    {
        clear();
    }

    /**
     * @brief Mediums.
     */
    std::vector<Medium*> mediums;

    /**
     * @brief Layers.
     */
    std::vector<Layer*> layers;

    /**
     * @brief Load from input stream.
     *
     * TODO
     */
    void load(std::istream& istr)
    {
        clear();

        bool is_medium = true;
        char line[128];
        long lineno = 1;

        while (istr.getline(line, sizeof(line))) {
            std::stringstream isstr(line);
            std::string str;
            isstr.setf(std::ios_base::skipws);
            isstr >> str;

            // Expecting medium?
            if (is_medium) {

                // Medium arguments.
                Float g = 0;
                Float mua = 0;
                Float mus = 0;
                Float eta = 1;

                // Keyword 'Vacuum'?
                if (str == "Vacuum") {
                    // Use defaults.
                }
                // Keyword 'Medium'?
                else if (str == "Medium") {
                    try {
                        // Read arguments.
                        // Note std::stod() throws if string is invalid.
                        while (isstr >> str) {
                            if (!str.compare(0, 2, "g=", 2)) {
                                g = std::stod(str.substr(2));
                            }
                            else 
                            if (!str.compare(0, 4, "mua=", 4)) {
                                mua = std::stod(str.substr(4));
                            }
                            else 
                            if (!str.compare(0, 4, "mus=", 4)) {
                                mus = std::stod(str.substr(4));
                            }
                            else 
                            if (!str.compare(0, 4, "eta=", 4)) {
                                eta = std::stod(str.substr(4));
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
                                .append(": on line ")
                                .append(std::to_string(lineno))
                                .append(", invalid argument '").append(str)
                                .append("'"));
                    }
                }
                else {
                    // Runtime error,
                    // keyword neither 'Vacuum' nor 'Medium'.
                    throw 
                        std::runtime_error(
                        std::string(__PRETTY_FUNCTION__)
                            .append(": on line ")
                            .append(std::to_string(lineno))
                            .append(", expected 'Vacuum' or 'Medium'"));
                }

                try {
                    // Push medium.
                    mediums.emplace_back(
                            new Medium(g, mua, mus, eta));
                }
                catch (const std::invalid_argument& except) {
                    throw 
                        std::runtime_error(
                        std::string(__PRETTY_FUNCTION__)
                            .append(": on line ")
                            .append(std::to_string(lineno))
                            .append(", caught exception: ")
                            .append(except.what()));
                }
            }
            else {
                if (str == "Layer") {
                    Float zheight = 0;
                    // No z=argument?
                    if (!(isstr >> str) ||
                        str.compare(0, 2, "z=", 2) != 0) {
                        // Runtime error.
                        throw
                            std::runtime_error(
                            std::string(__PRETTY_FUNCTION__)
                                .append(": on line ")
                                .append(std::to_string(lineno))
                                .append(", expected z="));
                    }
                    try {
                        // Read zheight.
                        // Note std::stod() throws if string is invalid.
                        zheight = std::stod(str.substr(2));
                    }
                    catch (...) {
                        // Runtime error.
                        throw
                            std::runtime_error(
                            std::string(__PRETTY_FUNCTION__)
                                .append(": on line ")
                                .append(std::to_string(lineno))
                                .append(", invalid argument '").append(str)
                                .append("'"));
                    }

                    if (layers.size() > 0 &&
                        !(layers.back()->zheight > zheight)) {
                        // Runtime error.
                        throw
                            std::runtime_error(
                            std::string(__PRETTY_FUNCTION__)
                                .append(": on line ")
                                .append(std::to_string(lineno))
                                .append(", layer z-height is out-of-order"));
                    }

                    // Keyword.
                    isstr >> str;

                    // Keyword 'Lambertian'?
                    if (str == "Lambertian") {

                        // Lambertian layer arguments.
                        Float fR = 0.5;
                        Float fT = 0.5;

                        try {
                            // Read arguments.
                            // Note std::stod() throws if string is invalid.
                            while (isstr >> str) {
                                if (str.compare(0, 3, "fR=", 3) == 0) {
                                    fR = std::stod(str.substr(3));
                                }
                                else 
                                if (str.compare(0, 3, "fT=", 3) == 0) {
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
                                    .append(": on line ")
                                    .append(std::to_string(lineno))
                                    .append(", invalid argument '").append(str)
                                    .append("'"));
                        }

                        try {
                            // Push layer.
                            layers.emplace_back(
                                    new LambertianLayer(fR, fT));
                        }
                        catch (const std::invalid_argument& except) {
                            throw 
                                std::runtime_error(
                                std::string(__PRETTY_FUNCTION__)
                                    .append(": on line ")
                                    .append(std::to_string(lineno))
                                    .append(", caught exception: ")
                                    .append(except.what()));
                        }
                    }
                    // Keyword 'Dielectric'?
                    else if (str == "Dielectric") {

                        // Dielectric layer arguments.
                        Float kR = 1;
                        Float kT = 1;
                        Float alpha = 0.5;

                        try {
                            // Read arguments.
                            // Note std::stod() throws if string is invalid.
                            while (isstr >> str) {
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
                                    .append(": on line ")
                                    .append(std::to_string(lineno))
                                    .append(", invalid argument '").append(str)
                                    .append("'"));
                        }

                        try {
                            // Push layer.
                            layers.emplace_back(
                                    new DielectricLayer(kR, kT, alpha));
                        }
                        catch (const std::invalid_argument& except) {
                            throw 
                                std::runtime_error(
                                std::string(__PRETTY_FUNCTION__)
                                    .append(": on line ")
                                    .append(std::to_string(lineno))
                                    .append(", caught exception: ")
                                    .append(except.what()));
                        }
                    }
                    else {
                        // Runtime error,
                        // keyword neither 'Lambertian' nor 'Dielectric'.
                        throw 
                            std::runtime_error(
                            std::string(__PRETTY_FUNCTION__)
                                .append(": on line ")
                                .append(std::to_string(lineno))
                                .append(", expected "
                                        "'Lambertian' or 'Dielectric' "
                                        "description"));
                    }

                    // Set zheight.
                    layers.back()->zheight = zheight;
                }
                else {
                    // Runtime error,
                    // keyword not 'Layer'.
                    throw 
                        std::runtime_error(
                        std::string(__PRETTY_FUNCTION__)
                            .append(": on line ")
                            .append(std::to_string(lineno))
                            .append(", expected 'Layer'"));
                }
            }

            // Flip.
            is_medium = !is_medium;

            // Increment line number.
            lineno++;
        }

        // Expecting medium?
        if (is_medium) {
            throw 
                std::runtime_error(
                std::string(__PRETTY_FUNCTION__)
                    .append(": missing bottom medium"));

        }

        // No layers?
        if (layers.empty()) {
            throw 
                std::runtime_error(
                std::string(__PRETTY_FUNCTION__)
                    .append(": no layers!"));
        }

        for (std::size_t pos = 0;
                         pos < layers.size(); pos++) {

            // Link pointers.
            mediums[pos + 0]->layer_below = layers[pos];
            mediums[pos + 1]->layer_above = layers[pos];
            layers[pos]->medium_above = mediums[pos + 0];
            layers[pos]->medium_below = mediums[pos + 1];
        }
    }

    /**
     * @brief Clear.
     */
    void clear()
    {
        // Delete mediums.
        for (auto medium : mediums) {
            delete medium;
        }
        mediums.clear();
        mediums.shrink_to_fit();

        // Delete layers.
        for (auto layer : layers) {
            delete layer;
        }
        layers.clear();
        layers.shrink_to_fit();
    }

    /**
     * @brief Value.
     * 
     * @param[inout] pcg
     * Generator.
     *
     * @param[in] wo0
     * Outgoing direction @f$ \omega_{o,0} @f$.
     *
     * @param[in] wi0
     * Incident direction @f$ \omega_{i,0} @f$.
     */
    Float value(
            Pcg32& pcg,
            const Vec3<Float>& wo0,
            const Vec3<Float>& wi0) 
    {
        // Initial layer.
        const Layer* layer;

        // Initial ray.
        Ray ray;
        ray.dir = -wo0;

        // Upper hemisphere?
        if (wo0[2] > 0) {
            layer = layers.front();
            // Move above top layer.
            ray.pos[2] = layer->zheight + 1;
            ray.medium = layer->medium_above;
        }
        else {
            layer = layers.back();
            // Move below bottom layer.
            ray.pos[2] = layer->zheight - 1;
            ray.medium = layer->medium_below;
        }

        Float f = 0;
        Float tau = 1;
        for (int bounce = 0;
                 bounce < MAX_BOUNCES; bounce++) {

            Vec3<Float> wo = -ray.dir;

            // Determine neighboring layers.
            const Layer* layer_above;
            const Layer* layer_below;
            // Ray exactly on layer?
            if (ray.pos[2] == layer->zheight) {
                layer_above = layer->medium_above->layer_above;
                layer_below = layer->medium_below->layer_below;
            }
            else {
                layer_above = ray.medium->layer_above;
                layer_below = ray.medium->layer_below;
            }
            assert(!layer_above || ray.pos[2] < layer_above->zheight);
            assert(!layer_below || ray.pos[2] > layer_below->zheight);

            // Intersect relevant layer.
            Hit hit;
            const Layer* layer_next = ray.dir[2] > 0 
                ? layer_above 
                : layer_below;
            if (!layer_next ||
                !layer_next->intersect(ray, hit)) {
                break;
            }
            layer = layer_next;

            // Sample medium.
            Float dmax = pr::length(hit.pos - ray.pos);
            Float d = ray.medium->transmittanceSample(pcg, tau, dmax);
            if (!(d == dmax)) {
    
                // Update ray.
                ray.pos = ray.pos + ray.dir * d;
                ray.dir = ray.medium->phaseSample(pcg, tau, wo);
            }
            else {

                // If at top and wi0 is upper hemisphere OR
                // if at bottom and wi0 is lower hemisphere, add to result.
                if ((layer == layers.front() && wi0[2] > 0) ||
                    (layer == layers.back()  && wi0[2] < 0)) {
                    f += tau * layer->bsdf(pcg, wo, wi0);
                }

                // Update ray.
                ray.pos = hit.pos;
                ray.dir = layer->bsdfSample(pcg, tau, wo);
                ray.medium = 
                    ray.dir[2] > 0 
                    ? layer->medium_above 
                    : layer->medium_below;
            }

            if (!(tau > 0)) {
                break;
            }
        }
        return f;
    }
};

int main(int arg, char** argv)
{
    // TODO Option parser
#if 0
    Assembly assembly;
    std::ifstream ifs("example.layers");
    assembly.load(ifs);
    ifs.close();
#endif
    return 0;
}
