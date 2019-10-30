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
#include <preform/multi.hpp>
#include <preform/multi_math.hpp>
#include <preform/multi_random.hpp>
#include <preform/microsurface.hpp>
#include <preform/kdtree.hpp>
#include <preform/thread_pool.hpp>
#include <preform/option_parser.hpp>
#include <preform/bash_format.hpp>

#ifndef DIELECTRIC_FM_ITERATIONS
#define DIELECTRIC_FM_ITERATIONS 8
#endif // #ifndef DIELECTRIC_FM_ITERATIONS

#ifndef MAX_BOUNCES
#define MAX_BOUNCES 64
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

/**
 * @brief Redundancy-reduced sample set.
 */
struct Rrss
{
public:

    /**
     * @brief Sample.
     */
    struct Sample
    {
    public:

        /**
         * @brief Direction.
         */
        Vec3<Float> dir;
        
        /**
         * @brief Direction probability density function.
         */
        Float dir_pdf = 0;

        /**
         * @brief Redundancy.
         */
        mutable Float redundancy = 1;

        /**
         * @brief Is enabled?
         */
        mutable bool is_enabled = true;

    public:

        /**
         * @brief Default constructor.
         */
        Sample() = default;

        /**
         * @brief Constructor.
         *
         * @param[in] dir
         * Direction.
         *
         * @param[in] dir_pdf
         * Direction probability density function.
         */
        Sample(const Vec3<Float>& dir, Float dir_pdf) :
            dir(dir),
            dir_pdf(dir_pdf)
        {
        }
    };

public:

    /**
     * @brief Constructor.
     */
    Rrss(const std::vector<Sample>& samples) : 
            samples_(samples),
            samples_remaining_(samples.size())
    {
        for (Sample& sample : samples_) {
            sample.redundancy = 1;
            sample.is_enabled = true;
        }

        sample_tree_.init(
        samples_.begin(), 
        samples_.end(),
        [&](const Sample& sample) -> 
            std::pair<Vec3<Float>, std::size_t> {
            return {
                sample.dir,
                std::size_t(&sample - &samples_[0])
            };
        });
    }

    /**
     * @brief Samples.
     */
    const std::vector<Sample>& samples() const
    {
        return samples_;
    }

public:

    /**
     * @brief Disable most redundant sample.
     */
    bool disableMostRedundant()
    {
        Sample* most_redundant = nullptr;
        for (Sample& sample : samples_) {

            // Is enabled?
            if (sample.is_enabled) {

                int num_enabled = 0;
                int num = 0;
                Float dir_pdf_sum = 0;

                // Visit samples in region.
                sample_tree_.nearby(
                sample.dir, 
                pr::sqrt(
                pr::numeric_constants<Float>::M_1_pi() /
                (sample.dir_pdf * samples_remaining_)),
                [&](auto node) {
                    Sample& node_sample = 
                    samples_[node->value.second];
                    if (node_sample.is_enabled) {
                        num_enabled++;
                    }
                    num++;
                    dir_pdf_sum += node_sample.dir_pdf;
                    return true;
                });

                // Update redundancy.
                sample.redundancy = 
                    num_enabled * 
                    num * (sample.dir_pdf / dir_pdf_sum);

                // Update most redundant sample.
                if (most_redundant == nullptr ||
                    most_redundant->redundancy < sample.redundancy) {
                    most_redundant = &sample;
                }
            }
        }

        if (most_redundant) {
            most_redundant->is_enabled = false;
            samples_remaining_--;
            return true;
        }
        else {
            return false;
        }
    }

    /**
     * @brief Disable most redundant samples.
     *
     * @param[in] num
     * Number of samples to disable.
     */
    bool disableMostRedundant(int num)
    {
        while (num-- > 0) {
            if (!disableMostRedundant()) {
                return false;
            }
        }
        return true;
    }

private:


    /**
     * @brief Samples.
     */
    std::vector<Sample> samples_;

    /**
     * @brief Samples remaining (i.e., still enabled).
     */
    std::size_t samples_remaining_ = 0;

    /**
     * @brief Sample tree.
     */
    pr::kdtree<Float, 3, std::size_t> sample_tree_;
};

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
            Float dmax) const // Float* dpdf = nullptr
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
    virtual bool isTransmissive() const = 0;

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
            zheight // Set exact.
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
     * @copydoc Layer::bsdfPdf()
     */
    Float bsdfPdf(
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>& wi) const
    {
        (void) pcg;
        return Vec3<Float>::cosine_hemisphere_pdf(pr::abs(wi[2])) *
              (pr::signbit(wo[2]) == pr::signbit(wi[2]) ? fR : fT) / (fR + fT);
    }

    /**
     * @copydoc Layer::bsdfPdfSample()
     */
    Vec3<Float> bsdfPdfSample(
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

    /**
     * @copydoc Layer::isTransmissive()
     */
    bool isTransmissive() const
    {
        return fT > 0;
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
            pr::microsurface_uniform_height> surf(
                    kR, kT, 
                    this->medium_above->eta /
                    this->medium_below->eta,
                    Vec2<Float>{alpha, alpha});

        return 
            surf.fm(
            [&pcg]() -> Float {
                return generateCanonical(pcg);
            },
            wo, wi, 0,
            DIELECTRIC_FM_ITERATIONS);
    }

    /**
     * @copydoc Layer::bsdfPdf()
     */
    Float bsdfPdf(
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>& wi) const
    {
        pr::microsurface_dielectric_bsdf<
            Float,
            pr::microsurface_trowbridge_reitz_slope,
            pr::microsurface_uniform_height> surf(
                    kR, kT, 
                    this->medium_above->eta /
                    this->medium_below->eta,
                    Vec2<Float>{alpha, alpha});

        return 
            surf.fm_pdf(
            [&pcg]() -> Float {
                return generateCanonical(pcg);
            },
            wo, wi, 0,
            DIELECTRIC_FM_ITERATIONS);
    }

    /**
     * @copydoc Layer::bsdfPdfSample()
     */
    Vec3<Float> bsdfPdfSample(
            Pcg32& pcg, 
            Float& tau,
            const Vec3<Float>& wo) const
    {
        pr::microsurface_dielectric_bsdf<
            Float,
            pr::microsurface_trowbridge_reitz_slope,
            pr::microsurface_uniform_height> surf(
                    kR, kT, 
                    this->medium_above->eta /
                    this->medium_below->eta,
                    Vec2<Float>{alpha, alpha});

        int k; 
        Vec3<Float> wi = 
            surf.fm_pdf_sample(
            [&pcg]() -> Float {
                return generateCanonical(pcg);
            },
            wo, k);
        if (!(kR == 1 &&
              kT == 1)) {
            Float fm_pdf;
            Float fm = 
                surf.fm(
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

    /**
     * @copydoc Layer::isTransmissive()
     */
    bool isTransmissive() const
    {
        return kT > 0;
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
    Float bsdf(
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>& wi) const 
    {
        // Initial layer.
        const Layer* layer;

        // Initial ray.
        Ray ray;
        ray.dir = -wo;

        // Upper hemisphere?
        if (wo[2] > 0) {
            // Move above top layer.
            layer = layers.front();
            ray.pos[2] = layer->zheight + 1;
            ray.medium = layer->medium_above;
        }
        else {
            // Move below bottom layer.
            layer = layers.back();
            ray.pos[2] = layer->zheight - 1;
            ray.medium = layer->medium_below;
        }

        Float f = 0;
        Float tau = 1;
        for (int bounce = 0;
                 bounce < MAX_BOUNCES; bounce++) {

            Vec3<Float> wk = -ray.dir;

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
                ray.dir = ray.medium->phaseSample(pcg, tau, wk);
            }
            else {

                // If at top and wi is upper hemisphere OR
                // if at bottom and wi is lower hemisphere, add to result.
                if ((layer == layers.front() && wi[2] > 0) ||
                    (layer == layers.back()  && wi[2] < 0)) {
                    f += tau * layer->bsdf(pcg, wk, wi);
                /*  fpdf += layer->bsdfPdf(pcg, wk, wi);  */
                }

                // Update ray.
                ray.pos = hit.pos;
                ray.dir = layer->bsdfPdfSample(pcg, tau, wk);
                ray.medium = 
                    ray.dir[2] > 0 
                    ? layer->medium_above 
                    : layer->medium_below;
            }

            // Throughput non-positive?
            // This test passes if throughput is zero, in which case
            // path is absorbed, but also if throughput is NaN. This shouldn't
            // happen, but we want to terminate if it does.
            if (!(tau > 0)) {
                break;
            }
        }
        return f;
    }

    void bsdf(
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>* wi, int num_wi,
            Float* f) const 
    {
        // Initial layer.
        const Layer* layer;

        // Initial ray.
        Ray ray;
        ray.dir = -wo;

        // Upper hemisphere?
        if (wo[2] > 0) {
            // Move above top layer.
            layer = layers.front();
            ray.pos[2] = layer->zheight + 1;
            ray.medium = layer->medium_above;
        }
        else {
            // Move below bottom layer.
            layer = layers.back();
            ray.pos[2] = layer->zheight - 1;
            ray.medium = layer->medium_below;
        }

        for (int j = 0; j < num_wi; j++) {
            f[j] = 0;
        }

        Float tau = 1;
        for (int bounce = 0;
                 bounce < MAX_BOUNCES; bounce++) {

            Vec3<Float> wk = -ray.dir;

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
                ray.dir = ray.medium->phaseSample(pcg, tau, wk);
            }
            else {

                // If at top and wi is upper hemisphere OR
                // if at bottom and wi is lower hemisphere, add to result.
                if (layer == layers.front()) {
                    for (int j = 0; j < num_wi; j++) {
                        if (wi[j][2] > 0) {
                            f[j] += tau * layer->bsdf(pcg, wk, wi[j]);
                        }
                    }
                }
                else 
                if (layer == layers.back()) {
                    for (int j = 0; j < num_wi; j++) {
                        if (wi[j][2] < 0) {
                            f[j] += tau * layer->bsdf(pcg, wk, wi[j]);
                        }
                    }
                }

                // Update ray.
                ray.pos = hit.pos;
                ray.dir = layer->bsdfPdfSample(pcg, tau, wk);
                ray.medium = 
                    ray.dir[2] > 0 
                    ? layer->medium_above 
                    : layer->medium_below;
            }

            // Throughput non-positive?
            // This test passes if throughput is zero, in which case
            // path is absorbed, but also if throughput is NaN. This shouldn't
            // happen, but we want to terminate if it does.
            if (!(tau > 0)) {
                break;
            }
        }
    }

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
    Float bsdfPdf(
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>& wi,
            Float* f = nullptr) const
    {
        // Initial layer.
        const Layer* layer;

        // Initial ray.
        Ray ray;
        ray.dir = -wo;

        // Upper hemisphere?
        if (wo[2] > 0) {
            // Move above top layer.
            layer = layers.front();
            ray.pos[2] = layer->zheight + 1;
            ray.medium = layer->medium_above;
        }
        else {
            // Move below bottom layer.
            layer = layers.back();
            ray.pos[2] = layer->zheight - 1;
            ray.medium = layer->medium_below;
        }

        if (f) {
            *f = 0;
        }

        Float fpdf = 0;
        Float tau = 1;
        for (int bounce = 0;
                 bounce < MAX_BOUNCES; bounce++) {

            Vec3<Float> wk = -ray.dir;

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
                ray.dir = ray.medium->phaseSample(pcg, tau, wk);
            }
            else {

                // If at top and wi is upper hemisphere OR
                // if at bottom and wi is lower hemisphere, add to result.
                if ((layer == layers.front() && wi[2] > 0) ||
                    (layer == layers.back()  && wi[2] < 0)) {
                    if (f) {
                        *f += tau * layer->bsdf(pcg, wk, wi);
                    }
                    // We're exactly sampling the BSDF-PDFs we're integrating,
                    // so the throughput is implicitly 1.
                    fpdf += layer->bsdfPdf(pcg, wk, wi);
                }

                // Update ray.
                ray.pos = hit.pos;
                ray.dir = layer->bsdfPdfSample(pcg, tau, wk);
                ray.medium = 
                    ray.dir[2] > 0 
                    ? layer->medium_above 
                    : layer->medium_below;
            }

            // Throughput non-positive?
            // This test passes if throughput is zero, in which case
            // path is absorbed, but also if throughput is NaN. This shouldn't
            // happen, but we want to terminate if it does.
            if (!(tau > 0)) {
                // Terminate.
                break;
            }
        }
        return fpdf;
    }

    /**
     * @brief BSDF probability density function sample.
     *
     * @param[inout] pcg
     * Generator.
     *
     * @param[in] wo
     * Outgoing direction.
     *
     * @returns
     * Incident direction.
     */
    Vec3<Float> bsdfPdfSample(
            Pcg32& pcg,
            const Vec3<Float>& wo) const
    {
        // Initial layer.
        const Layer* layer;

        // Initial ray.
        Ray ray;
        ray.dir = -wo;

        // Upper hemisphere?
        if (wo[2] > 0) {
            // Move above top layer.
            layer = layers.front();
            ray.pos[2] = layer->zheight + 1;
            ray.medium = layer->medium_above;
        }
        else {
            // Move below bottom layer.
            layer = layers.back();
            ray.pos[2] = layer->zheight - 1;
            ray.medium = layer->medium_below;
        }

        Float tau = 1;
        for (int bounce = 0;
                 bounce < MAX_BOUNCES; bounce++) {

            Vec3<Float> wk = -ray.dir;

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
            if (!layer_next) {
                // Exit.
                return ray.dir;
            }
            if (!layer_next->intersect(ray, hit)) {
                // Terminate.
                return Vec3<Float>{};
            }
            layer = layer_next;

            // Sample medium.
            Float dmax = pr::length(hit.pos - ray.pos);
            Float d = ray.medium->transmittanceSample(pcg, tau, dmax);
            if (!(d == dmax)) {
    
                // Update ray.
                ray.pos = ray.pos + ray.dir * d;
                ray.dir = ray.medium->phaseSample(pcg, tau, wk);
            }
            else {

                // Update ray.
                ray.pos = hit.pos;
                ray.dir = layer->bsdfPdfSample(pcg, tau, wk);
                ray.medium = 
                    ray.dir[2] > 0 
                    ? layer->medium_above 
                    : layer->medium_below;
            }

            // Throughput non-positive?
            // This test passes if throughput is zero, in which case
            // path is absorbed, but also if throughput is NaN. This shouldn't
            // happen, but we want to terminate if it does.
            if (!(tau > 0)) {
                // Terminate.
                return Vec3<Float>{};
            }
        }

        // Terminate.
        return Vec3<Float>{};
    }

    /**
     * @brief Is transmissive?
     *
     * If every layer is transmissive, then the assembly is transmissive. 
     * Alternatively, if even a single layer is not transmissive, then the
     * assembly is not transmissive.
     */
    bool isTransmissive() const
    {
        for (const Layer* layer : layers) {
            if (!layer->isTransmissive()) {
                return false;
            }
        }
        return true;
    }
};

int main(int argc, char** argv)
{
    pr::option_parser opt_parser("[OPTIONS] filename");

    int seed = 0;
    int num_iters = 16384;
    int num_wo = 16;
    int num_wi = 80;
    int num_threads = 0;
    std::string ifs_fname = "";
    std::string ofs_fname = "output.raw";

    // -s/--seed
    opt_parser.on_option(
    "-s", "--seed", 1,
    [&](char** argv) {
        try {
            seed = std::stoi(argv[0]);
        }
        catch (const std::exception&) {
            throw 
                std::runtime_error(
                std::string("-s/--seed expects 1 integer")
                    .append("(can't parse ").append(argv[0])
                    .append(")"));
        }
    })
    << "Specify seed. By default, 0.\n";

    // --num-iters
    opt_parser.on_option(
    nullptr, "--num-iters", 1,
    [&](char** argv) {
        try {
            num_iters = std::stoi(argv[0]);
            if (!(num_iters > 0)) {
                throw std::exception();
            }
        }
        catch (const std::exception&) {
            throw
                std::runtime_error(
                std::string("--num-iters expects 1 positive integer ")
                    .append("(can't parse ").append(argv[0])
                    .append(")"));
        }
    })
    << "Specify number of iterations per sample. By default, 16384.\n";

    // --num-wo
    opt_parser.on_option(nullptr, "--num-wo", 1,
    [&](char** argv) {
        try {
            num_wo = std::stoi(argv[0]);
            if (!(num_wo >= 4)) {
                throw std::exception();
            }
        }
        catch (const std::exception&) {
            throw
                std::runtime_error(
                std::string("--num-wo expects 1 integer >= 4 ")
                    .append("(can't parse ").append(argv[0])
                    .append(")"));
        }
    })
    << "Specify number of outgoing directions. By default, 16.\n"
       "This is the number of outgoing directions in the upper hemisphere,\n"
       "uniformly distributed in zenith. As the emergent BRDF/BSDF must be\n"
       "isotropic, the implementation does not sample in azimuth.\n";

    // --num-wi
    opt_parser.on_option(nullptr, "--num-wi", 1,
    [&](char** argv) {
        try {
            num_wi = std::stoi(argv[0]);
            if (!(num_wi >= 16)) {
                throw std::exception();
            }
        }
        catch (const std::exception&) {
            throw
                std::runtime_error(
                std::string("--num-wi expects 1 integer >= 16 ")
                    .append("(can't parse ").append(argv[0])
                    .append(")"));
        }
    })
    << "Specify number of incident directions. By default, 80.\n";

    // -o/--output
    opt_parser.on_option("-o", "--output", 1,
    [&](char** argv) {
        ofs_fname = argv[0];
    })
    << "Specify ASCII-format RAW output filename. By default, output.raw.\n";

    // -h/--help
    opt_parser.on_option("-h", "--help", 0,
    [&](char**) {
        std::cout << opt_parser << std::endl;
        std::exit(EXIT_SUCCESS);
    })
    << "Display this help and exit.\n";

    // Positional.
    opt_parser.on_positional(
    [&](char* argv) {
        ifs_fname = argv;
    });

    try {
        // Parse args.
        opt_parser.parse(argc, argv);
    }
    catch (const std::exception& exception) {
        std::cerr << "Unhandled exception in command line arguments!\n";
        std::cerr << "exception.what(): " << exception.what() << "\n";
        std::exit(EXIT_FAILURE);
    }

    if (ifs_fname == "") {
        std::cout << opt_parser << std::endl;
        std::exit(EXIT_SUCCESS);
    }

    std::vector<Float> thetao_array(num_wo);
    std::vector<Vec3<Float>> wo_array(num_wo);
    std::vector<Vec3<Float>> wi_array(num_wo * num_wi);
    std::vector<Float> f_array(num_wo * num_wi);

    for (int wo_index = 0; 
             wo_index < num_wo; wo_index++) {
        
        // Initialize outgoing angle.
        thetao_array[wo_index] = 
            wo_index / Float(num_wo) * 
            pr::numeric_constants<Float>::M_pi_2();

        // Initialize outgoing direction.
        wo_array[wo_index] = {
            pr::sin(thetao_array[wo_index]), 0, 
            pr::cos(thetao_array[wo_index])
        };
    }

    {
        std::ifstream ifs(ifs_fname);
        Assembly assembly;
        assembly.load(ifs);
        ifs.close();

        std::cout << '\n';
        std::cout.flush();
        std::mutex cout_mutex;
        int num_complete = 0;

        auto job = [&](int wo_index) {
            Pcg32 pcg = Pcg32(seed, wo_index);
            const Vec3<Float>& wo = wo_array[wo_index];
            {
                std::vector<Rrss::Sample> wi_input;
                Float wi_ffac = 0;
                wi_input.reserve(4 * num_wi);
                for (int wi_index = 0; 
                         wi_index < 4 * num_wi;) {
                    Vec3<Float> wi = assembly.bsdfPdfSample(pcg, wo);
                    if ((wi == 0).all()) {
                        continue;
                    }
                    Float wi_f = 0;
                    Float wi_fpdf = 0;
                    for (int iter = 0; 
                             iter < 512; iter++) {
                        Float f = 0;
                        Float fpdf = assembly.bsdfPdf(pcg, wo, wi, &f);
                    /*  wi_f += f * 
                            (Float(1) / Float(512));
                        wi_fpdf += fpdf *
                            (Float(1) / Float(512));  */
                        wi_f = wi_f + (f - wi_f) / (iter + 1);
                        wi_fpdf = wi_fpdf + (fpdf - wi_fpdf) / (iter + 1);
                    }
                    if (wi_f > 0 && 
                        wi_fpdf > 0) {
                    /*  wi_ffac += wi_f / wi_fpdf;  */
                        wi_ffac = wi_ffac + 
                                 (wi_f / wi_fpdf - wi_ffac) / (wi_index + 1);
                        wi_input.push_back({wi, wi_f});
                        wi_index++;
                    }
                }
            /*  wi_ffac /= 4 * num_wi;  */
                for (Rrss::Sample& wi_sample : wi_input) {
                    wi_sample.dir_pdf /= wi_ffac;
                }

                int wi_index = 0;
                Rrss rrss(wi_input);
                rrss.disableMostRedundant(3 * num_wi);
                for (const Rrss::Sample& wi_sample : rrss.samples()) {
                    if (wi_sample.is_enabled) {
                        wi_array[wo_index * num_wi + wi_index] = wi_sample.dir;
                        wi_index++;
                    }
                }
            }
            const Vec3<Float>* wi = &wi_array[wo_index * num_wi];
            Float* f = &f_array[wo_index * num_wi];
            std::vector<Float> f0tmp(num_wi);
            std::vector<Float> f1tmp(num_wi);
            std::vector<Float> f0(num_wi);
            std::vector<Float> f1(num_wi);
            for (int iter = 0; iter < num_iters / 2; iter++) {
                assembly.bsdf(pcg, wo, wi, num_wi, &f0tmp[0]);
                assembly.bsdf(pcg, wo, wi, num_wi, &f1tmp[0]);
                for (int j = 0; j < num_wi; j++) {
                    f0[j] = f0[j] + (f0tmp[j] - f0[j]) / (iter + 1);
                    f1[j] = f1[j] + (f1tmp[j] - f1[j]) / (iter + 1);
                }
                {
                    std::unique_lock<std::mutex> lock(cout_mutex);
                    num_complete += 2;
                    std::cout << '\r';
                    std::cout << 
                    pr::terminal_progress_bar{
                        double(num_complete) / 
                        double(num_wo * num_iters)};
                    std::cout.flush();
                }
            }
            for (int wi_index = 0; wi_index < num_wi; wi_index++) {
                f[wi_index] = (f0[wi_index] + f1[wi_index]) * Float(0.5);
                f[wi_index] /= pr::abs(wi[wi_index][2]);
            }
#if 0
            for (int wi_index = 0;
                     wi_index < num_wi; wi_index++) {
                const Vec3<Float>& wi = wi_array[wo_index * num_wi + wi_index];
                Float& f = f_array[wo_index * num_wi + wi_index];
                f = 0;
                for (int iter = 0; 
                         iter < num_iters; iter++) {
                /*  f += assembly.bsdf(pcg, wo, wi);  */
                    f = f + (assembly.bsdf(pcg, wo, wi) - f) / (iter + 1);
                }
            /*  f /= num_iters;  */
                f /= pr::abs(wi[2]);
                {
                    std::unique_lock<std::mutex> lock(cout_mutex);
                    num_complete++;
                    std::cout << '\r';
                    std::cout << 
                    pr::terminal_progress_bar{
                        double(num_complete) / 
                        double(num_wo * num_wi)};
                    std::cout.flush();
                }
            }
#endif
        };
    
        std::vector<std::future<void>> wait;
        wait.reserve(num_wo);
        pr::thread_pool pool;
        for (int wo_index = 0;
                 wo_index < num_wo; wo_index++) {
            wait.emplace_back(pool.submit(job, wo_index));
        }
        for (int wo_index = 0;
                 wo_index < num_wo; wo_index++) {
            wait[wo_index].wait();
        }   
    }

    {
        // Output filestream.
        std::ofstream ofs(ofs_fname);
        if (!ofs) {
            throw std::runtime_error("");
        }
        ofs << "RAWBH";
        ofs << "10A Layered-SQT\n";
        ofs << "1 0.5\n";
        ofs << num_wo;
        for (Float thetao : thetao_array) {
            ofs << ' ';
            ofs << thetao;
        }
        ofs << '\n';
        for (int wo_index = 0;
                 wo_index < num_wo; wo_index++) {
            for (int wi_index = 0;
                     wi_index < num_wi; wi_index++) {
                ofs << wo_index << ' ';
                ofs << wi_array[wo_index * num_wi + wi_index][0] << ' ';
                ofs << wi_array[wo_index * num_wi + wi_index][1] << ' ';
                ofs << wi_array[wo_index * num_wi + wi_index][2] << ' ';
                ofs << f_array[wo_index * num_wi + wi_index] << '\n';
            }
        }
    }

    std::exit(EXIT_SUCCESS);
    return 0;
}
