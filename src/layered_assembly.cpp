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
#include <layered-sqt/layered_assembly.hpp>
#include <layered-sqt/layer/lambertian.hpp>
#include <layered-sqt/layer/microsurface_lambertian.hpp>
#include <layered-sqt/layer/microsurface_dielectric.hpp>

namespace ls {

// Initialize from input stream.
void LayeredAssembly::init(std::istream& is)
{
    clear();

    long lineno = 1;
    bool expect_medium = true;
    std::string strln;
    while (std::getline(is, strln)) {
        try {
            std::istringstream iss(strln);
            std::string str;


            if (expect_medium) {
                // No 'Medium'?
                if (!(iss >> str) ||
                    !(str == "Medium")) {
                    throw 
                        std::runtime_error(
                        std::string(": expected 'Medium'"));
                }

                // Delegate.
                std::string arg;
                std::getline(iss, arg);
                mediums_.push_back(new Medium());
                mediums_.back()->init(arg);
            }
            else {
                // No 'Layer'?
                if (!(iss >> str) ||
                    !(str == "Layer")) {
                    throw 
                        std::runtime_error(
                        std::string(": expected 'Layer'"));
                }

                // No 'z='?
                Float zheight = 0;
                if (!(iss >> str) ||
                    str.compare(0, 2, "z=", 2) != 0) {
                    // Runtime error.
                    throw
                        std::runtime_error(
                        std::string(__PRETTY_FUNCTION__)
                            .append(": expected 'z='"));
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
                            .append(": invalid parameter '").append(str)
                            .append("'"));
                }
                if (layers_.size() > 0 &&
                    !(layers_.back()->zheight > zheight)) {
                    // Runtime error.
                    throw
                        std::runtime_error(
                        std::string(__PRETTY_FUNCTION__)
                            .append(": layer z-height is out-of-order"));
                }

                // Read identifier.
                iss >> str;
                if (str == "LambertianBsdf") {
                    layers_.push_back(new LambertianBsdfLayer());
                }
                else
                if (str == "MicrosurfaceLambertianBrdf") {
                    layers_.push_back(new MicrosurfaceLambertianBrdfLayer());
                }
                else
                if (str == "MicrosurfaceDielectricBsdf") {
                    layers_.push_back(new MicrosurfaceDielectricBsdfLayer());
                }
                else {
                    // Runtime error.
                    throw 
                        std::runtime_error(
                        std::string(__PRETTY_FUNCTION__)
                            .append(": expected layer identifier"));
                }

                // Delegate.
                std::string arg;
                std::getline(iss, arg);
                layers_.back()->init(arg);
                layers_.back()->zheight = zheight;
            }

            lineno++;
            expect_medium = !expect_medium;
        }
        catch (const std::exception& exception) {
            throw 
                std::runtime_error(
                std::string(__PRETTY_FUNCTION__)
                    .append(": on line ").append(std::to_string(lineno))
                    .append(", ").append(exception.what()));
        }
    }

    // Expect medium?
    if (expect_medium) {
        throw 
            std::runtime_error(
            std::string(__PRETTY_FUNCTION__)
                .append(": expected bottom medium"));

    }

    // No layers?
    if (layers_.empty()) {
        throw 
            std::runtime_error(
            std::string(__PRETTY_FUNCTION__)
                .append(": no layers"));
    }

    // Link pointers.
    for (std::size_t pos = 0;
                     pos < layers_.size(); pos++) {
        mediums_[pos + 0]->layer_below = layers_[pos];
        mediums_[pos + 1]->layer_above = layers_[pos];
        layers_[pos]->medium_above = mediums_[pos + 0];
        layers_[pos]->medium_below = mediums_[pos + 1];
    }
}

// Clear.
void LayeredAssembly::clear()
{
    // Destroy mediums.
    for (Medium* medium : mediums_) {
        delete medium;
    }
    mediums_.clear();
    mediums_.shrink_to_fit();

    // Destroy layers.
    for (Layer* layer : layers_) {
        delete layer;
    }
    layers_.clear();
    layers_.shrink_to_fit();
}

// Compute BSDF/BSDF-PDF.
void LayeredAssembly::compute(
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>* wi, int wi_count,
            Float* f,
            Float* f_pdf) const 
{
    assert(wi_count > 0);
    assert(wi && (f || f_pdf));

    // Initial layer.
    const Layer* layer;

    // Initial ray.
    Ray ray;
    ray.dir = -wo;

    // Upper hemisphere?
    if (wo[2] > 0) {
        // Move above top layer.
        layer = layers_.front();
        ray.pos[2] = layer->zheight + 1;
        ray.medium = layer->medium_above;
    }
    else {
        // Move below bottom layer.
        layer = layers_.back();
        ray.pos[2] = layer->zheight - 1;
        ray.medium = layer->medium_below;
    }

    // BSDF?
    if (f) {
        // Zero initialize.
        for (int j = 0; j < wi_count; j++) {
            f[j] = 0;
        }
    }

    // BSDF-PDF?
    if (f_pdf) {
        // Zero initialize.
        for (int j = 0; j < wi_count; j++) {
            f_pdf[j] = 0;
        }
    }

    Float tau = 1;
    for (int bounce = 0;
             bounce < 128; bounce++) {

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
            if (layer == layers_.front()) {
                for (int j = 0; j < wi_count; j++) {
                    if (pr::signbit(wi[j][2]) == 0) {
                        Float tmp_f_pdf;
                        Float tmp_f = layer->bsdf(pcg, wk, wi[j], &tmp_f_pdf);
                        if (f) {
                            f[j] += tau * tmp_f;
                        }
                        if (f_pdf) {
                            f_pdf[j] += tmp_f_pdf;
                        }
                    }
                }
            }
            else 
            if (layer == layers_.back()) {
                for (int j = 0; j < wi_count; j++) {
                    if (pr::signbit(wi[j][2]) == 1) {
                        Float tmp_f_pdf;
                        Float tmp_f = layer->bsdf(pcg, wk, wi[j], &tmp_f_pdf);
                        if (f) {
                            f[j] += tau * tmp_f;
                        }
                        if (f_pdf) {
                            f_pdf[j] += tmp_f_pdf;
                        }
                    }
                }
            }

            // Update ray.
            ray.pos = hit.pos;
            ray.dir = layer->bsdfSample(pcg, tau, wk);
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

// Compute BSDF/BSDF-PDF average.
void LayeredAssembly::computeAverage(
            int path_count,
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>* wi, int wi_count,
            Float* f,
            Float* f_pdf) const
{
    assert(wi_count > 0);
    assert(wi && (f || f_pdf));

    // BSDF?
    if (f) {
        // Zero initialize.
        for (int j = 0; j < wi_count; j++) {
            f[j] = 0;
        }
    }

    // BSDF-PDF?
    if (f_pdf) {
        // Zero initialize.
        for (int j = 0; j < wi_count; j++) {
            f_pdf[j] = 0;
        }
    }

    // Temporaries.
    Float* tmp = new Float[2 * wi_count];

    // Temporary BSDFs.
    Float* tmp_f = f ? tmp : nullptr;

    // Temporary BSDF-PDFs.
    Float* tmp_f_pdf = f_pdf ? tmp + wi_count : nullptr;

    for (int path = 0;
             path < path_count; path++) {

        compute(pcg, wo, wi, wi_count, tmp_f, tmp_f_pdf);

        // BSDF?
        if (f) {
            // Update.
            for (int j = 0; j < wi_count; j++) {
                f[j] = 
                f[j] + (tmp_f[j] - f[j]) / (path + 1);
            }
        }

        // BSDF-PDF?
        if (f_pdf) {
            // Update.
            for (int j = 0; j < wi_count; j++) {
                f_pdf[j] = 
                f_pdf[j] + (tmp_f_pdf[j] - f_pdf[j]) / (path + 1);
            }
        }
    }
    
    // Delete temporaries.
    delete[] tmp;
}

// Random scatter direction.
Vec3<Float> LayeredAssembly::randomScatterDirection(
                Pcg32& pcg, 
                const Vec3<Float>& wo) const
{
    for (;;) {

        // Initial layer.
        const Layer* layer;

        // Initial ray.
        Ray ray;
        ray.dir = -wo;

        // Upper hemisphere?
        if (wo[2] > 0) {
            // Move above top layer.
            layer = layers_.front();
            ray.pos[2] = layer->zheight + 1;
            ray.medium = layer->medium_above;
        }
        else {
            // Move below bottom layer.
            layer = layers_.back();
            ray.pos[2] = layer->zheight - 1;
            ray.medium = layer->medium_below;
        }

        Float tau = 1;
        for (int bounce = 0;
                 bounce < 128; bounce++) {

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

                // Update ray.
                ray.pos = hit.pos;
                ray.dir = layer->bsdfSample(pcg, tau, wk);
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

        // Try again.
    }

    // Unreachable.
    return Vec3<Float>{};
}

} // namespace ls
