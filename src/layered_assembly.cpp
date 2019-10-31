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
#include <layered-sqt/layer/dielectric_microsurface.hpp>

namespace ls {

// Load.
void LayeredAssembly::load(std::istream& is)
{
    clear();

    long lineno = 1;
    bool expect_medium = true;
    std::string strln;
    while (std::getline(is, strln)) {
        try {
            if (expect_medium) {
                mediums_.push_back(new Medium());
                mediums_.back()->load(strln);
            }
            else {
                std::stringstream iss(strln);
                std::string str;

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
                if (str == "DielectricMicrosurfaceBsdf") {
                    layers_.push_back(new DielectricMicrosurfaceBsdfLayer());
                }
                else {
                    // Runtime error.
                    throw 
                        std::runtime_error(
                        std::string(__PRETTY_FUNCTION__)
                            .append(": expected "
                                    "'LambertianBsdf' or "
                                    "'DielectricMicrosurfaceBsdf'"));
                }

                // Delegate.
                std::getline(iss, strln);
                layers_.back()->load(strln);
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

// BSDF.
void LayeredAssembly::compute(
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>* wi, int n,
            Float* f,
            Float* fpdf) const 
{
    assert(n && wi);
    assert(f || fpdf);

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

    for (int j = 0; j < n; j++) {
        // Zero-initialize BSDF.
        if (f) {
            f[j] = 0;
        }
        // Zero-initialize BSDF-PDF.
        if (fpdf) {
            fpdf[j] = 0;
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
                for (int j = 0; j < n; j++) {
                    if (wi[j][2] > 0) {

                        // Update BSDF.
                        if (f) {
                            f[j] += tau * layer->bsdf(pcg, wk, wi[j]);
                        }

                        // Update BSDF-PDF.
                        if (fpdf) {
                            fpdf[j] += layer->bsdfPdf(pcg, wk, wi[j]);
                        }
                    }
                }
            }
            else 
            if (layer == layers_.back()) {
                for (int j = 0; j < n; j++) {
                    if (wi[j][2] < 0) {

                        // Update BSDF.
                        if (f) {
                            f[j] += tau * layer->bsdf(pcg, wk, wi[j]);
                        }

                        // Update BSDF-PDF.
                        if (fpdf) {
                            fpdf[j] += layer->bsdfPdf(pcg, wk, wi[j]);
                        }
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

// Random scatter direction.
Vec3<Float> LayeredAssembly::randomScatterDirection(
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

} // namespace ls
