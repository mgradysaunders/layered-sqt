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
#include <sstream>
#include <layered-sqt/layered_assembly.hpp>
#include <layered-sqt/layer/null.hpp>
#include <layered-sqt/layer/lambertian.hpp>
#include <layered-sqt/layer/microsurface_lambertian.hpp>
#include <layered-sqt/layer/microsurface_dielectric.hpp>
#include <layered-sqt/layer/microsurface_conductive.hpp>
#include <layered-sqt/layer/oren_nayar_diffuse.hpp>
#include <layered-sqt/medium/henyey_greenstein.hpp>
#include <layered-sqt/medium/henyey_greenstein2.hpp>
#include <layered-sqt/medium/rayleigh.hpp>
#include <layered-sqt/medium/sggx.hpp>

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

                Float mua = 0;
                Float mus = 0;
                Float eta = 1;
                std::string medium_phase;
                try {
                    // Read parameters.
                    // Note std::stod() throws if string is invalid.
                    while (iss >> str) {
                        if (!str.compare(0, 4, "eta=", 4)) {
                            eta = std::stod(str.substr(4));
                        }
                        else
                        if (!str.compare(0, 4, "mua=", 4)) {
                            mua = std::stod(str.substr(4));
                        }
                        else 
                        if (!str.compare(0, 4, "mus=", 4)) {
                            mus = std::stod(str.substr(4));
                        }
                        else {
                            if (str == "HenyeyGreenstein" ||
                                str == "HenyeyGreenstein2" ||
                                str == "Rayleigh" ||
                                str == "Sggx") {
                                medium_phase = str;
                                break;
                            }
                            else {
                                // Trigger catch block.
                                throw std::exception();
                            }
                        }
                    }
                }
                catch (const std::exception&) {
                    std::runtime_error(
                    std::string(__PRETTY_FUNCTION__)
                        .append(": invalid argument '").append(str)
                        .append("'"));
                }

                const char* error_message = nullptr;
                if (!(eta > 0)) {
                    error_message = ": eta is non-positive";
                }
                else
                if (!(mua >= 0)) {
                    error_message = ": mua is negative";
                }
                else 
                if (!(mus >= 0)) {
                    error_message = ": mus is negative";
                }
                else
                if (mus > 0 && medium_phase.empty()) {
                    error_message = ": mus is positive, but no phase";
                }
                // Error?
                if (error_message) {
                    throw 
                        std::runtime_error(
                        std::string(__PRETTY_FUNCTION__).append(error_message));
                }

                // Medium phase identifier.
                if (medium_phase.empty()) {
                    mediums_.push_back(new Medium()); // Default.
                }
                else 
                if (medium_phase == "HenyeyGreenstein") {
                    mediums_.push_back(new HenyeyGreensteinMedium());
                }
                else 
                if (medium_phase == "HenyeyGreenstein2") {
                    mediums_.push_back(new HenyeyGreenstein2Medium());
                }
                else
                if (medium_phase == "Rayleigh") {
                    mediums_.push_back(new RayleighMedium());
                }
                else
                if (medium_phase == "Sggx") {
                    mediums_.push_back(new SggxMedium());
                }

                // Medium parameters.
                mediums_.back()->eta = eta;
                mediums_.back()->mua = mua;
                mediums_.back()->mus = mus;
                mediums_.back()->mu = mua + mus;

                // Medium phase parameters
                if (!medium_phase.empty()) {
                    // Delegate.
                    std::string arg;
                    std::getline(iss, arg);
                    mediums_.back()->init(arg);
                }
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
                if (str == "Null") {
                    layers_.push_back(new NullLayer());
                }
                else
                if (str == "Lambertian") {
                    layers_.push_back(new LambertianLayer());
                }
                else
                if (str == "MicrosurfaceLambertian") {
                    layers_.push_back(new MicrosurfaceLambertianLayer());
                }
                else
                if (str == "MicrosurfaceDielectric") {
                    layers_.push_back(new MicrosurfaceDielectricLayer());
                }
                else
                if (str == "MicrosurfaceConductive") {
                    layers_.push_back(new MicrosurfaceConductiveLayer());
                }
                else
                if (str == "OrenNayarDiffuse") {
                    layers_.push_back(new OrenNayarDiffuseLayer());
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

    const char* error_message = nullptr;

    // Expect medium?
    if (expect_medium) {
        error_message = ": expected bottom medium";
    }
    else 
    // Layers empty?
    if (layers_.empty()) {
        error_message = ": no layers";
    }
    else
    // Boundary medium has non-zero extinction?
    if (mediums_.front()->mu > 0 ||
        (mediums_.back()->mu > 0 && isTransmissive())) {
        error_message = ": boundary medium has non-zero extinction";
    }

    // Error?
    if (error_message) {
        throw
            std::runtime_error(
            std::string(__PRETTY_FUNCTION__)
                .append(error_message));
    }

    // Link pointers.
    for (std::size_t pos = 0; pos < layers_.size(); pos++) {
        mediums_[pos + 0]->layer_below = layers_[pos];
        mediums_[pos + 1]->layer_above = layers_[pos];
        layers_[pos]->medium_above = mediums_[pos + 0];
        layers_[pos]->medium_below = mediums_[pos + 1];
    }

    // Find non-null BSDF layer at top.
    for (auto itr = layers_.begin();
              itr < layers_.end(); itr++) {
        if (!(*itr)->isNull()) {
            layers_top_ = *itr;
            break;
        }
    }

    // Find non-null BSDF layer at bottom.
    for (auto itr = layers_.rbegin();
              itr < layers_.rend(); itr++) {
        if (!(*itr)->isNull()) {
            layers_bottom_ = *itr;
            break;
        }
    }

    // Validate layers.
    for (const Layer* layer : layers_) {
        layer->validate(); // May throw.
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
            Float* fs,
            Float* fs_pdf) const 
{
    assert(wi_count > 0);
    assert(wi && (fs || fs_pdf));

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
    if (fs) {
        // Zero initialize.
        for (int j = 0; j < wi_count; j++) {
            fs[j] = 0;
        }
    }

    // BSDF-PDF?
    if (fs_pdf) {
        // Zero initialize.
        for (int j = 0; j < wi_count; j++) {
            fs_pdf[j] = 0;
        }
    }

    Float tau = 1;
    for (int bounce = 0;
             bounce < 256; bounce++) {

        // Re-center since homogeneous.
        ray.pos[0] = 0;
        ray.pos[1] = 0;

        Vec3<Float> wk = -pr::normalize_fast(ray.dir);

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
        {
            const Layer* layer_next = ray.dir[2] > 0 
                ? layer_above 
                : layer_below;
            if (!layer_next ||
                !layer_next->intersect(ray, hit)) {
                break;
            }
            layer = layer_next;
        }

        // Sample medium.
        Float dmax = pr::length(hit.pos - ray.pos);
        Float d = ray.medium->transmittanceSample(pcg, tau, dmax);
        if (!(d == dmax)) {

            assert(ray.medium->layer_above);
            assert(ray.medium->layer_below);

            // Path absorbed?
            if (tau == 0) {
                // Terminate.
                break;
            }

            // Update hit.
            hit.pos = ray.pos + ray.dir * d;

            // Above top?
            if (!layers_top_ ||
                hit.pos[2] > layers_top_->zheight) {
                for (int j = 0; j < wi_count; j++) {

                    // Incident direction in upper hemisphere?
                    if (pr::signbit(wi[j][2]) == 0) {

                        // Phase.
                        Float ps = ray.medium->phase(pcg, wk, wi[j]);
                        Float tmp__fs = ps;
                        Float tmp__fs_pdf = ps;

                        // Transmittance.
                        Float tr = 1;
                        Float z0 = hit.pos[2];
                        const Layer* layer_next = ray.medium->layer_above;
                        while (layer_next) {
                            tr *= 
                                layer_next->medium_below->transmittance(
                               (layer_next->zheight - z0) / wi[j][2]);
                            z0 = layer_next->zheight;
                            layer_next = 
                            layer_next->medium_above->layer_above;
                        }

                        // BSDF?
                        tmp__fs *= tau * tr;
                        if (fs && pr::isfinite(tmp__fs)) {
                            fs[j] += tmp__fs;
                        }

                        // BSDF-PDF?
                        tmp__fs_pdf *= tr;
                        if (fs_pdf && pr::isfinite(tmp__fs_pdf)) {
                            fs_pdf[j] += tmp__fs_pdf;
                        }
                    }
                }
            }

            // Below bottom?
            if (!layers_bottom_ ||
                hit.pos[2] < layers_bottom_->zheight) {
                for (int j = 0; j < wi_count; j++) {

                    // Incident direction in lower hemisphere?
                    if (pr::signbit(wi[j][2]) == 1) {

                        // Phase.
                        Float ps = ray.medium->phase(pcg, wk, wi[j]);
                        Float tmp__fs = ps;
                        Float tmp__fs_pdf = ps;

                        // Transmittance.
                        Float tr = 1;
                        Float z0 = hit.pos[2];
                        const Layer* layer_next = ray.medium->layer_below;
                        while (layer_next) {
                            tr *= 
                                layer_next->medium_above->transmittance(
                               (layer_next->zheight - z0) / wi[j][2]);
                            z0 = layer_next->zheight;
                            layer_next = 
                            layer_next->medium_below->layer_below;
                        }

                        // BSDF?
                        tmp__fs *= tau * tr;
                        if (fs && pr::isfinite(tmp__fs)) {
                            fs[j] += tmp__fs;
                        }

                        // BSDF-PDF?
                        tmp__fs_pdf *= tr;
                        if (fs_pdf && pr::isfinite(tmp__fs_pdf)) {
                            fs_pdf[j] += tmp__fs_pdf;
                        }
                    }
                }
            }

            // Update ray.
            ray.pos = hit.pos;
            ray.dir = ray.medium->phaseSample(pcg, tau, wk);
        }
        else {

            // At top?
            if (layer == layers_top_) {
                for (int j = 0; j < wi_count; j++) {

                    // Incident direction in upper hemisphere?
                    if (pr::signbit(wi[j][2]) == 0) {

                        // BSDF/BSDF-PDF.
                        Float tmp__fs_pdf = 0;
                        Float tmp__fs = 
                                layer->bsdf(pcg, wk, wi[j], &tmp__fs_pdf);

                        // Transmittance.
                        Float tr = 1;
                        Float z0 = hit.pos[2];
                        const Layer* layer_next = layer_above;
                        while (layer_next) {
                            tr *= 
                                layer_next->medium_below->transmittance(
                               (layer_next->zheight - z0) / wi[j][2]);
                            z0 = layer_next->zheight;
                            layer_next = 
                            layer_next->medium_above->layer_above;
                        }

                        // BSDF?
                        tmp__fs *= tau * tr;
                        if (fs && pr::isfinite(tmp__fs)) {
                            fs[j] += tmp__fs;
                        }

                        // BSDF-PDF?
                        tmp__fs_pdf *= tr;
                        if (fs_pdf && pr::isfinite(tmp__fs_pdf)) {
                            fs_pdf[j] += tmp__fs_pdf;
                        }
                    }
                }
            }

            // At bottom?
            if (layer == layers_bottom_) {
                for (int j = 0; j < wi_count; j++) {

                    // Incident direction in lower hemisphere?
                    if (pr::signbit(wi[j][2]) == 1) {

                        // BSDF/BSDF-PDF.
                        Float tmp__fs_pdf = 0;
                        Float tmp__fs = 
                                layer->bsdf(pcg, wk, wi[j], &tmp__fs_pdf);

                        // Transmittance.
                        Float tr = 1;
                        Float z0 = hit.pos[2];
                        const Layer* layer_next = layer_below;
                        while (layer_next) {
                            tr *= 
                                layer_next->medium_above->transmittance(
                               (layer_next->zheight - z0) / wi[j][2]);
                            z0 = layer_next->zheight;
                            layer_next = 
                            layer_next->medium_below->layer_below;
                        }

                        // BSDF?
                        tmp__fs *= tau * tr;
                        if (fs && pr::isfinite(tmp__fs)) {
                            fs[j] += tmp__fs;
                        }

                        // BSDF-PDF?
                        tmp__fs_pdf *= tr;
                        if (fs_pdf && pr::isfinite(tmp__fs_pdf)) {
                            fs_pdf[j] += tmp__fs_pdf;
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
        if (!(tau > 0) || ray.dir[2] == 0) {
            break;
        }
    }
}

// Compute BSDF/BSDF-PDF average.
void LayeredAssembly::computeAverage(
            int path_count_before,
            int path_count,
            Pcg32& pcg,
            const Vec3<Float>& wo,
            const Vec3<Float>* wi, int wi_count,
            Float* fs,
            Float* fs_pdf) const
{
    assert(path_count > 
           path_count_before);

    assert(wi_count > 0);
    assert(wi && (fs || fs_pdf));

    // Temporaries.
    Float* tmp = new Float[4 * wi_count];

    // Temporary BSDFs.
    Float* tmp1__fs = fs ? tmp + 0 * wi_count : nullptr;
    Float* tmp2__fs = fs ? tmp + 1 * wi_count : nullptr;

    // Temporary BSDF-PDFs.
    Float* tmp1__fs_pdf = fs_pdf ? tmp + 2 * wi_count : nullptr;
    Float* tmp2__fs_pdf = fs_pdf ? tmp + 3 * wi_count : nullptr;

    for (int path = 0;
             path < path_count - path_count_before; path++) {

        compute(pcg, wo, wi, wi_count, tmp1__fs, tmp1__fs_pdf);

        // BSDF?
        if (fs) {
            // Update.
            for (int j = 0; j < wi_count; j++) {
                tmp2__fs[j] = 
                tmp2__fs[j] + 
                    (tmp1__fs[j] - tmp2__fs[j]) / (path + 1);
            }
        }

        // BSDF-PDF?
        if (fs_pdf) {
            // Update.
            for (int j = 0; j < wi_count; j++) {
                tmp2__fs_pdf[j] = 
                tmp2__fs_pdf[j] + 
                    (tmp1__fs_pdf[j] - tmp2__fs_pdf[j]) / (path + 1);
            }
        }
    }

    // Factor to combine averages.
    Float fac = 
        Float(path_count_before) / 
        Float(path_count);

    // BSDF?
    if (fs) {
        if (path_count_before == 0) {
            // Initialize.
            for (int j = 0; j < wi_count; j++) {
                fs[j] = tmp2__fs[j];
            }
        }
        else {
            // Update.
            for (int j = 0; j < wi_count; j++) {
                fs[j] = 
                fs[j] * fac + tmp2__fs[j] * (1 - fac);
            }
        }
    }

    // BSDF-PDF?
    if (fs_pdf) {
        if (path_count_before == 0) {
            // Initialize.
            for (int j = 0; j < wi_count; j++) {
                fs_pdf[j] = tmp2__fs_pdf[j];
            }
        }
        else {
            // Update.
            for (int j = 0; j < wi_count; j++) {
                fs_pdf[j] = 
                fs_pdf[j] * fac + tmp2__fs_pdf[j] * (1 - fac);
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
    for (int attempt = 0;
             attempt < 4096; attempt++) {

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

            Vec3<Float> wk = -pr::normalize_fast(ray.dir);

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

                // Reject nearly parallel rays.
                if (!(pr::fabs(ray.dir[2]) > Float(1e-3))) {
                    break;
                }
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

                // Path absorbed?
                if (tau == 0) {
                    // Terminate.
                    break;
                }

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
            if (!(tau > 0) || ray.dir[2] == 0) {
                // Terminate.
                break;
            }
        }

        // Try again.
    }

    throw std::runtime_error(
          std::string(__PRETTY_FUNCTION__)
              .append(": failed to sample direction in 4096 attempts"));

    // Unreachable.
    return Vec3<Float>{};
}

} // namespace ls
