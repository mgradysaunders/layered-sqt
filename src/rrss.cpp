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
#include <layered-sqt/rrss.hpp>

namespace ls {

// Constructor.
Rrss::Rrss(const std::vector<Sample>& samples) : 
        samples_(samples),
        samples_remaining_(samples.size())
{
    for (auto itr = samples_.begin(); itr < samples_.end();) {
        if (!(itr->dir_pdf >= pr::numeric_limits<Float>::min_invertible())) {
            itr = samples_.erase(itr);
            samples_remaining_--;
            // TODO Warn
        }
        else {
            itr->redundancy = 1;
            itr->is_enabled = true;
            ++itr;
        }
    }

    if (samples_.size() > 0) {
        sample_tree_.init(
        samples_.begin(), 
        samples_.end(),
        [](const Sample& sample) {
            return std::pair<Vec3<Float>, const Sample*>{sample.dir, &sample};
        });
    }
}

// Enabled sample directions.
std::vector<Vec3<Float>> Rrss::enabledSampleDirections() const
{
    std::vector<Vec3<Float>> sample_dirs;
    sample_dirs.reserve(samples_remaining_);
    for (const Sample& sample : samples_) {
        if (sample.is_enabled) {
            sample_dirs.push_back(sample.dir);
        }
    }
    return sample_dirs;
}

// Disable redundant samples.
bool Rrss::disable(std::size_t num)
{
    while (num-- > 0) {
        Sample* most_redundant = nullptr;
        for (Sample& sample : samples_) {

            // Is enabled?
            if (sample.is_enabled) {

                Float rad = 
                pr::sqrt(
                pr::numeric_constants<Float>::M_1_pi() /
                (sample.dir_pdf * samples_remaining_));

                if (!(rad > 0) || 
                    !pr::isfinite(rad)) {
                    sample.redundancy = 1;
                    // TODO Warn
                }
                else {
                    int nearby_enabled = 0;
                    int nearby = 0;
                    Float nearby_dir_pdf_sum = 0;

                    // Visit samples in region.
                    sample_tree_.nearby(
                    sample.dir, rad,
                    [&](auto nearby_node) {
                        const Sample& nearby_sample = 
                                     *nearby_node->value.second;
                        if (nearby_sample.is_enabled) {
                            nearby_enabled++;
                        }
                        nearby++;
                        nearby_dir_pdf_sum += nearby_sample.dir_pdf;
                        return true;
                    });

                    // Update redundancy.
                    sample.redundancy = 
                        nearby_enabled * 
                        nearby * (sample.dir_pdf / nearby_dir_pdf_sum);

                    if (pr::isnan(sample.redundancy)) {
                        sample.redundancy = 1;
                        // TODO Warn
                    }
                }

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
        }
        else {
            return false;
        }
    }
    return true;
}

} // namespace ls
