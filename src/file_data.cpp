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
#include <layered-sqt/file_data.hpp>

namespace ls {

// Initialize.
void FileData::init(int wo_count, int wi_count)
{
    clear();
    if (!(wo_count > 0 &&
          wi_count > 0)) {
        return;
    }

    slices_.resize(wo_count);
    int wo_index = 0;
    for (Slice& slice : slices_) {
        slice.init(
                wo_index / Float(wo_count) * 
                pr::numeric_constants<Float>::M_pi_2(),
                wi_count);
        wo_index++;
    }
}

// Clear.
void FileData::clear()
{
    for (Slice& slice : slices_) {
        slice.clear();
    }

    slices_.clear();
    slices_.shrink_to_fit();
}

// Write SQT RAW.
void FileData::writeSqtRaw(std::ostream& ostr) const
{
    bool is_bsdf = false;
    for (const Slice& slice : slices_) {
        if (slice.anyTransmitted()) {
            is_bsdf = true;
            break;
        }
    }
    ostr << "RAWB" << (is_bsdf ? 'S' : 'H');
    ostr << "10A Layered-SQT\n";
    ostr << "1 0.5\n";
    ostr << slices_.size();
    for (const Slice& slice : slices_) {
        ostr << ' ';
        ostr << slice.thetao_;
    }
    ostr << '\n';
    int wo_index = 0;
    for (const Slice& slice : slices_) {
        for (int wi_index = 0;
                 wi_index < slice.wi_count_; wi_index++) {
            ostr << wo_index << ' ';
            ostr << slice.wi_[wi_index][0] << ' ';
            ostr << slice.wi_[wi_index][1] << ' ';
            ostr << slice.wi_[wi_index][2] << ' ';
            ostr << slice.f_[wi_index] /
                    pr::fabs(slice.wi_[wi_index][2]) << '\n';
        }
    }
}

// Initialize.
void FileData::Slice::init(Float thetao, int wi_count)
{
    clear();

    // Setup outgoing direction.
    thetao_ = thetao;
    wo_ = {
        pr::sin(thetao_), Float(0),
        pr::cos(thetao_)
    };

    // Setup incident directions.
    if (wi_count > 0) {
        wi_count_ = wi_count;
        wi_ = new Vec3<Float>[wi_count_];
        f_  = new Float[wi_count_];
    }
}

// Clear.
void FileData::Slice::clear()
{
    // Deallocate incident directions.
    delete[] wi_;

    // Deallocate BSDF values.
    delete[] f_;

    // Nullify.
    *this = Slice();
}

// Sample incident directions.
void FileData::Slice::sampleIncidentDirs(
                LayeredAssembly& layered_assembly,
                ProgressBar& progress_bar,
                int rrss_oversampling,
                int rrss_path_count, 
                int rrss_path_count_per_iter)
{
    assert(rrss_oversampling >= 1 &&
           rrss_path_count >= 1 &&
           rrss_path_count_per_iter >= 1);

    // RRSS incident direction count.
    int rrss_wi_count = 
        rrss_oversampling * wi_count_;

    // RRSS incident directions.
    Vec3<Float>* rrss_wi = new Vec3<Float>[rrss_wi_count];

    // Initialize incident directions.
    for (int rrss_wi_index = 0;
             rrss_wi_index < rrss_wi_count; rrss_wi_index++) {
    
        rrss_wi[rrss_wi_index] = 
        layered_assembly.randomScatterDirection(path_pcg_, wo_); 
    }

    // RRSS BSDFs.
    Float* rrss_f = new Float[rrss_wi_count];

    // RRSS BSDF-PDFs.
    Float* rrss_f_pdf = new Float[rrss_wi_count];

    for (int curr_path_count = 0;
             curr_path_count < rrss_path_count;
             curr_path_count += rrss_path_count_per_iter) {

        int next_path_count = curr_path_count + rrss_path_count_per_iter;
        if (next_path_count > rrss_path_count) {
            next_path_count = rrss_path_count;
        }

        // Update BSDF/BSDF-PDF averages.
        layered_assembly.computeAverage(
                curr_path_count,
                next_path_count,
                path_pcg_, wo_,
                rrss_wi,
                rrss_wi_count,
                rrss_f, rrss_f_pdf);

        // Update progress bar.
        progress_bar.add(
                next_path_count - 
                curr_path_count);
    }

    // RRSS samples.
    std::vector<Rrss::Sample> rrss_samples;
    rrss_samples.reserve(rrss_wi_count);

    // RRSS BSDF integral.
    Float rrss_f_int = 0;

    for (int rrss_wi_index = 0;
             rrss_wi_index < rrss_wi_count; 
             rrss_wi_index++) {

        // Push sample.
        rrss_samples.push_back(
        Rrss::Sample{
            rrss_wi[rrss_wi_index],
            rrss_f [rrss_wi_index] // BSDF as initial PDF.
        });

        // BSDF integral term.
        Float f_int = 
            rrss_f[rrss_wi_index] /
            rrss_f_pdf[rrss_wi_index];

        // BSDF integral term okay?
        if (pr::isfinite(f_int)) {

            // Update BSDF integral estimate.
            rrss_f_int =
            rrss_f_int + (f_int - rrss_f_int) / 
                    (rrss_wi_index + 1);
        }
    }

    // BSDF integral okay?
    if (rrss_f_int > 0 &&
        pr::isfinite(rrss_f_int)) {

        // Normalize.
        for (Rrss::Sample& rrss_sample : rrss_samples) {
            rrss_sample.pdf /= rrss_f_int;
        }
    }
    else {

        // Default to BSDF-PDF.
        int rrss_wi_index = 0;
        for (Rrss::Sample& rrss_sample : rrss_samples) {
            rrss_sample.pdf = rrss_f_pdf[rrss_wi_index++];
        }
    }

    // Redundancy-reduced sample set.
    Rrss rrss(rrss_samples);
    if (!rrss.disableUntil(wi_count_)) {
        // Unreachable?
    }

    int wi_index = 0;
    for (int rrss_wi_index = 0;
             rrss_wi_index < rrss_wi_count; 
             rrss_wi_index++) {

        // Is sample enabled?
        if (rrss.samples()[rrss_wi_index].is_enabled) {

            // Use BSDF and incident direction.
            f_ [wi_index] = rrss_f [rrss_wi_index];
            wi_[wi_index] = rrss_wi[rrss_wi_index];
            wi_index++;
        }
    }

    // Remember path count.
    path_count_ = rrss_path_count;

    // Delete RRSS BSDF-PDFs.
    delete[] rrss_f_pdf;

    // Delete RRSS BSDFs.
    delete[] rrss_f;

    // Delete RRSS incident directions.
    delete[] rrss_wi;
}

// Update BSDF averages.
void FileData::Slice::updateBsdfAverages(
                LayeredAssembly& layered_assembly,
                ProgressBar& progress_bar,
                int path_count,
                int path_count_per_iter)
{
    assert(path_count >= 1 &&
           path_count_per_iter >= 1);

    int& this_path_count = path_count_;
    for (int target_path_count = this_path_count + path_count;
             target_path_count > this_path_count;) {

        // Next path count.
        int next_path_count = this_path_count + path_count_per_iter;
        if (next_path_count > target_path_count) {
            next_path_count = target_path_count;
        }

        // Update BSDF averages.
        layered_assembly.computeAverage(
                this_path_count, 
                next_path_count,
                path_pcg_, wo_,
                wi_, wi_count_, f_, 
                nullptr);

        // Update progress bar.
        progress_bar.add(
                next_path_count - 
                this_path_count);

        // Increment path count.
        this_path_count = next_path_count;
    }
}

} // namespace ls
