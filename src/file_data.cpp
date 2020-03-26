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
#include <random>
#include <preform/byte_order.hpp>
#include <layered-sqt/rrss.hpp>
#include <layered-sqt/tri.hpp>
#include <layered-sqt/file_data.hpp>

namespace ls {

// Basic initialize.
void FileData::basicInit(int wo_count, int wi_count, int seed)
{
    clear();
    if (!(wo_count > 0 &&
          wi_count > 0)) {
        return;
    }

    // Initialize slices.
    slices.resize(wo_count);

    int wo_index = 0;
    for (Slice& slice : slices) {

        // Outgoing angle.
        Float thetao = 
            wo_index / Float(wo_count) * 
            pr::numeric_constants<Float>::M_pi_2() * Float(0.99);

        // Initialize outgoing direction.
        slice.outgoing_angle = thetao;
        slice.outgoing_dir = {
            pr::sin(thetao), Float(0),
            pr::cos(thetao)
        };

        // Initialize incident directions.
        slice.incident_dirs.resize(wi_count);
        slice.bsdf_averages.resize(wi_count);
        slice.path_count = 0;

        // Initialize path PCG.
        slice.path_pcg = Pcg32(seed, wo_index);

        wo_index++;
    }
}

// Read LSQT-slice binary format.
void FileData::readLss(std::istream& istr)
{
    clear();

    // Wrap.
    pr::byte_stream_wrapper<std::istream> istr_le = 
    pr::byte_stream(istr, pr::byte_order::little);

    // Read magic word.
    std::uint32_t magic = 0;
    istr_le >> magic;
    if (magic != 0x4c535154UL) {
        throw std::runtime_error(__PRETTY_FUNCTION__);
    }

    // Read slice count.
    std::uint32_t slice_count = 0;
    istr_le >> slice_count;

    // Resize slices.
    slices.resize(slice_count);
    slices.shrink_to_fit();

    // Read slices.
    for (Slice& slice : slices) {
        slice.readLss(istr);
    }
}

// Write LSQT-slice binary format.
void FileData::writeLss(std::ostream& ostr) const
{
    // Wrap.
    pr::byte_stream_wrapper<std::ostream> ostr_le = 
    pr::byte_stream(ostr, pr::byte_order::little);

    // Write magic word.
    ostr_le << std::uint32_t(0x4c535154UL);

    // Write slice count.
    ostr_le << std::uint32_t(slices.size());

    // Write slices.
    for (const Slice& slice : slices) {
        slice.writeLss(ostr);
    }
}

// Write SQT RAW plain-text format.
void FileData::writeRaw(std::ostream& ostr, RawMode raw_mode) const
{
    bool is_bsdf = false;
    if (raw_mode == RAW_MODE_DEFAULT) {
        for (const Slice& slice : slices) {
            if (slice.anyTransmitted()) {
                is_bsdf = true;
                break;
            }
        }
    }
    ostr << "RAWB" << (is_bsdf ? 'S' : 'H');
    ostr << "10A Layered-SQT\n";
    ostr << "1 0.5\n";
    ostr << slices.size();
    for (const Slice& slice : slices) {
        ostr << ' ';
        ostr << slice.outgoing_angle;
    }
    ostr << '\n';
    int wo_index = 0;
    for (const Slice& slice : slices) {
        int wi_count = slice.incident_dirs.size();
        for (int wi_index = 0;
                 wi_index < wi_count; wi_index++) {
            // Get incident direction and BSDF average.
            auto incident_dir = slice.incident_dirs[wi_index];
            auto bsdf_average = slice.bsdf_averages[wi_index];

            // SQT expects non-cosine-weighted BSDF.
            bsdf_average /= pr::fabs(incident_dir[2]);
            if (!pr::isfinite(bsdf_average)) {
                continue;
            }

            // If BRDF only, ignore directions in lower hemisphere.
            // If BTDF only, ignore directions in upper hemisphere.
            if ((raw_mode == RAW_MODE_ONLY_BRDF && incident_dir[2] < 0.0f) ||
                (raw_mode == RAW_MODE_ONLY_BTDF_AS_BRDF &&
                 incident_dir[2] > 0.0f)) {
                continue;
            }
            // If BTDF only, flip to upper hemisphere.
            if (raw_mode == RAW_MODE_ONLY_BTDF_AS_BRDF) {
                incident_dir[2] = -incident_dir[2];
            }

            // Write sample.
            ostr << wo_index << ' ';
            ostr << +incident_dir[0] << ' ';
            ostr << +incident_dir[1] << ' ';
            ostr << +incident_dir[2] << ' ';
            ostr << bsdf_average << '\n';

            // Write Y-reflected sample. This assumes/enforces that the BSDF
            // is bilaterally symmetric. This should always be true given the
            // requirement that everything is homogeneous and isotropic. If
            // this is ever relaxed, this following code should be removed!
            ostr << wo_index << ' ';
            ostr << +incident_dir[0] << ' ';
            ostr << -incident_dir[1] << ' ';
            ostr << +incident_dir[2] << ' ';
            ostr << bsdf_average << '\n';
        }
        wo_index++;
    }
}

// Read LSQT-slice binary format.
void FileData::Slice::readLss(std::istream& istr)
{
    clear();

    // Wrap.
    pr::byte_stream_wrapper<std::istream> istr_le = 
    pr::byte_stream(istr, pr::byte_order::little);

    // Read outgoing angle.
    istr_le >> outgoing_angle;

    // Initialize outgoing direction.
    outgoing_dir = {
        pr::sin(outgoing_angle), Float(0),
        pr::cos(outgoing_angle)
    };

    // Read incident direction count.
    std::uint32_t wi_count = 0;
    istr_le >> wi_count;
    incident_dirs.resize(wi_count);
    bsdf_averages.resize(wi_count);

    // Read incident directions.
    for (auto& incident_dir : incident_dirs) {
        istr_le >> incident_dir[0];
        istr_le >> incident_dir[1];
        istr_le >> incident_dir[2];
    }

    // Read BSDF averages.
    for (auto& bsdf_average : bsdf_averages) {
        istr_le >> bsdf_average;
    }

    // Read path count.
    istr_le >> path_count;

    // Read path PCG.
    istr_le >> path_pcg.state_;
    istr_le >> path_pcg.inc_;
}

// Write LSQT-slice binary format.
void FileData::Slice::writeLss(std::ostream& ostr) const
{
    // Wrap.
    pr::byte_stream_wrapper<std::ostream> ostr_le = 
    pr::byte_stream(ostr, pr::byte_order::little);
    
    // Write outgoing angle.
    ostr_le << outgoing_angle;

    // Write incident direction count.
    ostr_le << std::uint32_t(incident_dirs.size());

    // Write incident directions.
    for (const auto& incident_dir : incident_dirs) {
        ostr_le << incident_dir[0];
        ostr_le << incident_dir[1];
        ostr_le << incident_dir[2];
    }

    // Write BSDF averages.
    for (const auto& bsdf_average : bsdf_averages) {
        ostr_le << bsdf_average;
    }

    // Write path count.
    ostr_le << path_count;

    // Write path PCG.
    ostr_le << path_pcg.state_;
    ostr_le << path_pcg.inc_;
}

// Smooth incident directions (presumably in upper hemisphere).
static 
void smoothIncidentDirs(std::vector<Vec3<Float>>& incident_dirs)
{
    if (incident_dirs.empty()) {
        return;
    }
    std::vector<Vec2<Float>> locs;
    locs.reserve(incident_dirs.size());

    // Initialize.
    for (const Vec3<Float>& incident_dir : incident_dirs) {
        locs.push_back(dirToPolar(incident_dir));
    }

    // Uniformly-spaced polar locations around edge of disk.
    for (int j = 0; j < 64; j++) {
        Float phi = j * pr::numeric_constants<Float>::M_pi() / 32;
        Float cos_phi = pr::cos(phi);
        Float sin_phi = pr::sin(phi);
        Vec2<Float> loc = {
            pr::numeric_constants<Float>::M_pi_2() * cos_phi,
            pr::numeric_constants<Float>::M_pi_2() * sin_phi
        };
        locs.push_back(loc);
    }

    std::vector<std::pair<Vec2<Float>, int>> smooth_locs(locs.size());
    for (auto& smooth_loc : smooth_locs) {
        smooth_loc.second = 0;
    }
    {
        // Triangulate.
        std::vector<Tri> tris;
        triangulate(locs, tris);

        // Compute unnormalized smooth locations.
        for (Tri tri : tris) {
            for (int k = 0; k < 3; k++) {
                int k0 = tri.indices[k];
                int k1 = tri.indices[(k + 1) % 3];
                int k2 = tri.indices[(k + 2) % 3];
                smooth_locs[k0].first += locs[k1] + locs[k2];
                smooth_locs[k0].second += 2;
            }
        }
    }

    // Overwrite directions.
    incident_dirs.clear();
    for (auto& smooth_loc : smooth_locs) {
        if (smooth_loc.second > 0) {
            smooth_loc.first /= smooth_loc.second; // Normalize.
            incident_dirs.push_back(polarToDir(smooth_loc.first));
        }
    }
}

// Compute incident directions.
void FileData::Slice::computeIncidentDirs(
                LayeredAssembly& layered_assembly,
                ProgressBar& progress_bar,
                int rrss_oversampling,
                int rrss_path_count, 
                int rrss_path_count_per_iter)
{
    assert(rrss_oversampling > 1 &&
           rrss_path_count >= 1 &&
           rrss_path_count_per_iter >= 1);
    if (incident_dirs.empty()) {
        return;
    }

    // Incident direction count.
    int wi_count = incident_dirs.size();

    // RRSS incident direction count.
    int rrss_wi_count = 
        rrss_oversampling * wi_count;

    // RRSS incident directions.
    Vec3<Float>* rrss_wi = nullptr;
    {
        // Incident directions.
        std::vector<Vec3<Float>> wi_upper;
        std::vector<Vec3<Float>> wi_lower;
        wi_upper.reserve(rrss_wi_count / 2);
        wi_lower.reserve(rrss_wi_count / 2);
        for (int rrss_wi_index = 0;
                 rrss_wi_index < rrss_wi_count; rrss_wi_index++) {

            // Sample incident direction.
            Vec3<Float> wi =
            layered_assembly.randomScatterDirection(path_pcg, outgoing_dir); 

            // Push.
            if (wi[2] > 0) {
                wi_upper.push_back(wi);
            }
            else {
                wi_lower.push_back(wi);
            }
        }

        // Smooth in each hemisphere.
        smoothIncidentDirs(wi_upper);
        smoothIncidentDirs(wi_lower);
        for (Vec3<Float>& wi : wi_lower) {
            wi[2] = -wi[2]; // Flip back to lower hemisphere.
        }

        // Overwrite RRSS incident direction count.
        rrss_wi_count = 
                wi_upper.size() + 
                wi_lower.size();

        // RRSS incident directions.
        rrss_wi = new Vec3<Float>[rrss_wi_count];
        int rrss_wi_index = 0;
        for (const Vec3<Float>& wi : wi_upper) {
            rrss_wi[rrss_wi_index] = wi, rrss_wi_index++;
        }
        for (const Vec3<Float>& wi : wi_lower) {
            rrss_wi[rrss_wi_index] = wi, rrss_wi_index++;
        }
    }

    // RRSS BSDFs.
    Float* rrss_fs = new Float[rrss_wi_count];

    // RRSS BSDF-PDFs.
    Float* rrss_fs_pdf = new Float[rrss_wi_count];

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
                path_pcg, outgoing_dir,
                rrss_wi,
                rrss_wi_count,
                rrss_fs, rrss_fs_pdf);

        // Update progress bar.
        progress_bar.add(
                next_path_count - 
                curr_path_count);
    }

    // RRSS BSDF integral.
    Float rrss_fs_int = 0;

    for (int rrss_wi_index = 0;
             rrss_wi_index < rrss_wi_count; 
             rrss_wi_index++) {

        // BSDF integral term.
        Float fs_int = 
            rrss_fs[rrss_wi_index] / 
            rrss_fs_pdf[rrss_wi_index] /
            pr::fabs(rrss_wi[rrss_wi_index][2]); // Non-cosine-weighted

        // BSDF integral term okay?
        if (pr::isfinite(fs_int)) {

            // Update BSDF integral estimate.
            rrss_fs_int =
            rrss_fs_int + (fs_int - rrss_fs_int) / 
                                   (rrss_wi_index + 1);
        }
    }

    // RRSS samples.
    std::vector<Rrss::Sample> rrss_samples;
    rrss_samples.reserve(rrss_wi_count);

    // BSDF integral okay?
    if (rrss_fs_int > 0 &&
        pr::isfinite(rrss_fs_int)) {

        for (int rrss_wi_index = 0;
                 rrss_wi_index < rrss_wi_count;
                 rrss_wi_index++) {

            // Set target PDF to normalized non-cosine-weighted BSDF.
            Rrss::Sample rrss_sample;
            rrss_sample.dir = rrss_wi[rrss_wi_index];
            rrss_sample.val = rrss_fs[rrss_wi_index];
            rrss_sample.pdf = rrss_fs[rrss_wi_index] /
                     pr::fabs(rrss_wi[rrss_wi_index][2]) / rrss_fs_int;

            // Push sample.
            if (pr::isfinite(rrss_sample.val) &&
                pr::isfinite(rrss_sample.pdf)) {
                rrss_samples.push_back(rrss_sample);
            }
        }
    }
    else {

        for (int rrss_wi_index = 0;
                 rrss_wi_index < rrss_wi_count;
                 rrss_wi_index++) {

            // Set target PDF to BSDF-PDF by default.
            Rrss::Sample rrss_sample;
            rrss_sample.dir = rrss_wi[rrss_wi_index];
            rrss_sample.val = rrss_fs[rrss_wi_index];
            rrss_sample.pdf = rrss_fs_pdf[rrss_wi_index];

            // Push sample.
            if (pr::isfinite(rrss_sample.val) &&
                pr::isfinite(rrss_sample.pdf)) {
                rrss_samples.push_back(rrss_sample);
            }
        }
    }

    // Delete RRSS BSDF-PDFs.
    delete[] rrss_fs_pdf;

    // Delete RRSS BSDFs.
    delete[] rrss_fs;

    // Delete RRSS incident directions.
    delete[] rrss_wi;

    // Redundancy-reduced sample set.
    Rrss rrss(rrss_samples);
    if (!rrss.disableUntil(wi_count)) {
        // Unreachable?
        // TODO Throw
    }

    // Initialize incident directions and BSDF averages.
    incident_dirs.clear();
    bsdf_averages.clear();
    for (const Rrss::Sample& rrss_sample : rrss.samples()) {
        // Is sample enabled?
        if (rrss_sample.is_enabled) {
            incident_dirs.push_back(rrss_sample.dir);
            bsdf_averages.push_back(rrss_sample.val);
        }
    }

    // Initialize path count.
    path_count = rrss_path_count;
}

// Compute BSDF averages.
void FileData::Slice::computeBsdfAverages(
                LayeredAssembly& layered_assembly,
                ProgressBar& progress_bar,
                int path_count_to_add,
                int path_count_to_add_per_iter)
{
    assert(path_count_to_add >= 1 &&
           path_count_to_add_per_iter >= 1);

    for (std::int64_t
            target_path_count = path_count + path_count_to_add;
            target_path_count > path_count;) {

        // Next path count.
        std::int64_t 
            next_path_count = path_count + path_count_to_add_per_iter;
        if (next_path_count > target_path_count) {
            next_path_count = target_path_count;
        }

        // Update BSDF averages.
        layered_assembly.computeAverage(
                path_count, next_path_count,
                path_pcg, 
                outgoing_dir, 
                incident_dirs.data(),
                incident_dirs.size(),
                bsdf_averages.data(), nullptr);

        // Update progress bar.
        progress_bar.add(next_path_count - path_count);

        // Increment path count.
        path_count = next_path_count;
    }
}

} // namespace ls
