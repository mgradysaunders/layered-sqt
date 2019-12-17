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
#include <random>
#include <preform/byte_order.hpp>
#include <layered-sqt/rrss.hpp>
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
            pr::numeric_constants<Float>::M_pi_2();

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
void FileData::writeSqtRaw(std::ostream& ostr) const
{
    bool is_bsdf = false;
    for (const Slice& slice : slices) {
        if (slice.anyTransmitted()) {
            is_bsdf = true;
            break;
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
            ostr << wo_index << ' ';
            ostr << slice.incident_dirs[wi_index][0] << ' ';
            ostr << slice.incident_dirs[wi_index][1] << ' ';
            ostr << slice.incident_dirs[wi_index][2] << ' ';
            ostr << slice.bsdf_averages[wi_index] /
                    pr::fabs(slice.incident_dirs[wi_index][2]) << '\n';
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
    std::uint64_t path_pcg_mem[2];
    istr_le >> path_pcg_mem;
    std::memcpy(
            &path_pcg, 
            &path_pcg_mem[0], sizeof(path_pcg_mem));
    static_assert(sizeof(path_pcg) == sizeof(path_pcg_mem));
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
    std::uint64_t path_pcg_mem[2];
    std::memcpy(
            &path_pcg_mem[0], 
            &path_pcg, sizeof(path_pcg));
    static_assert(sizeof(path_pcg) == sizeof(path_pcg_mem));
    ostr_le << path_pcg_mem;
}

// Compute incident directions.
void FileData::Slice::computeIncidentDirs(
                LayeredAssembly& layered_assembly,
                ProgressBar& progress_bar,
                int rrss_oversampling,
                int rrss_path_count, 
                int rrss_path_count_per_iter)
{
    assert(rrss_oversampling >= 1 &&
           rrss_path_count >= 1 &&
           rrss_path_count_per_iter >= 1);
    if (incident_dirs.empty()) {
        return;
    }

    // Incident direction count.
    int wi_count = incident_dirs.size();

    // Incident directions.
    Vec3<Float>* wi = incident_dirs.data();

    // RRSS incident direction count.
    int rrss_wi_count = 
        rrss_oversampling * wi_count;

    // RRSS incident directions.
    Vec3<Float>* rrss_wi = new Vec3<Float>[rrss_wi_count];

    // Initialize incident directions.
    for (int rrss_wi_index = 0;
             rrss_wi_index < rrss_wi_count; rrss_wi_index++) {
    
        rrss_wi[rrss_wi_index] = 
        layered_assembly.randomScatterDirection(path_pcg, outgoing_dir); 
    }

    // BSDFs.
    Float* f = bsdf_averages.data();

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
                path_pcg, outgoing_dir,
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
    if (!rrss.disableUntil(wi_count)) {
        // Unreachable?
    }

    // Initialize incident directions and BSDF averages.
    int wi_index = 0;
    for (int rrss_wi_index = 0;
             rrss_wi_index < rrss_wi_count; rrss_wi_index++) {

        // Is sample enabled?
        if (rrss.samples()[rrss_wi_index].is_enabled) {

            // Use incident direction and BSDF average.
            wi[wi_index] = rrss_wi[rrss_wi_index];
            f [wi_index] = rrss_f [rrss_wi_index]; wi_index++;
        }
    }

    // Initialize path count.
    path_count = rrss_path_count;

    // Delete RRSS BSDF-PDFs.
    delete[] rrss_f_pdf;

    // Delete RRSS BSDFs.
    delete[] rrss_f;

    // Delete RRSS incident directions.
    delete[] rrss_wi;
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