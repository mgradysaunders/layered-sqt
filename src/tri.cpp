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
#include <preform/delaunay.hpp>
#include <layered-sqt/tri.hpp>

namespace ls {

// Delaunay triangulation helper.
typedef pr::delaunay_triangulation<Float> DelaunayTriangulation;

// Initialize.
void TriInterpolator::init(
            const std::vector<Vec2<Float>>& locs,
            const std::vector<Float>& vals)
{
    if (locs.size() != vals.size() ||
        locs.size() <= 2) {
        if (locs.empty() &&
            vals.empty()) {
            return;
        }
        else {
            throw std::invalid_argument(__PRETTY_FUNCTION__);
        }
    }

    DelaunayTriangulation delaunay;
    delaunay.init(locs.begin(), locs.end());
    tris_.reserve(delaunay.triangles().size());
    for (const auto& index_tri : delaunay.triangles()) {
        Tri tri;
        for (int k = 0; k < 3; k++) {
            tri.vertices[k].loc = locs[index_tri[k]];
            tri.vertices[k].val = vals[index_tri[k]];
        }
        tris_.push_back(tri);
    }
}

// Value.
std::optional<Float> TriInterpolator::Tri::value(const Vec2<Float>& loc) const
{
    Vec2<float> h0 = vertices[0].loc - loc;
    Vec2<float> h1 = vertices[1].loc - loc;
    Vec2<float> h2 = vertices[2].loc - loc;
    float b0 = h1[1] * h2[0] - h1[0] * h2[1];
    float b1 = h2[1] * h0[0] - h2[0] * h0[1];
    float b2 = h0[1] * h1[0] - h0[0] * h1[1];
    if (b0 == 0.f ||
        b1 == 0.f ||
        b2 == 0.f) {
        b0 = double(h1[1]) * double(h2[0]) - double(h1[0]) * double(h2[1]);
        b1 = double(h2[1]) * double(h0[0]) - double(h2[0]) * double(h0[1]);
        b2 = double(h0[1]) * double(h1[0]) - double(h0[0]) * double(h1[1]);
    }
    if ((b0 < 0.f || b1 < 0.f || b2 < 0.f) &&
        (b0 > 0.f || b1 > 0.f || b2 > 0.f)) {
        return std::nullopt;
    }

    float q = b0 + b1 + b2;
    if (q == 0.f) {
        return std::nullopt;
    }
    b0 /= q;
    b1 /= q;
    b2 /= q;
    Float val = 
        Float(b0) * vertices[0].val + 
        Float(b1) * vertices[1].val + 
        Float(b2) * vertices[2].val;
    return std::make_optional(val);
}

// Value.
std::optional<Float> TriInterpolator::value(const Vec2<Float>& loc) const
{
    std::optional<Float> val = std::nullopt;
    for (Tri tri : tris_) {
        if ((val = tri.value(loc)).has_value()) {
            break;
        }
    }
    return val;
}

// Initialize.
void TriFileData::init(const FileData& file_data)
{
    // Initialize slices.
    slices_.clear();
    for (const FileData::Slice& file_data_slice : file_data.slices) {
        slices_.emplace_back();
        slices_.back().init(file_data_slice);
    }

    // Sort by outgoing angle.
    slices_.sort(
        [=](const Slice& lhs, const Slice& rhs) {
            return lhs.outgoing_dirz_ < 
                   rhs.outgoing_dirz_;
        });
}

// Value.
Float TriFileData::value(
                const Vec3<Float>& wo, 
                const Vec3<Float>& wi) const
{
    Float cos_thetao = wo[2];
    Float cos2_thetao = cos_thetao * cos_thetao;
    cos2_thetao = pr::fmin(cos2_thetao, Float(+1));
    cos2_thetao = pr::fmax(cos2_thetao, Float(-1));
    Float sin2_thetao = 1 - cos2_thetao;
    Float sin_thetao = pr::sqrt(sin2_thetao);
    Float cos_phio = wo[0] / sin_thetao;
    Float sin_phio = wo[1] / sin_thetao;
    if (!(pr::isfinite(cos_phio) && 
          pr::isfinite(sin_phio))) {
        cos_phio = 1;
        sin_phio = 0;
    }

    // Local directions.
/*  Vec3<Float> wo_local = {sin_thetao, 0, cos_thetao}; */
    Vec3<Float> wi_local = {
         cos_phio * wi[0] + sin_phio * wi[1],
        -sin_phio * wi[0] + cos_phio * wi[1],
        wi[2]
    };

    auto slice_itr = 
    std::lower_bound(
            slices_.begin(),
            slices_.end(),
            cos_thetao,
            [=](const Slice& lhs, Float rhs) {
                return lhs.outgoing_dirz_ < rhs;
            });
    if (slice_itr == slices_.begin() ||
        slice_itr == slices_.end()) {
        if (slice_itr == slices_.end()) {
            slice_itr--;
        }
        return slice_itr->value(wi_local).value_or(0);
    }
    else {
        const Slice& slice1 = *slice_itr--;
        const Slice& slice0 = *slice_itr;
        Float val0 = slice0.value(wi_local).value_or(0);
        Float val1 = slice1.value(wi_local).value_or(0);
        if (val0 == 0) return val1;
        if (val1 == 0) return val0;
        Float cos_thetao0 = slice0.outgoing_dirz_;
        Float cos_thetao1 = slice1.outgoing_dirz_;
        Float fac = (cos_thetao - cos_thetao0) / (cos_thetao1 - cos_thetao0);
        return (1 - fac) * val0 + fac * val1;
    }
}

// Initialize.
void TriFileData::Slice::init(const FileData::Slice& file_data_slice)
{
    // Outgoing direction Z-component.
    outgoing_dirz_ = file_data_slice.outgoing_dir[2];

    // Partition into upper/lower.
    std::vector<Vec2<Float>> locs_upper;
    std::vector<Vec2<Float>> locs_lower;
    std::vector<Float> vals_upper;
    std::vector<Float> vals_lower;

    // Iterate.
    auto incident_dir = file_data_slice.incident_dirs.begin();
    auto bsdf_average = file_data_slice.bsdf_averages.begin();
    for (; bsdf_average < file_data_slice.bsdf_averages.end();) {

        // In upper hemisphere?
        if ((*incident_dir)[2] > 0) {
            // Add to upper hemisphere.
            locs_upper.push_back(*incident_dir);
            vals_upper.push_back(*bsdf_average);
        }
        else {
            // Add to lower hemisphere.
            locs_lower.push_back(*incident_dir);
            vals_lower.push_back(*bsdf_average);
        }

        // Increment.
        incident_dir++;
        bsdf_average++;
    }

    // Initialize.
    bsdf_upper_.init(locs_upper, vals_upper);
    bsdf_lower_.init(locs_lower, vals_lower);
}

} // namespace ls
