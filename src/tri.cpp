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
#include <preform/image2.hpp>
#include <preform/image_filters.hpp>
#include <preform/float_interval.hpp>
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

    // Delaunay triangulation.
    DelaunayTriangulation delaunay;
    delaunay.init(locs.begin(), locs.end());

    // Add triangles.
    tris_.reserve(delaunay.triangles().size());
    for (const auto& index_tri : delaunay.triangles()) {
        Tri tri;
        for (int k = 0; k < 3; k++) {
            tri.vertices[k].loc = locs[index_tri[k]];
            tri.vertices[k].val = vals[index_tri[k]];
        }
        tris_.push_back(tri);
    }

    // Add triangle edges on convex hull.
    tri_edges_.reserve(delaunay.boundary_edges().size());
    for (const auto& boundary_edge : delaunay.boundary_edges()) {
        TriEdge tri_edge;
        tri_edge.vertices[0].loc = locs[boundary_edge.a];
        tri_edge.vertices[0].val = vals[boundary_edge.a];
        tri_edge.vertices[1].loc = locs[boundary_edge.b];
        tri_edge.vertices[1].val = vals[boundary_edge.b];
        tri_edges_.push_back(tri_edge);
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

// Interpolate value.
void TriInterpolator::TriEdge::value(
                const Vec2<Float>& loc, 
                Float& val,
                Float& min_dist2) const
{
    // Segment locations.
    const Vec2<Float>& seg_loc0 = vertices[0].loc;
    const Vec2<Float>& seg_loc1 = vertices[1].loc;
    Vec2<Float> seg_vec = seg_loc1 - seg_loc0;

    // Find nearest location on segment.
    Float u = 
        pr::dot(seg_vec, loc - seg_loc0) / 
        pr::dot(seg_vec, seg_vec);
    u = pr::fmin(u, Float(1));
    u = pr::fmax(u, Float(0));
    Vec2<Float> seg_loc = (1 - u) * seg_loc0 + u * seg_loc1;
    Float cur_dist2 = pr::dot(loc - seg_loc, loc - seg_loc);
    if (!(cur_dist2 < min_dist2)) {
        return;
    }
    else {
        // Update minimum distance square.
        min_dist2 = cur_dist2;

        // Update value.
        const Float& seg_val0 = vertices[0].val;
        const Float& seg_val1 = vertices[1].val;
        val = (1 - u) * seg_val0 + u * seg_val1;
    }
}

// Value.
Float TriInterpolator::value(const Vec2<Float>& loc) const
{
    // Interpolate.
    for (Tri tri : tris_) {
        std::optional<Float> val = std::nullopt;
        if ((val = tri.value(loc)).has_value()) {
            return val.value();
        }
    }

    // Location is outside convex hull.
    Float val = 0;
    Float min_dist2 = pr::numeric_limits<Float>::infinity();
    for (TriEdge tri_edge : tri_edges_) {
        tri_edge.value(loc, val, min_dist2);
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
    cos_thetao = pr::fmin(cos_thetao, Float(+1));
    cos_thetao = pr::fmax(cos_thetao, Float(-1));
    Float sin_thetao = pr::hypot(wo[0], wo[1]);
    Float cos_phio = wo[0] / sin_thetao;
    Float sin_phio = wo[1] / sin_thetao;
    if (sin_thetao < 0.000001) {
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
        return slice_itr->value(wi_local) * pr::fabs(wi[2]);
    }
    else {
        const Slice& slice1 = *slice_itr--;
        const Slice& slice0 = *slice_itr;
        Float val0 = slice0.value(wi_local);
        Float val1 = slice1.value(wi_local);
        Float cos_thetao0 = slice0.outgoing_dirz_;
        Float cos_thetao1 = slice1.outgoing_dirz_;
        Float fac = (cos_thetao - cos_thetao0) / (cos_thetao1 - cos_thetao0);
        return ((1 - fac) * val0 + fac * val1) * pr::fabs(wi[2]);
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

        Vec3<Float> loc = *incident_dir++;
        Float val = *bsdf_average++ / pr::fabs(loc[2]);
        if (!pr::isfinite(val)) {
            continue;
        }

        // In upper hemisphere?
        if (!pr::signbit(loc[2])) {
            // Add to upper hemisphere.
            locs_upper.push_back(loc);
            vals_upper.push_back(val);
        }
        else {
            // Add to lower hemisphere.
            locs_lower.push_back(loc);
            vals_lower.push_back(val);
        }
    }

    // Initialize.
    bsdf_upper_.init(locs_upper, vals_upper);
    bsdf_lower_.init(locs_lower, vals_lower);
}

// Intersect unit sphere.
static 
bool intersectUnitSphere(Ray ray, Hit& hit)
{
    // Float interval.
    typedef pr::float_interval<Float> FloatInterval;
    FloatInterval t0;
    FloatInterval t1;
    FloatInterval::solve_poly2(
            pr::dot(ray.pos, ray.pos) - 1,
            pr::dot(ray.dir, ray.pos) * 2,
            pr::dot(ray.dir, ray.dir),
            t0, t1);
    constexpr Float tmin = 0;
    constexpr Float tmax = pr::numeric_limits<Float>::infinity();
    if (!(t0.upper_bound() < tmax &&
          t1.lower_bound() > tmin)) {
        return false;
    }

    // Select root.
    FloatInterval t = t0;
    if (!(t.upper_bound() < tmax &&
          t.lower_bound() > tmin)) {
        t = t1;
        if (!(t.upper_bound() < tmax &&
              t.lower_bound() > tmin)) {
            return false;
        }
    }

    // Initialize hit.
    hit.pos = pr::normalize_fast(ray.pos + ray.dir * t.value());
    return true;
}

// Render sphere example.
void TriFileData::renderSphereExample(int image_dim, float* image_pixels) const
{
    if (!(image_dim > 0 && 
          image_pixels)) {
        throw std::invalid_argument(__PRETTY_FUNCTION__);
    }

    // 2-dimensional image with 1 channel.
    typedef pr::image2<Float, float, 1> Image;

    // 2-dimensional Mitchell filter.
    typedef pr::mitchell_filter2<Float> ImageFilter;

    // Initialize image.
    Image image;
    image.resize(Vec2<int>{image_dim, image_dim});

    // Initialize image filter.
    ImageFilter image_filter;
    Vec2<Float> image_filter_rad = {1, 1};

    // Initialize light directions.
    Vec3<Float> l0 = pr::normalize(Vec3<Float>{-5, +2, -4});
    Vec3<Float> l1 = pr::normalize(Vec3<Float>{+1, +2, +2});
    Vec3<Float> l2 = pr::normalize(Vec3<Float>{+2, -3, 0});

    // Iterate pixels.
    for (int i = 0; i < image_dim; i++)
    for (int j = 0; j < image_dim; j++) {

        // Iterate sub-pixel samples.
        for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++) {

            // Image location.
            Vec2<Float> image_loc = {
                Float(i) + Float(k) / 3,
                Float(j) + Float(l) / 3
            };

            // Ray end points.
            Vec3<Float> pos0 = {0, 0, -4};
            Vec3<Float> pos1 = {
                -3 * (image_loc[1] / image_dim - Float(0.5)),
                -3 * (image_loc[0] / image_dim - Float(0.5)),
                0
            };
            Ray ray;
            ray.pos = pos0;
            ray.dir = pr::normalize(pos1 - pos0);

            // Intersect unit sphere.
            Hit hit;
            if (intersectUnitSphere(ray, hit)) {

                // Reconstruct image.
                Mat3<Float> tbn = 
                Mat3<Float>::build_onb(hit.pos);
                Vec3<Float> wo = pr::dot(pr::transpose(tbn), -ray.dir);
                Vec3<Float> wi0 = pr::dot(pr::transpose(tbn), l0);
                Vec3<Float> wi1 = pr::dot(pr::transpose(tbn), l1);
                Vec3<Float> wi2 = pr::dot(pr::transpose(tbn), l2);
                Float f0 = value(wo, wi0) * 2;
                Float f1 = value(wo, wi1);
                Float f2 = value(wo, wi2) * Float(0.08);
                image.reconstruct(
                        {(f0 + f1 + f2) / 9},
                        image_loc,
                        image_filter_rad,
                        image_filter);
            }
        }
    }

    // Write image pixels.
    for (int i = 0; i < image_dim; i++)
    for (int j = 0; j < image_dim; j++) {
        *image_pixels++ = image(i, j)[0];
    }
}

} // namespace ls
