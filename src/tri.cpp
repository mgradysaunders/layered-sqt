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
#include <preform/delaunay.hpp>
#include <preform/image2.hpp>
#include <preform/image_filters.hpp>
#include <preform/float_interval.hpp>
#include <layered-sqt/tri.hpp>

namespace ls {

// Triangulation helper.
void triangulate(
            const std::vector<Vec2<Float>>& locs,
            std::vector<Tri>& tris)
{
    typedef pr::delaunay_triangulation<Float> DelaunayTriangulation;

    // Delaunay triangulation.
    DelaunayTriangulation delaunay;
    delaunay.init(locs.begin(), locs.end());

    // Add triangles.
    tris.clear();
    tris.reserve(delaunay.triangles().size());
    for (const auto& triangle : delaunay.triangles()) {
        Tri tri;
        tri.indices[0] = triangle[0];
        tri.indices[1] = triangle[1];
        tri.indices[2] = triangle[2];
        tris.push_back(tri);
    }

    // Set open triangles.
    for (const auto& edge : delaunay.boundary_edges()) {
        tris[delaunay.boundary_edge_to_triangle(edge)].is_open = true;
    }
}

// Direction to polar location.
Vec2<Float> dirToPolar(const Vec3<Float>& dir)
{
    Float cos_theta = pr::abs(dir[2]);
    Float sin_theta = pr::hypot(dir[0], dir[1]);
    cos_theta = pr::fmin(cos_theta, Float(1));
    sin_theta = pr::fmin(sin_theta, Float(1));
    Float theta = 
        cos_theta < sin_theta
            ? pr::acos(cos_theta)
            : pr::asin(sin_theta);
    Float cos_phi = dir[0] / sin_theta;
    Float sin_phi = dir[1] / sin_theta;
    if (!pr::isfinite(cos_phi) ||
        !pr::isfinite(sin_phi)) {
        cos_phi = 1;
        sin_phi = 0;
    }
    return {
        theta * cos_phi,
        theta * sin_phi
    };
}

// Polar location to direction.
Vec3<Float> polarToDir(const Vec2<Float>& loc)
{
    Float theta = pr::length(loc);
    if (theta > 0) {
        Float sin_theta = pr::sin(theta);
        Float cos_theta = pr::cos(theta);
        return {
            sin_theta * (loc[0] / theta),
            sin_theta * (loc[1] / theta),
            cos_theta
        };
    }
    else {
        return {0, 0, 1};
    }
}

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

    // Initialize vertices.
    vers_.reserve(locs.size());
    auto loc_itr = locs.begin();
    auto val_itr = vals.begin();
    for (; val_itr < vals.end();) {
        vers_.emplace_back(Ver{
            *loc_itr++,
            *val_itr++
        });
    }

    // Initialize triangles.
    triangulate(locs, tris_);
}

// Value.
Float TriInterpolator::value(const Vec2<Float>& loc) const
{
    // Interpolate.
    for (const Tri& tri : tris_) {
        const Ver& ver0 = vers_[tri.indices[0]];
        const Ver& ver1 = vers_[tri.indices[1]];
        const Ver& ver2 = vers_[tri.indices[2]];

        // Barycentric coordinates.
        Vec2<float> h0 = ver0.loc - loc;
        Vec2<float> h1 = ver1.loc - loc;
        Vec2<float> h2 = ver2.loc - loc;
        float b0 = h1[0] * h2[1] - h1[1] * h2[0];
        float b1 = h2[0] * h0[1] - h2[1] * h0[0];
        float b2 = h0[0] * h1[1] - h0[1] * h1[0];
        if (b0 == 0.f ||
            b1 == 0.f ||
            b2 == 0.f) {
            b0 = double(h1[0]) * double(h2[1]) - double(h1[1]) * double(h2[0]);
            b1 = double(h2[0]) * double(h0[1]) - double(h2[1]) * double(h0[0]);
            b2 = double(h0[0]) * double(h1[1]) - double(h0[1]) * double(h1[0]);
        }
        if (b0 < 0.f || b1 < 0.f || (!tri.is_open && b2 < 0.f)) {
            continue;
        }

        // Normalize.
        float q = b0 + b1 + b2;
        if (q == 0.f) {
            continue;
        }
        b0 /= q;
        b1 /= q;
        b2 /= q;

        // Is inside triangle?
        if (b2 > 0.f) {
            Float val = 
                Float(b0) * ver0.val + 
                Float(b1) * ver1.val + 
                Float(b2) * ver2.val;
           return val;
        }
        else {
            // Segment locations.
            const Vec2<Float>& seg_loc0 = ver0.loc;
            const Vec2<Float>& seg_loc1 = ver1.loc;
            Vec2<Float> seg_vec = seg_loc1 - seg_loc0;

            // Find nearest location on segment.
            Float u =
                pr::dot(seg_vec, loc - seg_loc0) / 
                pr::dot(seg_vec, seg_vec);
            u = pr::fmin(u, Float(1));
            u = pr::fmax(u, Float(0));
            Float val = (1 - u) * ver0.val + u * ver1.val; 
            return val;
        }
    }

    return 0;
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
            return lhs.outgoing_angle_ < 
                   rhs.outgoing_angle_;
        });
}

// Value.
Float TriFileData::value(
                const Vec3<Float>& tmp_wo, 
                const Vec3<Float>& tmp_wi_world) const
{
    Vec3<Float> wo = tmp_wo;
    Vec3<Float> wi_world = tmp_wi_world;
    if (wo[2] < 0) {
        wo[2] = -wo[2];
        wi_world[2] = -wi_world[2];
    }
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
    Float thetao = 
        pr::fabs(cos_thetao) < sin_thetao ?
        pr::acos(cos_thetao) :
        pr::asin(sin_thetao);

    // Local direction.
    Vec3<Float> wi = {
         cos_phio * wi_world[0] + sin_phio * wi_world[1],
        -sin_phio * wi_world[0] + cos_phio * wi_world[1],
         wi_world[2]
    };

    // Apply warping.
    Vec2<Float> wi_loc = dirToPolar(wi);

    auto slice_itr = 
    std::lower_bound(
            slices_.begin(),
            slices_.end(),
            thetao,
            [=](const Slice& lhs, Float rhs) {
                return lhs.outgoing_angle_ < rhs;
            });
    if (slice_itr == slices_.begin() ||
        slice_itr == slices_.end()) {
        if (slice_itr == slices_.end()) {
            slice_itr--;
        }
        return slice_itr->value(wi_loc, wi[2]);
    }
    else {
        const Slice& slice1 = *slice_itr--;
        const Slice& slice0 = *slice_itr;
        Float val0 = slice0.value(wi_loc, wi[2]);
        Float val1 = slice1.value(wi_loc, wi[2]);
        Float thetao0 = slice0.outgoing_angle_;
        Float thetao1 = slice1.outgoing_angle_;
        Float fac = (thetao - thetao0) / (thetao1 - thetao0);
        return (1 - fac) * val0 + fac * val1;
    }
}

// Initialize.
void TriFileData::Slice::init(const FileData::Slice& file_data_slice)
{
    // Outgoing angle.
    outgoing_angle_ = file_data_slice.outgoing_angle;

    // Partition into upper/lower.
    std::vector<Vec2<Float>> locs_upper;
    std::vector<Vec2<Float>> locs_lower;
    std::vector<Float> vals_upper;
    std::vector<Float> vals_lower;

    // Iterate.
    auto incident_dir = file_data_slice.incident_dirs.begin();
    auto bsdf_average = file_data_slice.bsdf_averages.begin();
    for (; bsdf_average < file_data_slice.bsdf_averages.end();) {

        Vec3<Float> wi = *incident_dir++;
        Float wi_val = *bsdf_average++ / pr::fabs(wi[2]);
        if (!pr::isfinite(wi_val)) {
            continue;
        }

        // Warping.
        Vec2<Float> wi_loc = dirToPolar(wi);

        // In upper hemisphere?
        if (!pr::signbit(wi[2])) {
            // Add to upper hemisphere.
            locs_upper.push_back(wi_loc);
            vals_upper.push_back(wi_val);
        }
        else {
            // Add to lower hemisphere.
            locs_lower.push_back(wi_loc);
            vals_lower.push_back(wi_val);
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
