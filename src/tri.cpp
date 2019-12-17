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
        throw std::invalid_argument(__PRETTY_FUNCTION__);
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

} // namespace ls
