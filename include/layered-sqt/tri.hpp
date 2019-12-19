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
#pragma once
#ifndef LAYERED_SQT_TRI_HPP
#define LAYERED_SQT_TRI_HPP

#include <optional>
#include <list>
#include <layered-sqt/common.hpp>
#include <layered-sqt/file_data.hpp>

namespace ls {

/**
 * @defgroup tri Triangulation
 *
 * `<layered-sqt/tri.hpp>`
 */
/**@{*/

/**
 * @brief Triangle interpolator.
 */
class TriInterpolator
{
public:

    /**
     * @brief Triangle.
     */
    class Tri
    {
    public:

        /**
         * @brief Vertex.
         */
        class Vertex 
        {
        public:

            /**
             * @brief Location.
             */
            Vec2<Float> loc;

            /**
             * @brief Value.
             */
            Float val = 0;
        };

        /**
         * @brief Vertices.
         */
        Vertex vertices[3];

        /**
         * @brief Interpolate value.
         *
         * @param[in] loc
         * Location.
         */
        std::optional<Float> value(const Vec2<Float>& loc) const;
    };

    /**
     * @brief Default constructor.
     */
    TriInterpolator() = default;

    /**
     * @brief Initialize.
     *
     * @param[in] locs
     * Locations.
     *
     * @param[in] vals
     * Values.
     */
    void init(const std::vector<Vec2<Float>>& locs,
              const std::vector<Float>& vals);

    /**
     * @brief Initialize directly.
     *
     * @param[in] tris
     * Triangles.
     */
    void init(const std::vector<Tri>& tris)
    {
        tris_ = tris;
    }

    /**
     * @brief Clear.
     */
    void clear()
    {
        tris_.clear();
    }

    /**
     * @brief Interpolate value.
     *
     * @param[in] loc
     * Location.
     */
    std::optional<Float> value(const Vec2<Float>& loc) const;

private:

    /**
     * @brief Triangles.
     */
    std::vector<Tri> tris_;
};

/**
 * @brief Triangulated file data.
 */
class TriFileData
{
public:

    /**
     * @brief Initialize.
     *
     * @param[in] file_data
     * File data.
     */
    void init(const FileData& file_data);

    /**
     * @brief Value.
     *
     * @param[in] wo
     * Outgoing direction @f$ \omega_o @f$.
     *
     * @param[in] wi
     * Incident direction @f$ \omega_i @f$.
     */
    Float value(
            const Vec3<Float>& wo, 
            const Vec3<Float>& wi) const;

    /**
     * @brief Triangulated file data slice.
     */
    class Slice
    {
    public:

        /**
         * @brief Initialize.
         *
         * @param[in] file_data_slice
         * File data slice.
         */
        void init(const FileData::Slice& file_data_slice);

        /**
         * @brief Value.
         *
         * @param[in] dir
         * Direction.
         */
        std::optional<Float> value(const Vec3<Float>& dir) const
        {
            return dir[2] > 0 ?
                    bsdf_upper_.value(Vec2<Float>(dir)) :
                    bsdf_lower_.value(Vec2<Float>(dir));
        }

    private:

        /**
         * @brief Outgoing direction Z-component.
         */
        Float outgoing_dirz_ = 0;

        /**
         * @brief Triangulated BSDF upper hemisphere.
         */
        TriInterpolator bsdf_upper_;

        /**
         * @brief Triangulated BSDF lower hemisphere.
         */
        TriInterpolator bsdf_lower_;

        // Friend.
        friend class TriFileData;
    };


private:

    /**
     * @brief Triangulated file data slices.
     */
    std::list<Slice> slices_;
};

/**@}*/

} // namespace ls

#endif // #ifndef LAYERED_SQT_TRI_HPP
