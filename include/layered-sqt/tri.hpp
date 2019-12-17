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
#include <layered-sqt/common.hpp>

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

/**@}*/

} // namespace ls

#endif // #ifndef LAYERED_SQT_TRI_HPP
