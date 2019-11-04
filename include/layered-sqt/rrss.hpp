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
#ifndef LAYERED_SQT_RRSS_HPP
#define LAYERED_SQT_RRSS_HPP

#include <preform/kdtree.hpp>
#include <layered-sqt/common.hpp>

namespace ls {

/**
 * @defgroup rrss Redundancy-reduced sample set
 *
 * `<layered-sqt/rrss.hpp>`
 */
/**@{*/

/**
 * @brief Redundancy-reduced sample set.
 */
class Rrss
{
public:

    /**
     * @brief Sample.
     */
    class Sample
    {
    public:

        /**
         * @brief Direction.
         */
        Vec3<Float> dir;
        
        /**
         * @brief Probability density.
         */
        Float pdf = 0;

        /**
         * @brief Redundancy.
         */
        Float redundancy = 1;

        /**
         * @brief Is enabled?
         */
        bool is_enabled = true;

    public:

        /**
         * @brief Default constructor.
         */
        Sample() = default;

        /**
         * @brief Constructor.
         *
         * @param[in] dir
         * Direction.
         *
         * @param[in] pdf
         * Probability density.
         */
        Sample(const Vec3<Float>& dir, Float pdf) : 
                dir(dir), 
                pdf(pdf)
        {
        }
    };

public:

    /**
     * @brief Constructor.
     */
    Rrss(const std::vector<Sample>& samples);

    /**
     * @brief Non-copyable.
     */
    Rrss(const Rrss&) = delete;

public:

    /**
     * @brief Samples.
     */
    const std::vector<Sample>& samples() const
    {
        return samples_;
    }

    /**
     * @brief Disable `num` redundant samples.
     */
    bool disable(std::size_t num);

    /**
     * @brief Disable redundant samples until `num` remain.
     */
    bool disableUntil(std::size_t num)
    {
        return 
            num >= samples_remaining_ ? 
            num == samples_remaining_ :
            disable(samples_remaining_ - num);
    }

private:

    /**
     * @brief Samples.
     */
    std::vector<Sample> samples_;

    /**
     * @brief Samples remaining (i.e., still enabled).
     */
    std::size_t samples_remaining_ = 0;

    /**
     * @brief Sample tree.
     */
    pr::kdtree<Float, 3, const Sample*> sample_tree_;
};

/**@}*/

} // namespace ls

#endif // #ifndef LAYERED_SQT_RRSS_HPP
