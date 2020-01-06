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
#pragma once
#ifndef LAYERED_SQT_PROGRESS_BAR_HPP
#define LAYERED_SQT_PROGRESS_BAR_HPP

#include <layered-sqt/common.hpp>

namespace ls {

/**
 * @defgroup progress_bar Progress bar
 *
 * `<layered-sqt/progress_bar.hpp>`
 */
/**@{*/

/**
 * @brief Progress bar.
 */
class ProgressBar
{
public:

    /**
     * @brief Constructor.
     */
    ProgressBar(std::int64_t path_count_total) :
                             path_count_total_(path_count_total)
    {
    }

    /**
     * @brief Add paths.
     */
    void add(int path_count);

    /**
     * @brief Finish, print newlines.
     */
    void finish();

private:

    /**
     * @brief Path count.
     */
    std::int64_t path_count_ = 0;

    /**
     * @brief Path count total.
     */
    std::int64_t path_count_total_ = 0;
};

/**@}*/

} // namespace ls

#endif // #ifndef LAYERED_SQT_PROGRESS_BAR_HPP
