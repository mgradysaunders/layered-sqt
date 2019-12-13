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
#ifndef LAYERED_SQT_FILE_DATA_HPP
#define LAYERED_SQT_FILE_DATA_HPP

#include <layered-sqt/common.hpp>
#include <layered-sqt/medium.hpp>
#include <layered-sqt/layer.hpp>
#include <layered-sqt/layered_assembly.hpp>
#include <layered-sqt/progress_bar.hpp>

namespace ls {

/**
 * @defgroup file_data File data
 *
 * `<layered-sqt/file_data.hpp>`
 */
/**@{*/

/**
 * @brief File data.
 */
class FileData
{
public:

    /**
     * @brief Default constructor.
     */
    FileData() = default;

    /**
     * @brief Destructor.
     */
    ~FileData()
    {
        clear();
    }

public:

    /**
     * @brief Slice.
     */
    struct Slice
    {
    public:

        // TODO Reorganize API?

        /**
         * @brief Initialize.
         */
        void init(Float thetao, int wi_count);

        /**
         * @brief Clear.
         */
        void clear();

        /**
         * @brief Sample incident directions.
         *
         * @param[in] layered_assembly
         * Layered assembly.
         *
         * @param[in] progress_bar
         * Progress bar.
         *
         * @param[in] rrss_oversampling
         * RRSS oversampling multiplier.
         *
         * @param[in] rrss_path_count
         * RRSS path count.
         *
         * @param[in] rrss_path_count_per_iter
         * RRSS path count per iteration. Set this to a smaller value
         * to update the progress bar more frequently.
         */
        void sampleIncidentDirs(
                LayeredAssembly& layered_assembly,
                ProgressBar& progress_bar,
                int rrss_oversampling,
                int rrss_path_count,
                int rrss_path_count_per_iter = 512);

        /**
         * @brief Update BSDF averages.
         *
         * @param[in] layered_assembly
         * Layered assembly.
         *
         * @param[in] progress_bar
         * Progress bar.
         *
         * @param[in] path_count
         * Path count to add.
         *
         * @param[in] path_count_per_iter
         * Path count per iteration. Set this to a smaller value to 
         * update the progress bar more frequently.
         */
        void updateBsdfAverages(
                LayeredAssembly& layered_assembly, 
                ProgressBar& progress_bar,
                int path_count,
                int path_count_per_iter = 4096);

        /**
         * @brief Any transmitted directions?
         */
        bool anyTransmitted() const
        {
            for (int wi_index = 0;
                     wi_index < wi_count_; wi_index++) {
                if (wi_[wi_index][2] < 0) {
                    return true;
                }
            }
            return false;
        }

    public:

        // TODO Non-private variable naming?

        /**
         * @brief Path PCG.
         */
        Pcg32 path_pcg_ = {};

        /**
         * @brief Path count in each BSDF average.
         */
        int path_count_ = 0;

        /**
         * @brief Outgoing angle @f$ \theta_o @f$.
         */
        Float thetao_ = 0;
    
        /**
         * @brief Outgoing direction @f$ \omega_o @f$.
         */
        Vec3<Float> wo_ = {};

        /**
         * @brief Incident direction count.
         */
        int wi_count_ = 0;
    
        /**
         * @brief Incident directions @f$ \omega_{i[k]} @f$.
         */
        Vec3<Float>* wi_ = nullptr;

        /**
         * @brief BSDF averages @f$ f_{[k]} @f$ for each incident direction.
         */
        Float* f_ = nullptr;

        // TODO BSDF average error/convergence estimates?
    };

public:

    // TODO Reorganize API?

    /**
     * @brief Initialize.
     */
    void init(int wo_count, int wi_count);

    /**
     * @brief Clear.
     */
    void clear();

    /**
     * @brief Get slices.
     */
    std::vector<Slice>& slices()
    {
        return slices_;
    }

    /**
     * @brief Get slices.
     */
    const std::vector<Slice>& slices() const
    {
        return slices_;
    }

//  void readLss(std::istream& istr);

//  void writeLss(std::ostream& ostr) const;

    /**
     * @brief Write SQT RAW.
     *
     * @param[inout] ostr
     * Output stream.
     */
    void writeSqtRaw(std::ostream& ostr) const;

private:

    /** 
     * @brief Slices.
     */
    std::vector<Slice> slices_;
};

/**@}*/

} // namespace ls

#endif // #ifndef LAYERED_SQT_FILE_DATA_HPP

