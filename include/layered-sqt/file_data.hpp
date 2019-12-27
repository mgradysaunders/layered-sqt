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

public:

    /**
     * @brief Slice.
     */
    struct Slice
    {
    public:

        /**
         * @brief Default constructor.
         */
        Slice() = default;

        /**
         * @brief Clear.
         */
        void clear()
        {
            *this = std::move(Slice());
        }

        /**
         * @brief Read LSQT-slice binary format.
         *
         * @param[inout] istr
         * Input stream.
         *
         * @note
         * This is only intended for use in `FileData::readLss()`.
         */
        void readLss(std::istream& istr);

        /**
         * @brief Write LSQT-slice binary format.
         *
         * @param[inout] ostr
         * Output stream.
         *
         * @note
         * This is only intended for use in `FileData::writeLss()`.
         */
        void writeLss(std::ostream& ostr) const;

        /**
         * @brief Compute incident directions.
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
        void computeIncidentDirs(
                LayeredAssembly& layered_assembly,
                ProgressBar& progress_bar,
                int rrss_oversampling,
                int rrss_path_count,
                int rrss_path_count_per_iter = 512);

        /**
         * @brief Compute BSDF averages.
         *
         * @param[in] layered_assembly
         * Layered assembly.
         *
         * @param[in] progress_bar
         * Progress bar.
         *
         * @param[in] path_count_to_add
         * Path count to add.
         *
         * @param[in] path_count_to_add_per_iter
         * Path count to add per iteration. Set this to a smaller value to 
         * update the progress bar more frequently.
         */
        void computeBsdfAverages(
                LayeredAssembly& layered_assembly, 
                ProgressBar& progress_bar,
                int path_count_to_add,
                int path_count_to_add_per_iter = 4096);

        /**
         * @brief Normalize BSDF to be non-absorbing.
         */
        void normalizeBsdf();

        /**
         * @brief Any transmitted directions?
         */
        bool anyTransmitted() const
        {
            for (const Vec3<Float>& wi : incident_dirs) {
                if (wi[2] < 0) {
                    return true;
                }
            }
            return false;
        }

    public:

        /**
         * @brief Outgoing angle @f$ \theta_o @f$ in radians.
         */
        Float outgoing_angle = 0;
    
        /**
         * @brief Outgoing direction @f$ \omega_o @f$.
         */
        Vec3<Float> outgoing_dir;

        /**
         * @brief Incident directions @f$ \omega_{i[k]} @f$.
         */
        std::vector<Vec3<Float>> incident_dirs;

        /**
         * @brief BSDF averages @f$ f_{[k]} @f$ for each incident direction.
         */
        std::vector<Float> bsdf_averages;

        /**
         * @brief Path count in each BSDF average.
         */
        std::int64_t path_count = 0;

        /**
         * @brief Path PCG.
         */
        Pcg32 path_pcg = {};
    };

public:

    /**
     * @brief Basic initialize.
     */
    void basicInit(
            int wo_count, 
            int wi_count,
            int seed = 0);

    /**
     * @brief Clear.
     */
    void clear()
    {
        slices.clear();
        slices.shrink_to_fit();
    }

    /**
     * @brief Read LSQT-slice binary format.
     *
     * @param[inout] istr
     * Input stream.
     */
    void readLss(std::istream& istr);

    /**
     * @brief Write LSQT-slice binary format.
     *
     * @param[inout] ostr
     * Output stream.
     */
    void writeLss(std::ostream& ostr) const;

    /**
     * @brief Write SQT RAW plain-text format.
     *
     * @param[inout] ostr
     * Output stream.
     */
    void writeSqtRaw(std::ostream& ostr) const;

public:

    /** 
     * @brief Slices.
     */
    std::vector<Slice> slices;
};

/**@}*/

} // namespace ls

#endif // #ifndef LAYERED_SQT_FILE_DATA_HPP

