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
 *
 * This encapsulates an algorithm for building a _redundancy-reduced 
 * sample set_ which, as the name suggests, is a set of sample directions
 * fit to a target solid-angle distribution with minimal redundancy&mdash;that
 * is, minimal unwarranted (and undesirable) clumping of sample directions.
 *
 * @par Algorithm
 *
 * Let @f$ \Omega = \{\omega_k\} @f$, where @f$ k = 1,2,\ldots,N @f$, be a
 * set of @f$ N @f$ random sample directions. We would like to select 
 * @f$ N' < N @f$ of these sample directions to form a subset 
 * @f$ \Omega' \subset \Omega @f$ which fits a solid-angle density @f$ p @f$ 
 * with minimal &ldquo;redundancy&rdquo;. In other words, we intend to discard
 * directions from @f$ \Omega @f$ such that the inherent density of the 
 * remaining directions better represents @f$ p @f$, a target density.
 *
 * We thus require a quantifiable measure of &ldquo;redundancy&rdquo; 
 * so that we can identify the most appropriate sample directions to discard.
 * We define such a measure by inspecting a region @f$ \mathcal{R}_k @f$
 * around each @f$ \omega_k @f$. For each @f$ \omega_k @f$, let 
 * @f$ m_k @f$ be the number of directions we actually find in 
 * @f$ \mathcal{R}_k @f$,
 * and let @f$ \bar{m}_k @f$ be the number of directions we 
 * would expect to find in @f$ \mathcal{R}_k @f$ if @f$ \Omega @f$ 
 * were ideally representative of @f$ p @f$,
 * @f[
 *      m_k = \sum_{k=1}^N \begin{cases}
 *          1 & \omega_k \in    \mathcal{R}_k
 *      \\  0 & \omega_k \notin \mathcal{R}_k
 *      \end{cases}
 * @f]
 * @f[
 *      \bar{m}_k 
 *      = N \operatorname{Pr}(\mathcal{R}_k)
 *      = N \int_{\mathcal{R}_k} p(\omega)\;\mathrm{d}\omega.
 * @f]
 * We define the redundancy @f$ \rho_k @f$ for each 
 * @f$ \omega_k @f$ as 
 * @f[
 *      \rho_k = \frac{m_k}{\bar{m}_k}.
 * @f]
 * When @f$ \rho_k \approx 1 @f$, we see about as many directions 
 * around @f$ \omega_k @f$ as we expect, and so @f$ \omega_k @f$ is 
 * properly representative. When @f$ \rho_k \gg 1 @f$, we see many more 
 * directions around @f$ \omega_k @f$ than we expect, and so @f$ \omega_k @f$ 
 * is redundant. Of course, we must still consider how to form the
 * @f$ \mathcal{R}_k @f$ and how to estimate the @f$ \bar{m}_k @f$. Before
 * continuing however, it is important to discuss how we might apply this 
 * measure to discard directions. To discard only the most redundant direction,
 * we can simply select and discard the direction with the greatest
 * @f$ \rho_k @f$. However, to discard @f$ \ell > 1 @f$ directions, we 
 * cannot simply select and discard the top-@f$ \ell @f$ directions with
 * greatest @f$ \rho_k @f$, as the act of discarding a direction implicitly
 * changes the @f$ \rho_k @f$ of the remaining directions. We must instead
 * discard the most redundant direction @f$ \ell @f$ times, updating the
 * @f$ \rho_k @f$ in between.
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
         * @brief Sample direction.
         */
        Vec3<Float> dir;

        /**
         * @brief Sample value.
         */
        Float val = 0;
        
        /**
         * @brief Target probability density.
         */
        Float pdf = 0;

        /*
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

#if 0
        /**
         * @brief Constructor.
         *
         * @param[in] dir
         * Direction.
         *
         * @param[in] pdf
         * Target probability density.
         */
        Sample(const Vec3<Float>& dir, 
               Float pdf, Float val) : dir(dir), pdf(pdf), val(val)
        {
        }
#endif
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
