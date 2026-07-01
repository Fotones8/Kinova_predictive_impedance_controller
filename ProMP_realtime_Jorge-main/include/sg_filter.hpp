#pragma once

/**
 * @file sg_filter.hpp
 * @brief Savitzky-Golay (SG) filter for first-derivative (velocity) estimation.
 *
 * Theory:
 *   For each sample i, a polynomial of degree P is fitted in least-squares
 *   sense to W consecutive samples centred on i.  The analytical derivative
 *   of that polynomial at i gives the velocity estimate.  The result is
 *   equivalent to a linear FIR filter (convolution), which is why the method
 *   is efficient and preserves peak amplitudes better than simple finite
 *   differences.
 *
 * Implementation:
 *   - Coefficients are computed per-point via the closed-form Gram matrix
 *     (AtA)^{-1} At, where A is the Vandermonde matrix of the window.
 *   - At boundaries the window is shifted (not shrunk) so the full W points
 *     are always used; the position of point i inside the window changes.
 *   - The Vandermonde matrix is built from the ACTUAL sample timestamps, so the
 *     derivative is exact for non-uniformly-sampled data.  The mean window
 *     spacing is used only as a numerical conditioning scale and does not
 *     affect the result.
 *
 * Reference:
 *   Savitzky A, Golay MJE. Smoothing and Differentiation of Data by
 *   Simplified Least Squares Procedures. Analytical Chemistry. 1964.
 */

#include <algorithm>
#include <cmath>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Cholesky>

namespace promp_rt {

/**
 * @brief Compute SG first-derivative convolution coefficients at a given point.
 *
 * Constructs the Vandermonde matrix A (window_size × poly_order+1) from the
 * ACTUAL time offset of each window sample relative to the evaluation point:
 *   A[i, k] = ((t_i − t_eval) / scale)^k,
 * then computes C = (AᵀA)⁻¹Aᵀ.  Row 1 of C gives the coefficients that, when
 * dotted with the raw signal values, yield dy/dx at the evaluation point in
 * units of 1/scale; dividing by `scale` converts that to per-second.
 *
 * Using the real timestamps (rather than integer sample indices) makes the
 * derivative exact for arbitrarily-spaced data — the previous integer-offset
 * formulation was only correct for uniform sampling.  The `scale` factor is a
 * pure conditioning aid: choosing it ≈ the mean window spacing keeps the
 * Vandermonde entries O(1), and the result is mathematically independent of
 * its value.  For uniform data with scale = dt this reduces exactly to the
 * classic integer-offset Savitzky-Golay coefficients.
 *
 * @param times        Time stamps of all samples (seconds).
 * @param lo           Index of the first window sample in `times`.
 * @param window_size  Number of points in the window (W ≥ poly_order+2).
 * @param eval_idx     Index in `times` of the point where the derivative is
 *                     evaluated (must lie in [lo, lo+window_size-1]).
 * @param poly_order   Polynomial degree (P).  Higher P preserves more signal
 *                     features but amplifies high-frequency noise.
 *                     Typical: 3 or 4.
 * @param scale        Positive conditioning scale (seconds); use the mean
 *                     window spacing.
 *
 * @return  Coefficient vector of length window_size.
 *          velocity ≈ dot(coeffs, signal_window)
 */
inline Eigen::VectorXd sg_deriv_coeffs_at(const std::vector<double>& times,
                                           int lo,
                                           int window_size,
                                           int eval_idx,
                                           int poly_order,
                                           double scale)
{
    // Build Vandermonde matrix from real time offsets: A[i,k] = ((t_i - t_eval)/scale)^k
    const double t_eval = times[static_cast<size_t>(eval_idx)];
    Eigen::MatrixXd A(window_size, poly_order + 1);
    for (int i = 0; i < window_size; ++i) {
        double x  = (times[static_cast<size_t>(lo + i)] - t_eval) / scale;
        double xi = 1.0;
        for (int k = 0; k <= poly_order; ++k) { A(i, k) = xi; xi *= x; }
    }
    // Normal equations: C = (AᵀA)⁻¹Aᵀ  → shape (P+1) × W
    Eigen::MatrixXd C =
        (A.transpose() * A).ldlt().solve(A.transpose());
    // Row 1 = first-derivative coefficients at x=0 (i.e. at t_eval) in units of
    // 1/scale; divide by scale for physical per-second units.
    return C.row(1) / scale;
}

/**
 * @brief Apply Savitzky-Golay filter to estimate the first derivative (velocity).
 *
 * For each sample i a window of exactly `window_size` points is chosen:
 *   - Interior points: window centred on i.
 *   - Boundary points: window shifted so it stays within [0, n-1].
 * This means boundary estimates are slightly one-sided but the window size
 * never decreases, which keeps noise suppression uniform.
 *
 * @param times        Time stamps of each sample (seconds).  May be
 *                     non-uniform; the mean spacing of each local window is
 *                     used for scaling.
 * @param values       Signal values (e.g., joint positions in rad).
 * @param window_size  Number of samples in the filter window (W).
 *                     MUST be odd and ≥ poly_order+2.
 *                     Larger W → more smoothing, more latency.
 *                     Recommended: 7–11 for typical 100–200 Hz motion data.
 *                     Automatically rounded up to nearest odd number if even.
 * @param poly_order   Polynomial degree (P) fitted within each window.
 *                     Higher P → more accurate peaks, less smoothing.
 *                     Recommended: 3 or 4.  Must satisfy P < W.
 *
 * @return  Velocity estimates at each sample (rad/s), same length as `values`.
 *          Returns a zero vector when fewer than 2 samples are provided.
 *
 * Tuning guide:
 *   - Increase window_size to suppress measurement noise at the cost of
 *     slightly blurring velocity peaks.
 *   - Increase poly_order to better preserve sharp acceleration changes,
 *     but do not exceed window_size − 2.
 *   - For 100 Hz data and smooth human motion: window=9, poly=4 is a safe default.
 */
/**
 * @brief Clamp/round a requested (window_size, poly_order) to a valid SG config.
 *
 * Ensures window_size is odd and ≥ poly_order+2, and shrinks it to fit @p n
 * (reducing poly_order to match).  Shared by the single-point and full-array
 * entry points so both use identical windowing.
 */
inline void sg_sanitize_config(int n, int& window_size, int& poly_order)
{
    if (window_size % 2 == 0) ++window_size;
    window_size = std::max(window_size, poly_order + 2);
    if (window_size % 2 == 0) ++window_size;

    if (window_size > n) {
        window_size = (n % 2 == 0) ? n - 1 : n; // largest odd number <= n
        window_size = std::max(window_size, 3);   // need at least 3 for deriv
        poly_order  = std::min(poly_order, window_size - 2);
        poly_order  = std::max(poly_order, 1);    // need at least degree 1
    }
}

/**
 * @brief Savitzky-Golay first derivative (velocity) at a SINGLE sample index.
 *
 * Identical windowing/boundary handling to sg_filter_deriv(), but evaluates only
 * one point — use this when velocities are needed at a sparse set of indices
 * (e.g. one observation per phase step) so cost is O(#queries · W) instead of
 * O(n · W) for the whole buffer.
 *
 * @param times        Time stamps of each sample (s).
 * @param values       Signal values (same length as times).
 * @param i            Index in [0, n-1] at which to evaluate the derivative.
 * @param window_size  SG window length (see sg_filter_deriv).
 * @param poly_order   SG polynomial order (see sg_filter_deriv).
 * @return  Velocity estimate at sample @p i (rad/s); 0 for degenerate input.
 */
inline double sg_deriv_at_index(const std::vector<double>& times,
                                const std::vector<double>& values,
                                int i,
                                int window_size = 9,
                                int poly_order  = 4)
{
    const int n = static_cast<int>(values.size());
    if (n < 2) return 0.0;
    if (n < 3) {
        const double dt = times[1] - times[0];
        return (dt > 0.0) ? (values[1] - values[0]) / dt : 0.0;
    }

    sg_sanitize_config(n, window_size, poly_order);

    const int half = window_size / 2;
    const int lo = std::clamp(i - half, 0, n - window_size);
    const int hi = lo + window_size - 1;
    const int pos_in_win = i - lo;

    // Mean spacing of this window, used only as a conditioning scale for the
    // Vandermonde; the derivative uses the actual timestamps, so it is exact for
    // non-uniform sampling.
    const double dt = (times[static_cast<size_t>(hi)] - times[static_cast<size_t>(lo)])
                      / static_cast<double>(window_size - 1);
    if (dt <= 0.0) return 0.0;

    const Eigen::VectorXd c =
        sg_deriv_coeffs_at(times, lo, window_size, lo + pos_in_win, poly_order, dt);

    double v = 0.0;
    for (int j = 0; j < window_size; ++j)
        v += c(j) * values[static_cast<size_t>(lo + j)];
    return v;
}

inline std::vector<double> sg_filter_deriv(const std::vector<double>& times,
                                            const std::vector<double>& values,
                                            int window_size = 9,
                                            int poly_order  = 4)
{
    const int n = static_cast<int>(values.size());
    if (n < 2) return std::vector<double>(static_cast<size_t>(n), 0.0);
    if (n < 3) {
        const double dt = times[1] - times[0];
        const double v  = (dt > 0.0) ? (values[1] - values[0]) / dt : 0.0;
        return {v, v};
    }

    std::vector<double> result(static_cast<size_t>(n), 0.0);
    for (int i = 0; i < n; ++i)
        result[static_cast<size_t>(i)] =
            sg_deriv_at_index(times, values, i, window_size, poly_order);
    return result;
}

} // namespace promp_rt
