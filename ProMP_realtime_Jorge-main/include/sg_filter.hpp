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
 *   - Time spacing inside each window is assumed uniform; the mean spacing of
 *     the window is used, which is exact for uniformly-sampled data and a
 *     good approximation otherwise.
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
 * @brief Compute SG first-derivative convolution coefficients.
 *
 * Constructs the Vandermonde matrix A (window_size × poly_order+1) where
 *   A[i, k] = (i − center)^k,
 * then computes C = (AᵀA)⁻¹Aᵀ.  Row 1 of C gives the coefficients that,
 * when dotted with the raw signal values in the window, yield dy/dx at
 * position `center`.  Dividing by `dt` converts from per-sample to per-second.
 *
 * @param window_size  Number of points in the window (W ≥ poly_order+2).
 * @param poly_order   Polynomial degree (P).  Higher P preserves more signal
 *                     features but amplifies high-frequency noise.
 *                     Typical: 3 or 4.
 * @param center       Index of the evaluation point within the window [0, W-1].
 * @param dt           Mean time spacing of the window (seconds).
 *
 * @return  Coefficient vector of length window_size.
 *          velocity ≈ dot(coeffs, signal_window)
 */
inline Eigen::VectorXd sg_deriv_coeffs_at(int window_size,
                                           int poly_order,
                                           int center,
                                           double dt)
{
    // Build Vandermonde matrix: A[i,k] = (i - center)^k
    Eigen::MatrixXd A(window_size, poly_order + 1);
    for (int i = 0; i < window_size; ++i) {
        double x  = static_cast<double>(i - center);
        double xi = 1.0;
        for (int k = 0; k <= poly_order; ++k) { A(i, k) = xi; xi *= x; }
    }
    // Normal equations: C = (AᵀA)⁻¹Aᵀ  → shape (P+1) × W
    Eigen::MatrixXd C =
        (A.transpose() * A).ldlt().solve(A.transpose());
    // Row 1 = first-derivative coefficients at x=0; divide by dt for physical units
    return C.row(1) / dt;
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
inline std::vector<double> sg_filter_deriv(const std::vector<double>& times,
                                            const std::vector<double>& values,
                                            int window_size = 9,
                                            int poly_order  = 4)
{
    int n = static_cast<int>(values.size());
    if (n < 2) return std::vector<double>(static_cast<size_t>(n), 0.0);

    // --- Special case: too few points for SG → simple finite differences ---
    if (n < 3) {
        double dt = times[1] - times[0];
        double v  = (dt > 0.0) ? (values[1] - values[0]) / dt : 0.0;
        return {v, v};
    }

    // Sanitise window: must be odd, >= poly_order+2.
    if (window_size % 2 == 0) ++window_size;
    window_size = std::max(window_size, poly_order + 2);
    if (window_size % 2 == 0) ++window_size;

    // Shrink to fit n (keep odd, keep >= 3, reduce poly_order accordingly).
    if (window_size > n) {
        window_size = (n % 2 == 0) ? n - 1 : n; // largest odd number <= n
        window_size = std::max(window_size, 3);   // need at least 3 for deriv
        poly_order  = std::min(poly_order, window_size - 2);
        poly_order  = std::max(poly_order, 1);    // need at least degree 1
    }

    int half = window_size / 2;
    std::vector<double> result(static_cast<size_t>(n), 0.0);

    for (int i = 0; i < n; ++i) {
        // Shift window to stay within bounds (boundary handling)
        int lo = std::clamp(i - half, 0, n - window_size);
        int hi = lo + window_size - 1;
        int pos_in_win = i - lo; // where point i sits in [0, window_size-1]

        // Mean dt for this window
        double dt = (times[static_cast<size_t>(hi)] - times[static_cast<size_t>(lo)])
                    / static_cast<double>(window_size - 1);
        if (dt <= 0.0) { result[static_cast<size_t>(i)] = 0.0; continue; }

        Eigen::VectorXd c = sg_deriv_coeffs_at(window_size, poly_order, pos_in_win, dt);

        double v = 0.0;
        for (int j = 0; j < window_size; ++j)
            v += c(j) * values[static_cast<size_t>(lo + j)];
        result[static_cast<size_t>(i)] = v;
    }
    return result;
}

} // namespace promp_rt
