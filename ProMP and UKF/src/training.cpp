/**
 * @file training.cpp
 * @brief ProMP_1 offline training from n_dof–DoF CSV demonstrations.
 *
 * ── Processing pipeline ──────────────────────────────────────────────────────
 *   1. Scan demo_folder for *.csv files (one file = one demonstration).
 *   2. Parse (n_dof + 1) columns: [time_s, pos0, …, pos(n_dof-1)].
 *      Optionally clip rows beyond max_demo_duration.
 *   3. Compute per-joint velocities from the RAW (pre-resample) position
 *      samples using the Savitzky-Golay filter → velocity quality scales
 *      with the original data frequency.
 *   4. Resample all 2*n_dof channels to a uniform n_resample_steps grid via
 *      linear interpolation on normalised phase.
 *   5. Train promp::ProMP (dims = 2*n_dof) via ridge regression + Gaussian
 *      weight distribution fitting.
 *   6. Persist model and metadata.
 *
 * ── Data-frequency handling ──────────────────────────────────────────────────
 * Different demos may have DIFFERENT sampling rates or irregular timestamps:
 *   • SG filter uses actual time stamps → correct for any rate.
 *   • Interpolation normalises time to [0,1] before resampling → rate-agnostic.
 *   • SG window (config.sg_window) should be set relative to the original
 *     data rate (see TrainingConfig::sg_window documentation).
 *
 * ── Dimension layout in the trained ProMP (dims = 2*n_dof) ──────────────────
 *   col 0: pos0 (rad)    col 1: vel0 (rad/s)
 *   col 2: pos1 (rad)    col 3: vel1 (rad/s)
 *   ...
 *   col 2k: pos_k        col 2k+1: vel_k
 */

#include "training.hpp"
#include "sg_filter.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <promp/promp.hpp>
#include <promp/trajectory.hpp>
#include <promp/io/serializer.hpp>

namespace fs = std::filesystem;

namespace promp_rt {

// ============================================================
// Internal helpers
// ============================================================

namespace {

/**
 * @brief Piecewise-linear 1-D interpolation.
 *
 * @param x_in   Source abscissae, sorted ascending.
 * @param y_in   Source ordinates (same length as x_in).
 * @param x_out  Query abscissae (any order).
 *
 * @return  Interpolated ordinates; values outside [x_in.front(), x_in.back()]
 *          are clamped to the boundary values.
 */
static std::vector<double> interp1(const std::vector<double>& x_in,
                                   const std::vector<double>& y_in,
                                   const std::vector<double>& x_out)
{
    std::vector<double> result;
    result.reserve(x_out.size());
    for (double xq : x_out) {
        if (xq <= x_in.front()) { result.push_back(y_in.front()); continue; }
        if (xq >= x_in.back())  { result.push_back(y_in.back());  continue; }
        auto it  = std::upper_bound(x_in.begin(), x_in.end(), xq);
        size_t hi = static_cast<size_t>(it - x_in.begin());
        size_t lo = hi - 1;
        double t  = (xq - x_in[lo]) / (x_in[hi] - x_in[lo]);
        result.push_back(y_in[lo] + t * (y_in[hi] - y_in[lo]));
    }
    return result;
}

/**
 * @brief Parse a single CSV demonstration file.
 *
 * Expected column layout (n_dof+1 columns):
 *   time_s,  pos0_rad,  pos1_rad,  ...,  pos(n_dof-1)_rad
 *
 * An optional text header is auto-detected by attempting to parse the
 * first token as a double; if that fails the row is treated as a header.
 *
 * @param[in]  path          File path.
 * @param[in]  n_dof         Expected number of position columns.
 * @param[out] times         Parsed time stamps (s), in original order.
 * @param[out] positions     positions[d] = raw position vector for joint d
 *                           (d = 0 … n_dof-1).
 * @param[in]  max_duration  Rows with time > max_duration are discarded.
 *                           ≤ 0 → keep all rows.
 *
 * @return true  when ≥ 5 valid rows are parsed.
 *         false on open failure or insufficient valid rows.
 *
 * ── Data-frequency note ────────────────────────────────────────────────────
 * There is no assumption on the time stamp spacing; the file may contain
 * non-uniform or varying-rate data.
 */
static bool parse_csv(const fs::path& path,
                      int n_dof,
                      std::vector<double>& times,
                      std::vector<std::vector<double>>& positions,
                      double max_duration)
{
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "[train] Cannot open: " << path << "\n";
        return false;
    }

    times.clear();
    positions.assign(static_cast<size_t>(n_dof), {});

    std::string line;
    bool header_skipped = false;

    while (std::getline(f, line)) {
        if (line.empty()) continue;

        std::istringstream ss(line);
        std::string tok;
        std::vector<std::string> tokens;
        while (std::getline(ss, tok, ',')) tokens.push_back(tok);
        if (static_cast<int>(tokens.size()) < 1 + n_dof) continue;

        double t_val{};
        try { t_val = std::stod(tokens[0]); }
        catch (...) {
            if (!header_skipped) { header_skipped = true; continue; }
            continue; // skip malformed row
        }
        header_skipped = true;

        if (max_duration > 0.0 && t_val > max_duration) break;

        bool row_ok = true;
        std::vector<double> p(static_cast<size_t>(n_dof));
        for (int d = 0; d < n_dof; ++d) {
            try { p[static_cast<size_t>(d)] = std::stod(tokens[static_cast<size_t>(d + 1)]); }
            catch (...) { row_ok = false; break; }
        }
        if (!row_ok) continue;

        times.push_back(t_val);
        for (int d = 0; d < n_dof; ++d)
            positions[static_cast<size_t>(d)].push_back(p[static_cast<size_t>(d)]);
    }

    if (static_cast<int>(times.size()) < 5) {
        std::cerr << "[train] Too few valid rows (<5) in: " << path << "\n";
        return false;
    }
    return true;
}

/**
 * @brief Build a 2*n_dof–dimensional promp::Trajectory from raw demo data.
 *
 * Steps:
 *   1. Normalise time to [0, 1] (raw data, before resampling).
 *   2. Compute velocity for each joint using sg_filter_deriv on the normalised
 *      time axis, then divide by the physical duration to convert to rad/s.
 *   3. Build a uniform query grid of n_steps points over [0, 1].
 *   4. Interpolate all 2*n_dof channels to the query grid.
 *   5. Pack into an Eigen matrix and wrap in a promp::Trajectory.
 *
 * ── Why SG before resampling ─────────────────────────────────────────────────
 * Running the SG filter on the ORIGINAL samples preserves the velocity
 * information at the original data rate.  Resampling first and then
 * differentiating would lose high-frequency content and produce artefacts
 * at the boundary of the resampled grid.
 *
 * @param times      Raw time stamps (s) of the demonstration.
 * @param positions  positions[d] = raw position vector for joint d.
 * @param n_dof      Number of joints.
 * @param n_steps    Number of output steps (uniform phase grid).
 * @param sg_window  SG filter window length (odd, ≥ sg_poly+2).
 *                   Should be chosen relative to the data rate:
 *                     sg_window ≈ data_rate_Hz × smoothing_time_s (odd).
 * @param sg_poly    SG polynomial order (< sg_window − 1).
 *
 * @return promp::Trajectory with dims = 2*n_dof and timesteps = n_steps.
 *         Column layout: [pos0, vel0, pos1, vel1, …, pos(n-1), vel(n-1)].
 */
static promp::Trajectory build_trajectory(const std::vector<double>& times,
                                          const std::vector<std::vector<double>>& positions,
                                          int n_dof,
                                          int n_steps,
                                          int sg_window,
                                          int sg_poly)
{
    // 1. Normalise time to [0, 1]
    const double t0   = times.front();
    double span       = times.back() - t0;
    if (span <= 0.0) span = 1.0;

    std::vector<double> t_norm(times.size());
    for (size_t i = 0; i < times.size(); ++i)
        t_norm[i] = (times[i] - t0) / span;

    // 2. Compute velocities in NORMALISED-time units, then divide by span
    //    to convert to physical rad/s.
    //
    //    Note: sg_filter_deriv returns dy/d(t_norm), so physical velocity =
    //          dy/dt = (dy/d(t_norm)) / span.
    std::vector<std::vector<double>> velocities(static_cast<size_t>(n_dof));
    for (int d = 0; d < n_dof; ++d) {
        const auto& pos_d = positions[static_cast<size_t>(d)];
        std::vector<double> v_norm = sg_filter_deriv(t_norm, pos_d, sg_window, sg_poly);
        velocities[static_cast<size_t>(d)].resize(v_norm.size());
        for (size_t i = 0; i < v_norm.size(); ++i)
            velocities[static_cast<size_t>(d)][i] = v_norm[i] / span;
    }

    // 3. Uniform query grid over normalised phase [0, 1]
    std::vector<double> t_common(static_cast<size_t>(n_steps));
    for (int i = 0; i < n_steps; ++i)
        t_common[static_cast<size_t>(i)] = static_cast<double>(i) / (n_steps - 1);

    // 4 & 5. Interpolate and pack into Eigen matrix (n_steps × 2*n_dof)
    //        Column layout: [pos0, vel0, pos1, vel1, …]
    Eigen::MatrixXd data(n_steps, 2 * n_dof);
    for (int d = 0; d < n_dof; ++d) {
        std::vector<double> pos_r = interp1(t_norm, positions[static_cast<size_t>(d)], t_common);
        std::vector<double> vel_r = interp1(t_norm, velocities[static_cast<size_t>(d)], t_common);
        for (int i = 0; i < n_steps; ++i) {
            data(i, d * 2)     = pos_r[static_cast<size_t>(i)];
            data(i, d * 2 + 1) = vel_r[static_cast<size_t>(i)];
        }
    }
    return promp::Trajectory(data);
}

/**
 * @brief Write the training metadata sidecar file (key=value format).
 *
 * All parameters needed by OnlinePredictor are stored here so that the
 * predictor is self-contained given only the model directory.
 *
 * @param path          Output path (typically model_dir/metadata.txt).
 * @param cfg           TrainingConfig used for this run.
 * @param mean_duration Mean demo duration across all loaded demonstrations (s).
 * @param n_demos       Number of demonstrations successfully loaded.
 *
 * @throws std::runtime_error if the file cannot be written.
 */
static void write_metadata(const fs::path& path,
                            const TrainingConfig& cfg,
                            double mean_duration,
                            int n_demos)
{
    std::ofstream f(path);
    if (!f.is_open())
        throw std::runtime_error("Cannot write metadata: " + path.string());

    f << "n_dof="              << cfg.n_dof               << "\n"
      << "n_steps="            << cfg.n_resample_steps     << "\n"
      << "num_basis="          << cfg.num_basis_functions  << "\n"
      << "std_bf="             << cfg.std_bf               << "\n"
      << "mean_duration_s="    << mean_duration            << "\n"
      << "sg_window="          << cfg.sg_window            << "\n"
      << "sg_poly_order="      << cfg.sg_poly_order        << "\n"
      << "obs_pos_noise_var="  << cfg.obs_pos_noise_var    << "\n"
      << "obs_vel_noise_var="  << cfg.obs_vel_noise_var    << "\n"
      << "goal_pos_noise_var=" << cfg.goal_pos_noise_var   << "\n"
      << "goal_vel_noise_var=" << cfg.goal_vel_noise_var   << "\n"
      << "n_demos="            << n_demos                  << "\n";
}

} // anonymous namespace

// ============================================================
// Public API
// ============================================================

/**
 * @brief Train a ProMP from a folder of CSV demonstrations and save the model.
 *
 * ── Input ─────────────────────────────────────────────────────────────────────
 * @param demo_folder
 *   Directory containing *.csv files.  Each file is one demonstration.
 *   Column format: time_s, pos0_rad, …, pos(n_dof-1)_rad.
 *   Files are processed in lexicographic order.
 *   Demos with < 5 valid rows are silently skipped.
 *   Different files may have different sampling rates (handled transparently).
 *
 * @param model_dir
 *   Output directory (created if absent).  Writes:
 *     model.promp   – serialised ProMP model (promp::io::save_promp format).
 *     metadata.txt  – all training parameters for self-contained inference.
 *
 * @param config
 *   Training hyper-parameters.  Key fields for a first run:
 *     n_dof                – must match the number of position columns in CSVs.
 *     num_basis_functions  – start at 10; increase for complex trajectories.
 *     n_resample_steps     – 100 is sufficient for most applications.
 *     sg_window            – set to ~data_rate_Hz × 0.09 s (≈ 9 for 100 Hz).
 *     max_demo_duration    – set to expected motion time to trim tail noise.
 *
 * @return true on success; false if no valid demos are found or an I/O
 *         error occurs (error message printed to stderr).
 *
 * ── Data-frequency compatibility ─────────────────────────────────────────────
 * Demos at different sampling rates can be mixed freely.  The SG filter runs
 * on original timestamps so velocity quality reflects the original rate.
 * The common phase grid (n_resample_steps) is always uniform regardless of
 * input rates.
 */
bool train_from_folder(const std::string& demo_folder,
                       const std::string& model_dir,
                       const TrainingConfig& config)
{
    const int n_dof = config.n_dof;

    fs::path demo_path(demo_folder);
    if (!fs::exists(demo_path) || !fs::is_directory(demo_path)) {
        std::cerr << "[train] Demo folder not found: " << demo_folder << "\n";
        return false;
    }

    // Collect CSV files in lexicographic order for reproducibility.
    std::vector<fs::path> csv_files;
    for (const auto& entry : fs::directory_iterator(demo_path))
        if (entry.path().extension() == ".csv")
            csv_files.push_back(entry.path());
    std::sort(csv_files.begin(), csv_files.end());

    if (csv_files.empty()) {
        std::cerr << "[train] No *.csv files found in: " << demo_folder << "\n";
        return false;
    }

    std::vector<promp::Trajectory> trajectories;
    std::vector<double> durations;
    trajectories.reserve(csv_files.size());

    for (const auto& csv : csv_files) {
        std::vector<double> times;
        std::vector<std::vector<double>> positions;

        if (!parse_csv(csv, n_dof, times, positions, config.max_demo_duration))
            continue;

        const double duration = times.back() - times.front();
        if (duration <= 0.0) continue;

        trajectories.push_back(build_trajectory(
            times, positions,
            n_dof,
            config.n_resample_steps,
            config.sg_window,
            config.sg_poly_order));
        durations.push_back(duration);

        std::cout << "[train] Loaded: " << csv.filename()
                  << "  (" << times.size() << " pts, "
                  << duration << " s)\n";
    }

    if (trajectories.empty()) {
        std::cerr << "[train] No valid trajectories loaded.\n";
        return false;
    }

    // Mean demo duration – used by the predictor to map real time → phase step.
    double mean_duration = 0.0;
    for (double d : durations) mean_duration += d;
    mean_duration /= static_cast<double>(durations.size());

    std::cout << "[train] Training ProMP on " << trajectories.size()
              << " demos  (n_dof=" << n_dof
              << ", dims=" << 2 * n_dof
              << ", mean_duration=" << mean_duration << " s)...\n";

    // Train: ridge regression per demo → weight vectors → Gaussian fit (μ_w, Σ_w).
    promp::ProMP model(trajectories, config.num_basis_functions, config.std_bf);

    // Persist.
    fs::create_directories(model_dir);
    const fs::path model_file = fs::path(model_dir) / "model.promp";
    promp::io::save_promp(model_file.string(), model);
    write_metadata(fs::path(model_dir) / "metadata.txt",
                   config, mean_duration,
                   static_cast<int>(trajectories.size()));

    std::cout << "[train] Saved → " << model_file << "\n";
    std::cout << "[train] Done.  basis=" << config.num_basis_functions
              << "  steps=" << config.n_resample_steps
              << "  dims="  << model.get_dims() << "\n";
    return true;
}

} // namespace promp_rt
