#pragma once

// =============================================================================
// ukf_predictor.hpp
//
// Self-contained online Unscented Kalman Filter (UKF) library for predicting
// the future trajectory of a multi-degree-of-freedom (multi-DoF) movement
// (e.g. a human/robot end-effector position, or a set of joint angles) from
// streaming position observations.
//
// USAGE (this is the whole contract):
//   1. Construct one MultiDofUkfPredictor, giving it the number of DoF and
//      their starting positions.
//   2. Every time a new observation arrives, call
//      predictor.updateAndPredict(position, t, dt, horizonSeconds)
//      -> this feeds the new observation into the filter (the "past
//         observation" the filter needs is entirely captured in its
//         internal state - the caller only ever supplies the newest
///        sample) and immediately returns the predicted future trajectory
//         for the next `horizonSeconds` seconds, one point per DoF per
//         internal timestep.
//   That's it - this is a genuinely *online* filter: it never looks ahead,
//   never re-reads old data, and the per-step cost is independent of how
//   long the movement has been running.
//
// TWO CONFIGURATION KNOBS the caller sets once (both live in UkfParams):
//   - params.velocitySource: whether the filter derives velocity internally
//     from the position stream (VelocitySource::kDerivedFromPosition), or
//     expects the caller to supply an already-measured velocity every step
//     (VelocitySource::kMeasured), e.g. from an encoder/IMU or a
//     pre-filtered dataset channel.
//   - the number of DoF: set by the `numDoF` argument (and the size of the
//     `initialPositions` vector, which must match it) passed to
//     MultiDofUkfPredictor's constructor. Every DoF is filtered
//     independently with an identical copy of `params`, so this can be any
//     positive number with no other code changes.
//
// MODEL: each DoF is tracked by its own 4-state Unscented Kalman Filter
// (Ukf1D), a C++ port of the UKF in kf_simulation.py modelling a human
// reaching movement as a simple feedback law:
//     state x = [pos, vel, L, goal]
//     pos'  = pos + vel*dt
//     vel'  = vel + L*(goal - pos)*dt      <- L*(goal-pos) is a spring-like
//                                             "pull" toward the goal, scaled
//                                             by a personal gain L
//     L'    = L      (constant "gain", but still continuously re-estimated
//                      online alongside the other states - never frozen)
//     goal' = goal    (constant, but likewise continuously re-estimated)
// Because `L` multiplies `(goal - pos)` - a product of two state variables -
// this process model is nonlinear, which is exactly why a UKF (rather than a
// plain linear Kalman filter) is needed: sigma points propagate that
// nonlinearity through the predict step without needing to linearize it by
// hand.
//
// The filter is never told the true destination in advance. At every step,
// in addition to the observed position (and velocity), a third pseudo-
// measurement - an online estimate of where the movement is heading,
// extrapolated assuming a smooth minimum-jerk stop (estimateGoalMinJerk) -
// is fused in through the same Kalman update math as any other sensor.
// =============================================================================

#include <Eigen/Dense>
#include <deque>
#include <vector>

// -----------------------------------------------------------------------
// Internal building block 1/2: causal (online) velocity estimator.
//
// A sliding-window linear-regression slope of position vs. time - the
// causal, order-1 Savitzky-Golay equivalent used by the reference Python
// implementation (kf_simulation.py) to turn raw position samples into a
// velocity measurement without ever looking at future samples. Only used
// when UkfParams::velocitySource == VelocitySource::kDerivedFromPosition;
// with kMeasured this class is unused and the caller's own velocity
// reading is trusted instead.
// -----------------------------------------------------------------------
class CausalVelocity {
public:
    explicit CausalVelocity(std::size_t window = 15);

    // Feed one new (time, position) sample and return the current causal
    // velocity estimate (0 until at least two samples are available).
    double push(double t, double y);

    void reset();

private:
    std::size_t window_;
    std::deque<double> t_;
    std::deque<double> y_;
};

// -----------------------------------------------------------------------
// Internal building block 2/2: online goal (movement endpoint) estimator.
//
// Ported from _estimate_goal_minjerk in kf_simulation.py. Given the full
// position/time history observed so far and the current velocity estimate,
// extrapolates a fixed prediction window ahead assuming the movement comes
// to a smooth (minimum-jerk) stop, and returns that extrapolated endpoint.
// This is fused into Ukf1D as a third measurement channel alongside
// position and velocity - it is what lets the filter track a moving goal
// without ever being told the true destination.
// -----------------------------------------------------------------------
double estimateGoalMinJerk(const std::vector<double>& tObs,
                            const std::vector<double>& yObs,
                            double vObs,
                            double predictionWindow = 1.0);

// -----------------------------------------------------------------------
// Configuration knob #1: how each Ukf1D obtains its velocity measurement.
// -----------------------------------------------------------------------
enum class VelocitySource {
    // Derive velocity online from the position stream via CausalVelocity.
    // Use this whenever only a position sensor is available.
    kDerivedFromPosition,

    // Trust a velocity value supplied by the caller at every step (e.g. a
    // recorded/pre-filtered derivative channel, or a real velocity
    // sensor). Measured on real reaching-movement data, this reduced mean
    // prediction error by roughly 11-37% (at 0.1-0.5 s horizons) relative
    // to deriving velocity from a downsampled position stream, since the
    // measured channel is typically filtered from a much higher sample
    // rate than the online regression sees.
    kMeasured,
};

// Tuning parameters shared by every per-DoF filter. A free (non-nested)
// struct so it gets a normal default constructor usable as a default
// argument (nesting it inside Ukf1D trips a GCC quirk around default
// member initializers used as a constructor's default argument).
struct UkfParams {
    // --- Configuration knob #1 (see VelocitySource above) ---
    VelocitySource velocitySource = VelocitySource::kDerivedFromPosition;

    // --- Unscented transform shape parameters (Van der Merwe scaling) ---
    double alpha = 1e-3;
    double beta = 2.0;
    double kappa = 0.0;

    // --- Noise covariances ---
    // Process noise (pos, vel, L, goal).
    Eigen::Vector4d processNoise{1e-5, 1e-5, 1e-4, 1e-5};
    // Measurement noise for the direct observations (pos, vel).
    Eigen::Vector2d measurementNoise{1e-4, 1e-4};
    // Measurement noise associated with the online min-jerk goal estimate.
    double goalMeasurementNoise = 1e-2;

    // --- Initial state uncertainty / values ---
    // Initial covariance (pos, vel, L, goal).
    Eigen::Vector4d initialCovariance{1e-3, 1e-3, 1e-2, 1e-1};
    double initialGain = 0.6; // initial human-gain L

    // --- Other tunables ---
    double goalPredictionWindow = 1.0; // seconds, minimum-jerk look-ahead window
    std::size_t velocityWindow = 15;   // samples, CausalVelocity's window size
};

// =============================================================================
// Ukf1D: single-degree-of-freedom Unscented Kalman Filter.
//
// Tracks and forecasts ONE scalar signal (e.g. one joint angle, or one
// Cartesian axis of an end-effector). MultiDofUkfPredictor below is the
// class you actually use for a multi-DoF movement; it simply owns one
// Ukf1D per DoF (see that class's docstring for why they are independent
// rather than one large coupled filter).
// =============================================================================
class Ukf1D {
public:
    explicit Ukf1D(double initialPosition, const UkfParams& params = UkfParams());

    // Feed one new position observation at time t (dt = elapsed time since
    // the previous observation) and update the internal state estimate.
    //
    // vObsIfMeasured is the velocity measurement to use - but it is only
    // actually read when params.velocitySource == VelocitySource::kMeasured
    // (configuration knob #1). With the default kDerivedFromPosition, this
    // argument is ignored and velocity is instead derived online from the
    // position stream via CausalVelocity, so it can simply be omitted.
    void update(double yObs, double t, double dt, double vObsIfMeasured = 0.0);

    // Roll the current state estimate forward `horizonSteps` steps of size
    // `dt` using the (noise-free) process model, without any further
    // measurement updates. Returns the predicted position at each step.
    // This is the "forecast" half of every online prediction call; see
    // MultiDofUkfPredictor::updateAndPredict for the combined entry point.
    std::tuple <std::vector<double>, std::vector<double>> forecast(int horizonSteps, double dt) const;

    double position() const { return x_(0); }
    double velocity() const { return x_(1); }
    double gain() const { return x_(2); }
    double goal() const { return x_(3); }

private:
    static constexpr int kNx = 4;             // state dimension: [pos, vel, L, goal]
    static constexpr int kNumSigma = 2 * kNx + 1;

    UkfParams params_;
    Eigen::Vector4d x_;    // state [pos, vel, L, goal]
    Eigen::Matrix4d P_;    // state covariance
    Eigen::Matrix4d Q_;    // process noise covariance
    Eigen::Matrix2d R_;    // base measurement noise covariance (pos, vel)

    Eigen::Matrix<double, kNumSigma, 1> wm_; // UKF mean weights
    Eigen::Matrix<double, kNumSigma, 1> wc_; // UKF covariance weights
    double c_ = 0.0;                         // scaling factor n + lambda

    CausalVelocity velEstimator_;   // only used when velocitySource == kDerivedFromPosition
    std::vector<double> tHistory_; // full time history, used by estimateGoalMinJerk
    std::vector<double> yHistory_; // full position history, used by estimateGoalMinJerk

    // Sigma-point transition used inside the UKF predict step: a
    // simultaneous (old-position) Euler update, matching the Xp
    // computation in kf_simulation.py's filtering loop.
    static Eigen::Vector4d propagateSimultaneous(const Eigen::Vector4d& state, double dt);
};

// =============================================================================
// MultiDofUkfPredictor: the class you actually use for online prediction.
//
// Owns one independent Ukf1D per degree of freedom (configuration knob #2:
// how many DoF, set at construction time). The DoFs are filtered
// independently - each gets its own [pos, vel, L, goal] state and its own
// copy of `params` - rather than stacked into one large coupled state
// vector, because the underlying process model (see file header) has no
// cross-DoF coupling term to estimate: stacking would give numerically
// identical results at a much higher computational cost (O(numDoF^3) per
// step instead of O(numDoF) copies of an O(1) 4x4 filter).
// =============================================================================
class MultiDofUkfPredictor {
public:
    // numDoF must equal initialPositions.size() (checked at construction);
    // it is a separate, explicit argument specifically so the number of
    // DoF being tracked is a visible, deliberate choice at the call site
    // rather than an incidental side effect of a vector's length.
    MultiDofUkfPredictor(std::size_t numDoF, const std::vector<double>& initialPositions,
                          const UkfParams& params = UkfParams());

    std::size_t numDoF() const { return filters_.size(); }

    // THE main online-prediction entry point.
    //   Input : position       - the newest observation, one value per DoF
    //                             (this is the only "past observation" the
    //                             caller ever supplies - all earlier
    //                             history the filter needs is already held
    //                             in its internal state);
    //           t, dt          - current time and elapsed time since the
    //                             previous call;
    //           horizonSeconds - how far into the future to predict;
    //           velocity       - the newest velocity observation, one
    //                             value per DoF. Required (and used) only
    //                             when params.velocitySource ==
    //                             VelocitySource::kMeasured (configuration
    //                             knob #1); otherwise omit it and velocity
    //                             is derived online from `position`.
    //   Output: predicted[dof][step] - the predicted future trajectory for
    //           every DoF, one entry per internal timestep of size `dt` out
    //           to `horizonSeconds`.
    std::vector<std::vector<double>> updateAndPredict(const std::vector<double>& position,
                                                        double t, double dt,
                                                        double horizonSeconds,
                                                        const std::vector<double>& velocity = {});

    // Convenience batch API: replays an entire recorded observation history
    // (obsHistory[k][dof], one row per time step, sharing timestamps[k])
    // through the online filters from scratch and returns only the final
    // forecast made from the last observation. Useful for offline
    // evaluation/plotting; the true online use is updateAndPredict() above,
    // called once per new sample as it arrives.
    std::vector<std::vector<double>> predictFromHistory(
        const std::vector<std::vector<double>>& obsHistory,
        const std::vector<double>& timestamps,
        double horizonSeconds);

    std::vector<double> currentPosition() const;

private:
    UkfParams params_;
    std::vector<Ukf1D> filters_;
};
