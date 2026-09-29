#include "ukf_predictor.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

// =============================================================================
// CausalVelocity
// =============================================================================

CausalVelocity::CausalVelocity(std::size_t window) : window_(window) {}

void CausalVelocity::reset() {
    t_.clear();
    y_.clear();
}

double CausalVelocity::push(double t, double y) {
    t_.push_back(t);
    y_.push_back(y);
    if (t_.size() > window_) {
        t_.pop_front();
        y_.pop_front();
    }

    const std::size_t n = t_.size();
    if (n < 2) return 0.0;

    // Ordinary least-squares slope of y vs. t over the current window:
    // slope = sum((t-meanT)*(y-meanY)) / sum((t-meanT)^2).
    double meanT = 0.0, meanY = 0.0;
    for (std::size_t i = 0; i < n; ++i) { meanT += t_[i]; meanY += y_[i]; }
    meanT /= static_cast<double>(n);
    meanY /= static_cast<double>(n);

    double sxy = 0.0, sxx = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double dt = t_[i] - meanT;
        sxy += dt * (y_[i] - meanY);
        sxx += dt * dt;
    }
    if (sxx < 1e-12) return 0.0;
    return sxy / sxx;
}

// =============================================================================
// estimateGoalMinJerk
// =============================================================================

// Fit y = c2*t^2 + c1*t + c0 via least squares, solving the normal
// equations (A^T A) c = A^T y with Eigen's 3x3 solver. Internal helper used
// only to estimate local acceleration for the min-jerk extrapolation below.
static void quadFit(const std::vector<double>& t, const std::vector<double>& y,
                     double& c2, double& c1, double& c0) {
    const std::size_t n = t.size();
    Eigen::Matrix3d ata = Eigen::Matrix3d::Zero();
    Eigen::Vector3d aty = Eigen::Vector3d::Zero();

    for (std::size_t i = 0; i < n; ++i) {
        Eigen::Vector3d row(t[i] * t[i], t[i], 1.0);
        ata += row * row.transpose();
        aty += row * y[i];
    }

    const Eigen::Vector3d c = ata.ldlt().solve(aty);
    c2 = c(0);
    c1 = c(1);
    c0 = c(2);
}

double estimateGoalMinJerk(const std::vector<double>& tObs,
                            const std::vector<double>& yObs,
                            double vObs,
                            double predictionWindow) {
    if (tObs.size() < 3) {
        return yObs.empty() ? 0.0 : yObs.back();
    }

    const double x0 = yObs.back();
    const double v0 = vObs;
    double a0 = 0.0;

    // Estimate local acceleration from a short quadratic fit to the most
    // recent samples (needs at least 5 points to be well-conditioned).
    if (tObs.size() >= 5) {
        const std::size_t nFit = std::min<std::size_t>(9, tObs.size());
        std::vector<double> tRel(nFit), yRecent(nFit);
        const std::size_t start = tObs.size() - nFit;
        const double t0 = tObs[start];
        for (std::size_t i = 0; i < nFit; ++i) {
            tRel[i] = tObs[start + i] - t0;
            yRecent[i] = yObs[start + i];
        }
        double c2, c1, c0;
        quadFit(tRel, yRecent, c2, c1, c0);
        a0 = 2.0 * c2;
    }

    // Closed-form minimum-jerk extrapolation to a smooth stop over a fixed
    // window T: goal ~= x0 + 0.6*v0*T + 0.1*a0*T^2 (see kf_simulation.py).
    const double T = std::max(predictionWindow, 0.01);
    const double T2 = T * T;
    return x0 + 0.6 * v0 * T + 0.1 * a0 * T2;
}

// =============================================================================
// Ukf1D
// =============================================================================

Ukf1D::Ukf1D(double initialPosition, const UkfParams& params)
    : params_(params),
      x_(initialPosition, 0.0, params.initialGain, initialPosition),
      P_(params.initialCovariance.asDiagonal()),
      Q_(params.processNoise.asDiagonal()),
      R_(params.measurementNoise.asDiagonal()),
      velEstimator_(params.velocityWindow) {
    // Scaled unscented transform weights (Van der Merwe), n = kNx = 4:
    //   lambda = alpha^2 * (n + kappa) - n
    //   c      = n + lambda                  (sigma-point spread scale)
    //   wm[0]  = lambda / c,           wc[0] = wm[0] + (1 - alpha^2 + beta)
    //   wm[i]  = wc[i] = 1 / (2c)      for i = 1..2n
    const double n = static_cast<double>(kNx);
    const double lambda = params_.alpha * params_.alpha * (n + params_.kappa) - n;
    c_ = n + lambda;

    wm_(0) = lambda / c_;
    wc_(0) = wm_(0) + (1.0 - params_.alpha * params_.alpha + params_.beta);
    for (int i = 1; i < kNumSigma; ++i) {
        wm_(i) = 1.0 / (2.0 * c_);
        wc_(i) = wm_(i);
    }
}

Eigen::Vector4d Ukf1D::propagateSimultaneous(const Eigen::Vector4d& state, double dt) {
    const double pos = state(0);
    const double vel = state(1);
    const double L = state(2);
    const double goal = state(3);
    return {pos + vel * dt, vel + L * (goal - pos) * dt, L, goal};
}

void Ukf1D::update(double yObs, double t, double dt, double vObsIfMeasured) {
    // Configuration knob #1: pick the velocity measurement source.
    const double vEst = (params_.velocitySource == VelocitySource::kMeasured)
                             ? vObsIfMeasured
                             : velEstimator_.push(t, yObs);

    // ---- 1) Causal measurement construction -------------------------------
    // Position, the chosen velocity, and an online minimum-jerk estimate of
    // where the movement is heading (the "goal" pseudo-measurement).
    tHistory_.push_back(t);
    yHistory_.push_back(yObs);
    const double goalEst = estimateGoalMinJerk(tHistory_, yHistory_, vEst,
                                                params_.goalPredictionWindow);

    // ---- 2) Sigma points ----------------------------------------------------
    // X = [x, x + col_i(S), x - col_i(S)], where S is the lower-triangular
    // Cholesky factor of the scaled covariance c*P (2*kNx+1 sigma points
    // total, matching the standard unscented transform construction).
    const Eigen::Matrix4d S = Eigen::LLT<Eigen::Matrix4d>(P_ * c_).matrixL();
    Eigen::Matrix<double, kNx, kNumSigma> X;
    X.col(0) = x_;
    for (int i = 0; i < kNx; ++i) {
        X.col(1 + i) = x_ + S.col(i);
        X.col(1 + kNx + i) = x_ - S.col(i);
    }

    // ---- 3) Predict step (process model) -------------------------------
    // Propagate every sigma point through the nonlinear process model, then
    // recombine into a predicted mean/covariance via the UKF weights.
    Eigen::Matrix<double, kNx, kNumSigma> Xp;
    for (int i = 0; i < kNumSigma; ++i) Xp.col(i) = propagateSimultaneous(X.col(i), dt);

    const Eigen::Vector4d xPred = Xp * wm_;

    Eigen::Matrix4d pPred = Q_;
    for (int i = 0; i < kNumSigma; ++i) {
        const Eigen::Vector4d dx = Xp.col(i) - xPred;
        pPred += wc_(i) * (dx * dx.transpose());
    }

    // ---- 4) Measurement update -------------------------------------------
    // Observe z = [pos, vel, goal]: the min-jerk goal estimate is fused in
    // exactly like any other noisy sensor reading, via the same
    // cross-covariance/Kalman-gain machinery as the direct pos/vel
    // observations.
    constexpr int kNz = 3;
    Eigen::Matrix<double, kNz, kNumSigma> Z;
    Z.row(0) = Xp.row(0);
    Z.row(1) = Xp.row(1);
    Z.row(2) = Xp.row(3);

    const Eigen::Vector3d zPred = Z * wm_;

    Eigen::Matrix3d rAug = Eigen::Matrix3d::Zero();
    rAug(0, 0) = R_(0, 0);
    rAug(1, 1) = R_(1, 1);
    rAug(2, 2) = params_.goalMeasurementNoise;

    Eigen::Matrix3d pzz = rAug;
    Eigen::Matrix<double, kNx, kNz> pxz = Eigen::Matrix<double, kNx, kNz>::Zero();
    for (int i = 0; i < kNumSigma; ++i) {
        const Eigen::Vector4d dx = Xp.col(i) - xPred;
        const Eigen::Vector3d dz = Z.col(i) - zPred;
        pzz += wc_(i) * (dz * dz.transpose());
        pxz += wc_(i) * (dx * dz.transpose());
    }

    const Eigen::Matrix<double, kNx, kNz> K = pxz * pzz.inverse();

    const Eigen::Vector3d zMeas(yObs, vEst, goalEst);
    const Eigen::Vector3d innovation = zMeas - zPred;

    x_ = xPred + K * innovation;
    P_ = pPred - K * pzz * K.transpose();

    x_(2) = std::max(x_(2), 0.01); // human gain L must stay positive
}

std::tuple <std::vector<double>, std::vector<double>> Ukf1D::forecast(int horizonSteps, double dt) const {
    std::vector<double> pred(static_cast<std::size_t>(horizonSteps));
    std::vector<double> pred_vel(static_cast<std::size_t>(horizonSteps));
    Eigen::Vector4d xf = x_;
    for (int j = 0; j < horizonSteps; ++j) {
        // Sequential (semi-implicit) Euler rollout, matching the forecast
        // loop in kf_simulation.py: velocity uses the *already updated*
        // position for this step.
        xf(0) = xf(0) + xf(1) * dt;
        xf(1) = xf(1) + xf(2) * (xf(3) - xf(0)) * dt;
        pred[static_cast<std::size_t>(j)] = xf(0);
        pred_vel[static_cast<std::size_t>(j)] = xf(1);
    }
    return make_tuple(pred,pred_vel);
}

// =============================================================================
// MultiDofUkfPredictor
// =============================================================================

MultiDofUkfPredictor::MultiDofUkfPredictor(std::size_t numDoF,
                                            const std::vector<double>& initialPositions,
                                            const UkfParams& params)
    : params_(params) {
    if (initialPositions.size() != numDoF) {
        throw std::invalid_argument(
            "MultiDofUkfPredictor: initialPositions.size() must equal numDoF");
    }
    filters_.reserve(numDoF);
    for (double p0 : initialPositions) filters_.emplace_back(p0, params_);
}

std::vector<std::vector<double>> MultiDofUkfPredictor::updateAndPredict(
    const std::vector<double>& position, double t, double dt, double horizonSeconds,
    const std::vector<double>& velocity) {
    if (position.size() != filters_.size()) {
        throw std::invalid_argument("updateAndPredict: position size must match numDoF()");
    }
    if (params_.velocitySource == VelocitySource::kMeasured && velocity.size() != filters_.size()) {
        throw std::invalid_argument(
            "updateAndPredict: velocitySource is kMeasured, so velocity must be supplied "
            "with one value per DoF");
    }

    //std::cout << position.size() << " Before horizonSteps ";

    const int horizonSteps = std::max(1, static_cast<int>(std::round(horizonSeconds / dt)));

    std::vector<std::vector<double>> predicted(filters_.size()*2);
    for (std::size_t i = 0; i < filters_.size(); ++i) {
        const double vObs = (params_.velocitySource == VelocitySource::kMeasured) ? velocity[i] : 0.0;
        filters_[i].update(position[i], t, dt, vObs);
        tie(predicted[i], predicted[i+filters_.size()]) = filters_[i].forecast(horizonSteps, dt);
    }
    return predicted;
}

std::vector<std::vector<double>> MultiDofUkfPredictor::predictFromHistory(
    const std::vector<std::vector<double>>& obsHistory,
    const std::vector<double>& timestamps,
    double horizonSeconds) {
    if (obsHistory.empty()) throw std::invalid_argument("predictFromHistory: empty history");
    if (obsHistory.size() != timestamps.size()) {
        throw std::invalid_argument("predictFromHistory: history/timestamps size mismatch");
    }

    // Reset every filter to the first observation in the supplied history.
    const std::size_t numDoF = obsHistory.front().size();
    filters_.clear();
    filters_.reserve(numDoF);
    for (double p0 : obsHistory.front()) filters_.emplace_back(p0, params_);

    std::vector<std::vector<double>> lastForecast;
    for (std::size_t k = 1; k < obsHistory.size(); ++k) {
        const double dt = timestamps[k] - timestamps[k - 1];
        lastForecast = updateAndPredict(obsHistory[k], timestamps[k], dt, horizonSeconds);
    }
    return lastForecast;
}

std::vector<double> MultiDofUkfPredictor::currentPosition() const {
    std::vector<double> pos(filters_.size());
    for (std::size_t i = 0; i < filters_.size(); ++i) pos[i] = filters_[i].position();
    return pos;
}
