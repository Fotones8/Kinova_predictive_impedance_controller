# promp_realtime

Real-time **n-DoF** trajectory prediction using the **ProMP_1** library
(INRIA / hucebot/promp).

The method trains a `2×n_dof`-dimensional Probabilistic Movement Primitive
(joint positions + joint velocities) from CSV demonstrations.
Velocities are computed automatically from measured positions using the
**Savitzky-Golay filter**.  At runtime, incoming sensor measurements are
used to update the trajectory prediction via Bayesian conditioning.

The number of degrees of freedom (`n_dof`) is set at training time and stored
in the model — no recompile is needed to change it.

---

## Folder structure

```
promp_realtime/
├── CMakeLists.txt          # Build system
├── README.md               # This file
├── include/
│   ├── training.hpp        # Training API  (TrainingConfig, train_from_folder)
│   ├── prediction.hpp      # Prediction API (OnlinePredictor, PredictionResult)
│   └── sg_filter.hpp       # Header-only Savitzky-Golay filter utility
├── src/
│   ├── training.cpp        # Training implementation
│   ├── prediction.cpp      # Online prediction implementation
│   └── main.cpp            # CLI entry point (promp_rt)
└── lib/
    └── promp_1/            # Bundled ProMP_1 library (no install needed)
```

---

## Dependencies

| Library | Min version | How to obtain |
|---|---|---|
| **CMake** | 3.17 | [cmake.org](https://cmake.org/download/) or package manager |
| **C++ compiler** | GCC ≥ 9 / Clang ≥ 10 / MSVC 2019 | System compiler |
| **Eigen3** | 3.3 | Install via package manager — see below |
| **ProMP_1** | — | **Bundled** in `lib/promp_1/`, no install needed |

**The only library you need to install is Eigen3.**

---

## Installing Eigen3

Eigen3 is a header-only linear algebra library.

### Ubuntu / Debian

```bash
sudo apt install libeigen3-dev
```

### macOS (Homebrew)

```bash
brew install eigen
```

### Windows — vcpkg

```powershell
# Install vcpkg if not already present
git clone https://github.com/microsoft/vcpkg.git C:\vcpkg
C:\vcpkg\bootstrap-vcpkg.bat

# Install Eigen3
C:\vcpkg\vcpkg install eigen3:x64-windows
```

Then pass the vcpkg toolchain file to CMake:

```powershell
cmake .. -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake
```

### Windows — conda / mamba

```bash
conda install -c conda-forge eigen
```

### Custom location

If you have Eigen3 at a non-standard path (the directory that contains the
`Eigen/` folder):

```bash
cmake .. -DEIGEN3_ROOT=/path/to/eigen-3.x.x
```

---

## Build

```bash
# Ubuntu / Debian — install Eigen3 first (once):
sudo apt install libeigen3-dev

cd promp_realtime
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

The executable `promp_rt` is produced inside `build/`.

### macOS

```bash
brew install eigen    # once

cd promp_realtime
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(sysctl -n hw.logicalcpu)
```

### Windows (Visual Studio)

```powershell
cd promp_realtime
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release `
         -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build . --config Release --parallel
```

---

### CMake variables

| Variable | Default | Description |
|---|---|---|
| `EIGEN3_ROOT` | _(empty)_ | Path to a custom Eigen3 directory (the parent of the `Eigen/` folder). Leave empty to use the system installation. |
| `PROMP1_SRC_DIR` | `lib/promp_1` | Path to the bundled promp_1 source. Override only if you move the library. |

---

## Demo CSV format

Each demonstration is a **single CSV file**.  The required columns are:

```
time_s, pos0_rad, pos1_rad, ..., pos(n_dof-1)_rad
```

Example for 3-DoF:

```
time_s,pos0_rad,pos1_rad,pos2_rad
0.000,0.1234,0.4500,0.0100
0.010,0.1280,0.4420,0.0115
...
```

| Column | Description |
|---|---|
| `time_s` | Elapsed time in seconds (monotonically increasing) |
| `pos0_rad … pos(n_dof-1)_rad` | Joint positions in radians |

**Velocities are not required** in the CSV — they are computed internally
from positions using the Savitzky-Golay filter.

An optional header row is auto-detected and skipped.
All `*.csv` files in the demo folder are loaded.
Different demonstrations may have different lengths or sampling rates.

---

## Internal trajectory representation

After loading, each demonstration is converted to a **`2×n_dof`-dimensional**
trajectory with interleaved position and velocity channels:

```
dim 0: pos0 (rad)     dim 1: vel0 (rad/s)
dim 2: pos1 (rad)     dim 3: vel1 (rad/s)
...
dim 2k: pos_k (rad)   dim 2k+1: vel_k (rad/s)
```

The ProMP is trained on these `2×n_dof`-D trajectories.
At inference, each sensor observation is augmented with an online
velocity estimate (from the SG filter) before Bayesian conditioning.

---

## Usage

### 1. Training

```bash
promp_rt train <demo_folder> <model_dir> [options]
```

**Options:**

| Flag | Default | Description |
|---|---|---|
| `--dof N` | 3 | Number of degrees of freedom (joints) in the input CSV |
| `--basis N` | 10 | Number of RBF basis functions. Increase for complex shapes. |
| `--std S` | -1 (auto) | RBF std-dev; ≤ 0 = auto-tuned to `1/(basis²)` |
| `--steps N` | 100 | Common trajectory length after resampling |
| `--max-duration D` | -1 (all) | Clip demonstrations beyond D seconds |
| `--sg-window W` | 9 | SG filter window length (odd). Larger = smoother velocity. |
| `--sg-poly P` | 4 | SG polynomial order. Must satisfy P < W. |
| `--obs-pos-noise V` | 1e-4 | Position observation noise variance for inference conditioning |
| `--obs-vel-noise V` | 1e-2 | Velocity observation noise variance (relaxed; SG estimation error) |
| `--goal-pos-noise V` | 1e-6 | Goal position noise variance (tight constraint) |
| `--goal-vel-noise V` | 1e-2 | Goal velocity noise variance (soft; natural stop assumed) |

**Example (3-DoF):**

```bash
promp_rt train ./demos ./model --dof 3 --basis 10 --steps 100 --max-duration 1.4
```

**Example (6-DoF):**

```bash
promp_rt train ./demos ./model --dof 6 --basis 15 --sg-window 11
```

**Output** in `./model/`:
- `model.promp` — serialised ProMP (weight mean μ_w, covariance Σ_w, basis parameters)
- `metadata.txt` — training configuration and mean demo duration (loaded automatically at inference)

---

### 2. Prediction (real-time streaming)

```bash
promp_rt predict <model_dir> --mode <mode> [options]
```

The predictor reads `n_dof` from `metadata.txt` automatically — no `--dof`
flag is needed at inference time.

**Conditioning modes:**

| Mode | Description |
|---|---|
| `past_traj` | Condition on all past observations (positions + velocities) |
| `past_traj_target` | Also condition on a known final target position |

**Options:**

| Flag | Default | Description |
|---|---|---|
| `--mode` | `past_traj` | Conditioning mode |
| `--target t0 t1 … t(n-1)` | zeros | Final target positions (rad) for `past_traj_target`; must supply exactly `n_dof` values |
| `--duration D` | (from training) | Override expected motion duration (s) |
| `--future-steps N` | 100 | Resolution of the returned future trajectory |
| `--max-obs N` | 0 (all) | Limit past observations per prediction tick (0 = all) |

**Input format (stdin):** one observation per line:

```
<time_s> <pos0_rad> <pos1_rad> … <pos(n_dof-1)_rad>
```

**Output format (stdout):** one block per observation:

```
PRED <n_points> <current_phase>
<time_s> <pos0_mean> <pos0_std> <vel0_mean> <vel0_std>  <pos1_mean> …
...
```

stdout is flushed after each block for pipe-compatible real-time use.

**Example – past trajectory only (3-DoF):**

```bash
promp_rt predict ./model --mode past_traj --future-steps 50
```

**Example – past trajectory + target (3-DoF):**

```bash
promp_rt predict ./model --mode past_traj_target \
    --target 0.0 -0.5 0.0 --duration 1.4
```

**Example – piping from a live sensor:**

```bash
robot_sensor | promp_rt predict ./model --mode past_traj_target \
    --target 0.0 -0.5 0.0
```

---

## Embedding the API in your application

```cpp
#include "training.hpp"
#include "prediction.hpp"

// ── Training (once, offline) ─────────────────────────────────────────────
promp_rt::TrainingConfig cfg;
cfg.n_dof               = 3;    // number of joints
cfg.num_basis_functions = 10;
cfg.sg_window           = 9;    // Savitzky-Golay window for velocity
cfg.sg_poly_order       = 4;
promp_rt::train_from_folder("./demos", "./model", cfg);

// ── Prediction (online, at each sensor tick) ─────────────────────────────
promp_rt::OnlinePredictor predictor("./model", /*max_obs=*/0);
predictor.reset(/*expected_duration_s=*/1.4);

int n = predictor.get_n_dof();   // read n_dof from the loaded model

// Sensor loop:
while (running) {
    double t = get_elapsed_time_s();
    std::vector<double> pos(n);
    for (int i = 0; i < n; ++i) pos[i] = read_joint_rad(i);

    std::vector<double> target = {0.0, -0.5, 0.0};   // one value per joint
    auto result = predictor.update_and_predict(
        t, pos,
        promp_rt::ConditioningMode::PAST_TRAJ_TARGET,
        /*n_future_steps=*/100,
        target);

    if (result.valid) {
        for (int step = 0; step < (int)result.future_times_s.size(); ++step) {
            double t_pred = result.future_times_s[step];
            for (int dof = 0; dof < n; ++dof) {
                double pos_mean = result.mean_pos(step, dof);
                double vel_mean = result.mean_vel(step, dof);
                double pos_sd   = result.std_pos(step, dof);
            }
        }
    }
}
```

---

## Conditioning modes — mathematical detail

### `PAST_TRAJ`

Each accumulated observation `(t_k, pos_k, vel_k)` is mapped to phase step
`s_k = round(t_k / T_expected × (N−1))` and applied as a Bayesian via-point:

```
via_k = [pos0_k, vel0_k, pos1_k, vel1_k, ..., pos(n-1)_k, vel(n-1)_k]  ∈ ℝ^(2n)
Σ_obs = diag(σ²_pos, σ²_vel, σ²_pos, σ²_vel, ...)

L    = Σ_w Φ(s_k) (Σ_obs + Φ(s_k)ᵀ Σ_w Φ(s_k))⁻¹
μ_w ← μ_w + L(via_k − Φ(s_k)ᵀ μ_w)
Σ_w ← Σ_w − L Φ(s_k)ᵀ Σ_w
```

### `PAST_TRAJ_TARGET`

Same as above, plus a goal condition at the final step:

```
goal = [pos0_target, 0, pos1_target, 0, ..., pos(n-1)_target, 0]  ∈ ℝ^(2n)
Σ_goal = diag(σ²_goal_pos, σ²_goal_vel, ...)
```

Position channels use the tight `goal_pos_noise_var`; velocity channels
use the soft `goal_vel_noise_var` (assumes the motion naturally comes to rest).

The trained model is **never modified** — a fresh copy is created for each
`predict()` call.

---

## Tuning guide

| Parameter | Effect | Start here |
|---|---|---|
| `--dof N` | Joints in data | Match your sensor |
| `--basis N` | Model capacity | 8–15 for smooth motions, up to 20 for complex |
| `--sg-window W` | Velocity smoothness | 7–11 (odd); scale with sensor rate |
| `--sg-poly P` | Velocity peak accuracy | 3–4, must be < W |
| `--obs-pos-noise V` | Tightness of position conditioning | 1e-5 – 1e-3 |
| `--obs-vel-noise V` | Trust in estimated velocities | 1e-3 – 0.1 |
| `--goal-pos-noise V` | Strictness of target constraint | 1e-7 – 1e-5 |
| `--max-obs N` | Per-tick computation budget | 20–50 for high-rate sensors |
| `--duration D` | Phase alignment | Set to actual motion duration |

### Savitzky-Golay window vs. sensor rate

```
sg_window ≈ data_rate_Hz × desired_smoothing_time_s   (round to odd)
```

| Sensor rate | Smoothing time | Recommended window |
|---|---|---|
| 50 Hz | 90 ms | 5 |
| 100 Hz | 90 ms | 9 (default) |
| 200 Hz | 90 ms | 19 |
| 500 Hz | 90 ms | 45 |
| 1000 Hz | 90 ms | 91 |

### Data-frequency compatibility

The predictor is fully rate-agnostic:

- **Phase mapping** is time-ratio based: `phase = elapsed_s / duration_s` — no sample rate assumption.
- **SG velocity filter** uses actual timestamps — works with non-uniform or variable-rate data.
- **Resampling** (training) uses linear interpolation on normalised phase, so the training grid is always uniform regardless of input rate.
- **Mixed rates** across demonstrations are supported; each demo is processed independently before merging.

Only `sg_window` needs to be scaled relative to the sensor rate (see table above).
