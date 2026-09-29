// Single-sensor load-cell reader for Windows.
//
// Reads the 8-byte binary frames sent by load_cell_one_sensor.ino,
// applies tare + fractional decimation + optional EMA, and prints (and
// optionally logs) force in Newtons.
//
// Build (MinGW / MSYS2):
//     g++ -std=c++17 -O2 read_force.cpp -o read_force.exe
//
// Build (MSVC, "x64 Native Tools" prompt):
//     cl /std:c++17 /O2 /EHsc read_force.cpp
//
// Linux/macOS: needs porting to termios in place of the WinAPI section
// -- the frame parsing, tare, and decimation are portable.

#include <windows.h>

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <deque>
#include <fstream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

// ============================================================================
//  USER CONFIG
// ============================================================================

// Serial link  (must match the Arduino sketch)
static const char* PORT            = "COM3";      // "COM10" and up: OK either way
static const int   BAUD            = 250000;
static const int   ARDUINO_RATE_HZ = 1000;

// The knob you'll actually tweak: any value <= ARDUINO_RATE_HZ.
static const double READ_RATE_HZ = 100.0;

// Calibration:  mass_g = A*raw + B,   force_N = mass_g * G
static const double A = 19.96;
static const double B = 10.32;
static const double G = 0.00981;

// Tare
static const int TARE_WARMUP  = 500;   // raw samples discarded first
static const int TARE_SAMPLES = 1000;  // raw samples averaged after warmup

// EMA filter on emitted outputs (1.0 = off)
static const double FILTER_ALPHA = 1.0;

// Logging / display
static const bool   LOG_TO_FILE    = false;
static const double PRINT_INTERVAL = 0.1;   // seconds between console refreshes
static const double RATE_ALPHA     = 0.3;   // EMA on displayed rate

// Live plot (requires gnuplot; see comment near GnuplotLive)
static const bool   USE_PLOT       = true;
static const double PLOT_WINDOW_S  = 10.0;  // rolling x-axis width, seconds
static const double PLOT_INTERVAL  = 0.2;   // seconds between plot redraws

// ============================================================================
//  PROTOCOL CONSTANTS  --  must match the Arduino sketch
// ============================================================================

static const uint8_t HEADER0    = 0xAA;
static const uint8_t HEADER1    = 0x55;
//static const int     FRAME_SIZE = 8;   // 2 hdr + 4 counter + 2 raw
static const int FRAME_SIZE = 10;

// ============================================================================

using clk = std::chrono::steady_clock;

static double now_s() {
    return std::chrono::duration<double>(clk::now().time_since_epoch()).count();
}

static double raw_to_newton(uint16_t raw) {
    return (A * (double)raw + B) * G;
}

// Open a Windows COM port. Returns INVALID_HANDLE_VALUE on failure.
static HANDLE open_serial(const char* port, int baud) {
    // COM10 and above require the "\\.\COMx" prefix to work with CreateFile.
    std::string name = port;
    if (name.rfind("\\\\.\\", 0) != 0) name = std::string("\\\\.\\") + name;

    HANDLE h = CreateFileA(name.c_str(),
                           GENERIC_READ | GENERIC_WRITE,
                           0, nullptr, OPEN_EXISTING, 0, nullptr);
    if (h == INVALID_HANDLE_VALUE) return h;

    DCB dcb{};
    dcb.DCBlength = sizeof(dcb);
    if (!GetCommState(h, &dcb)) { CloseHandle(h); return INVALID_HANDLE_VALUE; }
    dcb.BaudRate = baud;
    dcb.ByteSize = 8;
    dcb.Parity   = NOPARITY;
    dcb.StopBits = ONESTOPBIT;
    dcb.fBinary  = TRUE;
    dcb.fDtrControl = DTR_CONTROL_ENABLE;   // triggers Arduino auto-reset
    if (!SetCommState(h, &dcb)) { CloseHandle(h); return INVALID_HANDLE_VALUE; }

    COMMTIMEOUTS to{};
    to.ReadIntervalTimeout        = MAXDWORD;
    to.ReadTotalTimeoutConstant   = 100;    // ms
    to.ReadTotalTimeoutMultiplier = 0;
    SetCommTimeouts(h, &to);

    SetupComm(h, 65536, 4096);              // larger RX buffer helps at 1 kHz
    return h;
}

// ----------------------------------------------------------------------------
//  Live plot via a gnuplot pipe.
//
//  Prerequisite: gnuplot installed (see GNUPLOT_EXE below).
//    Windows: install from https://sourceforge.net/projects/gnuplot/
//    Verify:  run `"C:\Program Files\gnuplot\bin\gnuplot.exe" --version`.
//
//  We invoke gnuplot by its absolute path rather than relying on PATH: a
//  PATH change made via the registry only reaches processes started after
//  the change, so a shell/IDE opened earlier wouldn't see a PATH-only fix.
//
//  Redraws run on their own thread at PLOT_INTERVAL, decoupled from the
//  serial-read loop via a mutex-protected snapshot. Without this, a redraw
//  that gnuplot can't keep up with blocks the fprintf/fflush call, which
//  would stall serial reads too -- the whole program looked "frozen".
// ----------------------------------------------------------------------------
static const char* GNUPLOT_EXE = "C:\\Program Files\\gnuplot\\bin\\gnuplot.exe";

class GnuplotLive {
public:
    GnuplotLive(double window_s, double interval_s)
        : window_s_(window_s), interval_s_(interval_s) {
        std::string cmd = std::string("\"") + GNUPLOT_EXE + "\"";
        pipe_ = _popen(cmd.c_str(), "w");
        if (!pipe_) {
            fprintf(stderr,
                    "Warning: could not launch gnuplot -- live plot disabled.\n"
                    "         Expected it at %s\n", GNUPLOT_EXE);
            return;
        }
        // Native win32 GDI terminal -- much lighter than the wxt/qt terminals,
        // which struggled to keep up and caused the pipe writes to block.
        fprintf(pipe_, "set terminal windows\n");
        fprintf(pipe_, "set title 'Force (N)'\n");
        fprintf(pipe_, "set xlabel 't (s)'\n");
        fprintf(pipe_, "set ylabel 'F (N)'\n");
        fprintf(pipe_, "set grid\n");
        fflush(pipe_);
        worker_ = std::thread(&GnuplotLive::run, this);
    }

    ~GnuplotLive() {
        if (pipe_) {
            stop_ = true;
            if (worker_.joinable()) worker_.join();
            fprintf(pipe_, "exit\n");
            _pclose(pipe_);
        }
    }

    bool ok() const { return pipe_ != nullptr; }

    // Called from the acquisition loop -- just appends under a short lock,
    // never touches the gnuplot pipe, so it can't stall serial reads.
    void push(double t, double f) {
        if (!pipe_) return;
        std::lock_guard<std::mutex> lk(mu_);
        data_.push_back({t, f});
        while (!data_.empty() && data_.front().t < t - window_s_) {
            data_.pop_front();
        }
        latest_t_ = t;
    }

private:
    struct Pt { double t, f; };

    void run() {
        while (!stop_) {
            std::this_thread::sleep_for(std::chrono::duration<double>(interval_s_));

            std::vector<Pt> snapshot;
            double t_now;
            {
                std::lock_guard<std::mutex> lk(mu_);
                if (data_.size() < 2) continue;
                snapshot.assign(data_.begin(), data_.end());
                t_now = latest_t_;
            }

            fprintf(pipe_, "set xrange [%f:%f]\n", t_now - window_s_, t_now);
            fprintf(pipe_, "plot '-' with lines lw 1.5 title 'F'\n");
            for (const auto& p : snapshot) fprintf(pipe_, "%f %f\n", p.t, p.f);
            fprintf(pipe_, "e\n");
            fflush(pipe_);
        }
    }

    std::deque<Pt>    data_;
    double            window_s_;
    double            interval_s_;
    double            latest_t_ = 0.0;
    FILE*             pipe_ = nullptr;
    std::mutex        mu_;
    std::thread       worker_;
    std::atomic<bool> stop_{false};
};

int main() {
    printf("Opening %s at %d baud...\n", PORT, BAUD);
    HANDLE h = open_serial(PORT, BAUD);
    if (h == INVALID_HANDLE_VALUE) {
        fprintf(stderr,
                "Failed to open %s (error %lu). Wrong port name, or is the\n"
                "Arduino IDE Serial Monitor / another program still holding it?\n",
                PORT, (unsigned long)GetLastError());
        return 1;
    }

    // Wait for the Arduino to finish auto-reset after opening the port.
    std::this_thread::sleep_for(std::chrono::seconds(2));
    PurgeComm(h, PURGE_RXCLEAR | PURGE_TXCLEAR);

    // Fractional decimation setup
    double read_hz            = READ_RATE_HZ;
    if (read_hz > ARDUINO_RATE_HZ) {
        printf("Warning: READ_RATE_HZ (%g) > ARDUINO_RATE_HZ (%d); capping.\n",
               read_hz, ARDUINO_RATE_HZ);
        read_hz = ARDUINO_RATE_HZ;
    }
    const double samples_per_output = (double)ARDUINO_RATE_HZ / read_hz;
    printf("Reading at %g Hz (~%.3f Arduino samples per output).\n",
           read_hz, samples_per_output);

    // CSV log
    std::ofstream log;
    if (LOG_TO_FILE) {
        std::time_t t = std::time(nullptr);
        char fname[64];
        std::strftime(fname, sizeof(fname),
                      "force_%Y%m%d_%H%M%S.csv", std::localtime(&t));
        log.open(fname);
        log << "t_s,F_N,counter\n";
        printf("Logging to %s\n", fname);
    }

    printf("Warming up (%d samples) before taring...\n", TARE_WARMUP);

    // Live plot (optional)
    GnuplotLive plot(PLOT_WINDOW_S, PLOT_INTERVAL);
    const bool plot_enabled = USE_PLOT && plot.ok();

    // --- State -------------------------------------------------------------
    std::vector<double> tare_buf;
    std::vector<double> tare_buf2;
    tare_buf.reserve(TARE_SAMPLES);
    tare_buf2.reserve(TARE_SAMPLES);

    int      tare_warmup_left  = TARE_WARMUP;
    bool     tared             = false;
    double   tare_offset       = 0.0;
    double   tare_offset2 = 0.0;

    bool     filt_init         = false;
    double   filt              = 0.0;
    double filt2 = 0.0;

    uint32_t first_counter     = 0;
    uint32_t last_counter      = 0;
    bool     have_last_counter = false;
    uint64_t n_drops           = 0;

    double   acc               = 0.0;
    double acc2 = 0.0;

    int      n_in_block        = 0;
    uint64_t n_since_start     = 0;
    double   next_emit_at      = samples_per_output;

    double   next_print        = now_s() + PRINT_INTERVAL;
    int      n_out_since_print = 0;
    double   rate_hz           = 0.0;

    std::vector<uint8_t> buf;
    buf.reserve(4096);
    uint8_t chunk[1024];

    // --- Main loop ---------------------------------------------------------
    while (true) {
        DWORD got = 0;
        if (!ReadFile(h, chunk, sizeof(chunk), &got, nullptr)) {
            fprintf(stderr, "\nSerial read error.\n");
            break;
        }
        if (got > 0) buf.insert(buf.end(), chunk, chunk + got);

        // Drain complete frames from the buffer.
        size_t pos = 0;
        while (buf.size() - pos >= (size_t)FRAME_SIZE) {
            if (buf[pos] != HEADER0 || buf[pos + 1] != HEADER1) {
                ++pos;                          // resync, one byte at a time
                continue;
            }
            const uint8_t* p = buf.data() + pos;
            uint32_t counter =  (uint32_t)p[2]
                              | ((uint32_t)p[3] <<  8)
                              | ((uint32_t)p[4] << 16)
                              | ((uint32_t)p[5] << 24);
            uint16_t raw     =  (uint16_t)p[6] | ((uint16_t)p[7] << 8);
            uint16_t raw2 = (uint16_t)p[8] | ((uint16_t)p[9] << 8);

            pos += FRAME_SIZE;

            double f_N = raw_to_newton(raw);
            double f2 = raw_to_newton(raw2);

            // Drop detection ------------------------------------------------
            if (have_last_counter) {
                uint32_t expected = last_counter + 1;
                if (counter != expected) {
                    uint32_t missed = counter - expected;
                    n_drops += missed;
                    printf("!! Missed %u sample(s) (counter jumped %u -> %u)\n",
                           (unsigned)missed, (unsigned)expected,
                           (unsigned)counter);
                }
            }
            last_counter      = counter;
            have_last_counter = true;

            // Tare phase ---------------------------------------------------
            if (!tared) {
                if (tare_warmup_left > 0) {
                    if (--tare_warmup_left == 0) {
                        printf("Taring: averaging %d unloaded samples...\n",
                               TARE_SAMPLES);
                    }
                    continue;
                }
                tare_buf.push_back(f_N);
                tare_buf2.push_back(f2);
                if ((int)tare_buf.size() >= TARE_SAMPLES) {
                    double m = 0.0;
                    for (double v : tare_buf) m += v;
                    m /= tare_buf.size();
                    double s2 = 0.0;
                    for (double v : tare_buf) s2 += (v - m) * (v - m);
                    double s = std::sqrt(s2 / tare_buf.size());
                    tare_offset = m;

                    m = 0.0;
                    for (double v : tare_buf2) m += v;
                    m /= tare_buf2.size();
                    s2 = 0.0;
                    for (double v : tare_buf2) s2 += (v - m) * (v - m);
                    s = std::sqrt(s2 / tare_buf2.size());
                    tare_offset2 = m;


                    tared       = true;
                    printf("Tared: F offset = %+.4f N   (std = %.4f N)\n",
                           tare_offset, s);
                    first_counter = counter + 1;
                    next_print    = now_s() + PRINT_INTERVAL;
                    printf("%8s  %10s  %10s  %10s  %6s\n",
                           "t (s)", "F1 (N)", "F2 (N)", "rate (Hz)", "drops");
                }
                continue;
            }

            // Live sample: accumulate then emit ----------------------------
            acc            += f_N;
            acc2 += f2;
            n_in_block     += 1;
            n_since_start  += 1;
            if ((double)n_since_start < next_emit_at) continue;

            double f_avg = acc / n_in_block - tare_offset;
            double f_avg2 = acc2 / n_in_block - tare_offset2;
            next_emit_at += samples_per_output;
            acc          = 0.0;
            acc2 = 0.0;
            n_in_block   = 0;

            if (!filt_init) { filt = f_avg; filt2 = f_avg2; filt_init = true; }
            else
            {
                filt += FILTER_ALPHA * (f_avg - filt);
                filt2 += FILTER_ALPHA * (f_avg2 - filt2);
            }

            double t = (double)(counter - first_counter) / ARDUINO_RATE_HZ;

            if (log.is_open()) {
                char line[128];
                std::snprintf(line, sizeof(line), "%.6f,%.4f,%.4f,%u\n",
                              t, filt, filt2, (unsigned)counter);
                log << line;
            }

            if (plot_enabled) plot.push(t, filt);

            ++n_out_since_print;
            double now = now_s();
            if (now >= next_print) {
                double inst = n_out_since_print / PRINT_INTERVAL;
                rate_hz = (rate_hz == 0.0)
                            ? inst
                            : (1.0 - RATE_ALPHA) * rate_hz + RATE_ALPHA * inst;
                printf("%8.2f  %10.3f  %10.3f  %10.1f  %6llu\n",
                       t, filt, filt2, rate_hz, (unsigned long long)n_drops);
                std::fflush(stdout);
                next_print       += PRINT_INTERVAL;
                n_out_since_print = 0;
            }
        }
        // Compact: drop the bytes we've consumed.
        if (pos > 0) buf.erase(buf.begin(), buf.begin() + pos);
    }

    CloseHandle(h);
    return 0;
}
