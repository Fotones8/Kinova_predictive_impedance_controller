//---------------------------------------------------------
// The main function for impedance controller
//---------------------------------------------------------
// Description: The main function for the impedance controller used in the PHRI project
// Copyright: Yihan Liu 2024
//---------------------------------------------------------
#define _HAS_STD_BYTE 0

#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <fstream>
#include <thread>
#include <iomanip>
#include <atomic>
#include <chrono>
#include <cmath>
#include <mutex>


#include <KDetailedException.h>

#include <BaseClientRpc.h>
#include <BaseCyclicClientRpc.h>
#include <ActuatorConfigClientRpc.h>
#include <SessionClientRpc.h>
#include <SessionManager.h>

#include <RouterClient.h>
#include <TransportClientTcp.h>
#include <TransportClientUdp.h>

#include <google/protobuf/util/json_util.h>


#if defined(_MSC_VER)
#include <Windows.h>
#else
#include <unistd.h>
#endif
#include <time.h>

#include <Eigen/Dense>
#include <Jacobian.h>
#include <cmath>
#include <Fwd_kinematics.h>
#include <Dynamics.h>
#include <Controller.h>
#include <Filter.h>
#include <conio.h>

#include "training.hpp"
#include "prediction.hpp"
#include "ukf_predictor.hpp"

#include <timeapi.h>

namespace k_api = Kinova::Api;
using namespace Jacobian;
using namespace Fwd_kinematics;
using namespace Dynamics;
using namespace Controller;
using namespace Filter;

#define IP_ADDRESS "192.168.1.10"
#define PORT 10000
#define PORT_REAL_TIME 10001
#define ACTUATOR_COUNT 7
#define CONTROL_FREQUENCY 450

// FORCE SENSORS
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
static const char* SPORT            = "COM3";      // "COM10" and up: OK either way
static const int   BAUD            = 250000;
static const int   ARDUINO_RATE_HZ = 1000;

// The knob you'll actually tweak: any value <= ARDUINO_RATE_HZ.
static const double READ_RATE_HZ = 400.0;

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
static const double PRINT_INTERVAL = 0.025;   // seconds between console refreshes
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
    to.ReadTotalTimeoutConstant   = 1;    // ms
    to.ReadTotalTimeoutMultiplier = 0;
    SetCommTimeouts(h, &to);

    DWORD latency = 1;  // ms
    DWORD bytes_returned;
    DeviceIoControl(h,
        0x00166A08,          // IOCTL_SERIAL_SET_LATENCY (FTDI-specific)
        &latency, sizeof(latency),
        nullptr, 0,
        &bytes_returned, nullptr);

    SetupComm(h, 65536, 4096);              // larger RX buffer helps at 1 kHz
    return h;
}

//  SHARED SENSOR DATA
//  Declared globally so all threads can access it.

struct SensorReading {
    double force_N = 0.0;
};

// The shared variable
static SensorReading g_sensor;

// The mutex that protects it — always lock this before reading or writing
static std::mutex    g_sensor_mu;

// Optional: a lightweight flag to tell the sensor thread to stop cleanly
static std::atomic<bool> g_sensor_stop{false};

void sensor_thread_func()
{
    timeBeginPeriod(1);
    // ── YOUR SENSOR INITIALISATION CODE GOES HERE ────────────────────────
    // e.g. open serial port, tare, warm up, etc.

        printf("Opening %s at %d baud...\n", SPORT, BAUD);
    HANDLE h = open_serial(SPORT, BAUD);
    if (h == INVALID_HANDLE_VALUE) {
        fprintf(stderr,
                "Failed to open %s (error %lu). Wrong port name, or is the\n"
                "Arduino IDE Serial Monitor / another program still holding it?\n",
                SPORT, (unsigned long)GetLastError());
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
    /*
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
    */

    printf("Warming up (%d samples) before taring...\n", TARE_WARMUP);

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
                    //printf("%8s  %10s  %10s  %10s  %6s\n",
                    //       "t (s)", "F1 (N)", "F2 (N)", "rate (Hz)", "drops");
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

            /*
            if (log.is_open()) {
                char line[128];
                std::snprintf(line, sizeof(line), "%.6f,%.4f,%.4f,%u\n",
                              t, filt, filt2, (unsigned)counter);
                log << line;
            }
            */

            {
                std::lock_guard<std::mutex> lk(g_sensor_mu);
                g_sensor.force_N = filt+filt2;
            }
            /*
            ++n_out_since_print;
            double now = now_s();
            if (now >= next_print) {
                double inst = n_out_since_print / PRINT_INTERVAL;
                rate_hz = (rate_hz == 0.0)
                            ? inst
                            : (1.0 - RATE_ALPHA) * rate_hz + RATE_ALPHA * inst;
                //printf("%8.2f  %10.3f  %10.3f  %10.1f  %6llu\n",
                //       t, filt, filt2, rate_hz, (unsigned long long)n_drops);
                //std::fflush(stdout);
                next_print       += PRINT_INTERVAL;
                n_out_since_print = 0;


            }
            */
        }
        // Compact: drop the bytes we've consumed.
        if (pos > 0) buf.erase(buf.begin(), buf.begin() + pos);
    }
    timeEndPeriod(1);
    CloseHandle(h);

    // ── YOUR SENSOR CLEANUP CODE GOES HERE ───────────────────────────────
    // e.g. close serial port, close file handles, etc.
}


double read_sensor()
{
    // Always take a local copy under the lock, then use the copy outside.
    // Never hold the mutex while doing heavy computation.
    SensorReading local;
    {
        std::lock_guard<std::mutex> lk(g_sensor_mu);
        local = g_sensor;
    }

    return local.force_N;
}

namespace {

constexpr UINT WM_STATUS_UPDATE = WM_APP + 1;

struct StatusUpdate {
    std::wstring message;
    COLORREF bgColor;
    COLORREF textColor;
};

std::atomic<HWND> g_hwnd{nullptr};
std::wstring g_message = L"Starting the system";
COLORREF g_bgColor   = RGB(255, 255, 255);
COLORREF g_textColor = RGB(0, 0, 0);
HBRUSH g_bgBrush = nullptr;

LRESULT CALLBACK StatusWndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    switch (msg) {
        case WM_STATUS_UPDATE: {
            // Runs on the window's own thread (messages are only ever
            // dispatched here, never called directly from another thread).
            auto* update = reinterpret_cast<StatusUpdate*>(lParam);

            g_message   = update->message;
            g_bgColor   = update->bgColor;
            g_textColor = update->textColor;
            delete update; // we own it; caller heap-allocated it for us

            if (g_bgBrush) DeleteObject(g_bgBrush);
            g_bgBrush = CreateSolidBrush(g_bgColor);
            SetClassLongPtrW(hwnd, GCLP_HBRBACKGROUND, (LONG_PTR)g_bgBrush);

            InvalidateRect(hwnd, nullptr, TRUE);
            return 0;
        }
        case WM_PAINT: {
            PAINTSTRUCT ps;
            HDC hdc = BeginPaint(hwnd, &ps);

            RECT rect;
            GetClientRect(hwnd, &rect);
            FillRect(hdc, &rect, g_bgBrush);

            SetTextColor(hdc, g_textColor);
            SetBkMode(hdc, TRANSPARENT);

            HFONT font = CreateFontW(128, 0, 0, 0, FW_BOLD, FALSE, FALSE, FALSE,
                                      DEFAULT_CHARSET, OUT_DEFAULT_PRECIS,
                                      CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY,
                                      DEFAULT_PITCH, L"Arial");
            HFONT oldFont = (HFONT)SelectObject(hdc, font);

            DrawTextW(hdc, g_message.c_str(), -1, &rect,
                      DT_CENTER | DT_VCENTER | DT_SINGLELINE);

            SelectObject(hdc, oldFont);
            DeleteObject(font);
            EndPaint(hwnd, &ps);
            return 0;
        }
        case WM_CLOSE:
            // Ignore the user clicking the window's own close button — this
            // is meant to be a persistent status display for the duration
            // of the program. Remove this case if you want it closable.
            return 0;
        case WM_DESTROY:
            PostQuitMessage(0);
            return 0;
        default:
            return DefWindowProcW(hwnd, msg, wParam, lParam);
    }
}

void status_window_thread_main() {
    g_bgBrush = CreateSolidBrush(g_bgColor);

    const wchar_t* className = L"PersistentStatusWindow";
    WNDCLASSW wc = {};
    wc.lpfnWndProc   = StatusWndProc;
    wc.hInstance     = GetModuleHandleW(nullptr);
    wc.lpszClassName = className;
    wc.hbrBackground = g_bgBrush;
    RegisterClassW(&wc);

    HWND hwnd = CreateWindowExW(
        WS_EX_TOPMOST, className, L"Status", WS_OVERLAPPEDWINDOW & ~WS_THICKFRAME,
        CW_USEDEFAULT, CW_USEDEFAULT, 2500, 1000,
        nullptr, nullptr, wc.hInstance, nullptr);

    g_hwnd.store(hwnd);

    ShowWindow(hwnd, SW_SHOW);
    UpdateWindow(hwnd);

    MSG msg;
    while (GetMessageW(&msg, nullptr, 0, 0)) {
        TranslateMessage(&msg);
        DispatchMessageW(&msg);
    }
}

} // namespace

// ----------------------------------------------------------------------------
// start_status_window — call ONCE, early in main(), before other threads
// start pushing updates. Spawns the window's owning thread and detaches it;
// the window lives for the rest of the process.
// ----------------------------------------------------------------------------
void start_status_window() {
    std::thread(status_window_thread_main).detach();
    // Give the window a moment to actually exist before anyone posts to it.
    // (A more robust version could use a condition_variable instead of a
    // fixed sleep, but a few ms is fine for a status display at startup.)
    while (g_hwnd.load() == nullptr) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
}

// ----------------------------------------------------------------------------
// update_status_window — call from ANY thread, ANY time. Non-blocking:
// PostMessageW queues the update and returns immediately; the actual
// redraw happens asynchronously on the window's own thread.
// ----------------------------------------------------------------------------
void update_status_window(const std::wstring& message, COLORREF bgColor, COLORREF textColor) {
    HWND hwnd = g_hwnd.load();
    if (!hwnd) return; // start_status_window() wasn't called, or hasn't finished yet

    auto* update = new StatusUpdate{message, bgColor, textColor};
    PostMessageW(hwnd, WM_STATUS_UPDATE, 0, reinterpret_cast<LPARAM>(update));
}


// FLAGS FOR ALGORITHM CONTROL
std::atomic<int> flag_subject_number{0};
std::atomic<int> flag_mode{1};
std::atomic<int> flag_movement{1};
std::atomic<int> flag_goal{0};
std::atomic<int> flag_running{0};
std::atomic<int> flag_reset{0};
std::atomic<int> flag_initial_position{0};


// ----------------------------------------------------------------------------
// run_single_movement — everything from "go to initial position" through
// "target reached" for ONE trial. Shared by both manual mode (called once
// per outer-loop iteration) and auto mode (called 5x in a row per iteration),
// so the initial-position wait / reset / "press any button to start" /
// running wait / status messages only need to be written once.
// ----------------------------------------------------------------------------
std::wstring message = L"Do not hold,\n waiting for calibration";
string input;

void run_single_movement() {
    message = L"Move to the initial position";
    update_status_window(message, RGB(0, 0, 0), RGB(255, 255, 255));

    flag_initial_position.store(0);
    while (flag_initial_position.load() == 0) {
        // Waiting, displaying GO
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }

    // Reset the system to start new movement
    flag_reset.store(1);
    message = L"Prepare to go to target " + std::to_wstring(flag_goal.load());
    update_status_window(message, RGB(0, 0, 200), RGB(255, 255, 255));

    // Kept exactly as before, in both manual and auto mode: still a real
    // console prompt each movement, even though mode/movement/goal aren't
    // re-asked in auto mode.
    std::cout << "Press any button to start: \n";
    cin >> input;
    flag_running.store(1);

    message = L"START MOVING TO " + std::to_wstring(flag_goal.load());
    update_status_window(message, RGB(0, 200, 0), RGB(255, 255, 255));
    while (flag_running.load() == 1) {
        // Waiting, displaying GO
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }

    message = L"Target reached, reset to initial position";
    update_status_window(message, RGB(173, 216, 230), RGB(255, 255, 255));
}

// ----------------------------------------------------------------------------
// ask_yes_no — small helper so the auto-mode prompt matches the same
// "keep asking until valid" style as the rest of the console prompts here.
// ----------------------------------------------------------------------------
bool ask_yes_no(const std::string& prompt) {
    while (true) {
        std::cout << prompt;
        cin >> input;
        if (input == "y" || input == "Y") return true;
        if (input == "n" || input == "N") return false;
        std::cout << "Invalid input, try again." << std::endl;
    }
}

void console_input()
{

    message = L"Do not hold,\n waiting for calibration";

    update_status_window(message, RGB(0, 0, 0), RGB(255, 255, 255));


    // Keep asking until we get a valid choice.
    while (true) {
        std::cout << "Select subject number with an int: \n";
        cin >> input;
        try
        {
            flag_subject_number.store(stoi(input));
            break;
        }catch (...) {
            cout << "Invalid input, try again." << endl;
        }
    }
    message = L"Hold the handle";
    update_status_window(message, RGB(255, 255, 255), RGB(0, 0, 0));


    //Main loop
    while (true){

    const bool autoMode = ask_yes_no("Auto mode? (y/n): ");


    // Mode is always asked, in both manual and auto mode.
    while (true) {
        std::cout << "Select mode with an int: \n";
        std::cout << " 1 or 6 = admittance \n 2 = ProMP prediction variable impedance \n "
                     "3 = Kalman prediction variable impedance \n 4 = ProMP prediction "
                     "fixed impedance \n 5 = Kalman prediction fixed impedance ";
        cin >> input;
        try {
            flag_mode.store(stoi(input));
            break;
        } catch (...) {
            cout << "Invalid input, try again." << endl;
        }
    }

    // Manual mode: movement is asked once per run, as before.
    if (!autoMode) {
        while (true) {
            std::cout << "Select movement with an int: ";
            cin >> input;
            try {
                flag_movement.store(stoi(input));
                break;
            } catch (...) {
                cout << "Invalid input, try again." << endl;
            }
        }
    }

    if (autoMode) {
        // 25 runs total: movement increments every 5 runs (block 0 -> movement
        // 1, block 1 -> movement 2, ...), goal cycles 1-5 within each block
        // per numbers[]. Both movement and goal are derived from the run
        // index, so neither is asked interactively in auto mode.
        const int numbers[25] = {2, 1, 3, 4, 5,
                                  3, 1, 2, 4, 5,
                                  2, 4, 5, 1, 3,
                                  4, 5, 2, 3, 1,
                                  1, 3, 4, 2, 5};
        for (int i = 0; i < 25; ++i) {
            const int movement = 10+10*(i / 5) + numbers[i];
            std::cout << "Doing movement: "<< movement<< "\n";
            flag_movement.store(movement);
            flag_goal.store(numbers[i]);
            run_single_movement();
        }
    } else {
        while (true) {
            std::cout << "Select goal with an int: ";
            cin >> input;
            try {
                flag_goal.store(stoi(input));
                break;
            } catch (...) {
                cout << "Invalid input, try again." << endl;
            }
        }
        run_single_movement();
    }
    /*
    // Keep asking until we get a valid choice.
    while (true) {
        std::cout << "Select mode with an int: \n";
        std::cout << " 1 = admittance \n 2 = ProMP prediction variable impedance \n 3 = Kalman prediction variable impedance \n 4 = ProMP prediction fixed impedance \n 5 = Kalman prediction fixed impedance ";
        cin >> input;
        try
        {
            flag_mode.store(stoi(input));
            break;
        }catch (...) {
            cout << "Invalid input, try again." << endl;
        }
    }
        message = L"Move to the initial position";
        update_status_window(message, RGB(255, 255, 255), RGB(0, 0, 0));

    while (true) {
        std::cout << "Select movement with an int: ";
        cin >> input;
        try
        {
            flag_movement.store(stoi(input));
            break;
        }catch (...) {
            cout << "Invalid input, try again." << endl;
        }
    }

    while (true) {
        std::cout << "Select goal with an int: ";
        cin >> input;
        try
        {
            flag_goal.store(stoi(input));
            break;
        }catch (...) {
            cout << "Invalid input, try again." << endl;
        }
    }
        flag_initial_position.store(0);

        while (flag_initial_position.load() == 0) {
            //Waiting, displaying GO

            std::this_thread::sleep_for(std::chrono::milliseconds(5));
        }

    //Reset the system to start new movement
    flag_reset.store(1);
        message = L"Prepare to go to target "+ std::to_wstring(flag_goal.load());
        update_status_window(message, RGB(0, 0, 200), RGB(255, 255, 255));

        // Displaying GOAL=TARGET
    std::cout << "Press any button to start: \n";
    cin >> input;
        flag_running.store(1);

        message = L"START MOVING TO "+ std::to_wstring(flag_goal.load());
        update_status_window(message, RGB(0, 200, 0), RGB(255, 255, 255));
        while (flag_running.load() == 1) {
            //Waiting, displaying GO

            std::this_thread::sleep_for(std::chrono::milliseconds(5));
        }

        message = L"Target reached, reset to initial position";
        update_status_window(message, RGB(173, 216, 230), RGB(255, 255, 255));

    */
    }
}



// Maximum allowed waiting time during actions
constexpr auto TIMEOUT_PROMISE_DURATION = chrono::seconds{20};

// Create an event listener that will set the promise action event to the exit value
// Will set promise to either END or ABORT
// Use finish_promise.get_future.get() to wait and get the value
function<void(k_api::Base::ActionNotification)>
create_event_listener_by_promise(promise<k_api::Base::ActionEvent>& finish_promise)
{
    return [&finish_promise] (k_api::Base::ActionNotification notification)
    {
        const auto action_event = notification.action_event();
        switch(action_event)
        {
            case k_api::Base::ActionEvent::ACTION_END:
            case k_api::Base::ActionEvent::ACTION_ABORT:
                finish_promise.set_value(action_event);
                break;
            default:
                break;
        }
    };
}

//-----------------------------------------------------------
// Function of impedance control
//-----------------------------------------------------------
// 3 inputs:
// base; base_cyclic; actuator_config = paramters for robotic communication
// p_d, dp_d, ddp_d = desired position, velocity and acceleration
// K_d_diag = diagonal values for desired stiffness matrix

// 1 output:
// status of controller: True or False
//------------------------------------------------------------
bool impedance_control(k_api::Base::BaseClient* base, k_api::BaseCyclic::BaseCyclicClient* base_cyclic,
                       k_api::ActuatorConfig::ActuatorConfigClient* actuator_config,
                       VectorXd& p_d, VectorXd& dp_d, VectorXd& ddp_d, VectorXd& K_d_diag,
                       std::string model_dir, promp_rt::ConditioningMode mode, double override_duration,
                       double record_dt_s, int max_obs, std::vector<double> target_pos) {
    std::cout << "Starting impedance control" << "\n";
    bool return_status = true;

    // Clearing faults
    try {
        base->ClearFaults();
    }
    catch (...) {
        cout << "Unable to clear robot faults" << endl;
        return false;
    }

    k_api::BaseCyclic::Feedback base_feedback;
    k_api::BaseCyclic::Command base_command;

    vector<float> commands;

    auto servoing_mode = k_api::Base::ServoingModeInformation();

    int timer_count = 0;

    // KINOVA feedback (joint space variables)
    VectorXd q(7), dq(7), last_dq(7), ddq(7), torque(7), current(7);

    // Task space variables
    VectorXd u(7), p(6), dp(6), ddp(6);
    MatrixXd T_B7(4,4);  // Rotation matrix

    //EXTERNAL FORCE (TASK SPACE)
    VectorXd f_ext(6);


    // Define the D_n and K_n for nullspace
    VectorXd K_n_diag(7);
    K_n_diag << 0, 5, 40, 0, 0, 0, 20;

    // Time for one control iterative
    double dt = 1.0 / CONTROL_FREQUENCY;
    const double iteration_time = (1.0 / CONTROL_FREQUENCY) * 1000;

    // Buffer to save the last velocity
    last_dq << 0,0,0,0,0,0,0;

    cout << "Initializing the arm for torque control ^^!" << endl;
    try
    {
        // Set the base in low-level servoing mode
        servoing_mode.set_servoing_mode(k_api::Base::ServoingMode::LOW_LEVEL_SERVOING);
        base->SetServoingMode(servoing_mode);
        base_feedback = base_cyclic->RefreshFeedback();

        // Initialize each actuator to their current position
        for (int i = 0; i < ACTUATOR_COUNT; i++)
        {
            commands.push_back(base_feedback.actuators(i).position());

            // Save the current actuator position, to avoid a following error
            base_command.add_actuators()->set_position(base_feedback.actuators(i).position());
        }

        // Send a first frame
        base_feedback = base_cyclic->Refresh(base_command);

        // Set actuatorS in torque mode now that the command is equal to measure
        auto control_mode_message = k_api::ActuatorConfig::ControlModeInformation();
        //control_mode_message.set_control_mode(k_api::ActuatorConfig::ControlMode::POSITION);
        control_mode_message.set_control_mode(k_api::ActuatorConfig::ControlMode::TORQUE);
        for (int id = 1; id < ACTUATOR_COUNT+1; id++)
        {
            actuator_config->SetControlMode(control_mode_message, id);
        }
        auto information_servo_mode = base->GetServoingMode();
        std::cout << "Servoing mode = "
          << information_servo_mode.servoing_mode()
          << std::endl;

        this_thread::sleep_for(chrono::milliseconds(40));

        // Clock to record the time period between two control/measure loop
        auto start_measure = chrono::high_resolution_clock::now();
        auto start_control = chrono::high_resolution_clock::now();

        // Control loop
        cout << "Starting control loop" << endl;

        std::string filename = "C:/Imperial/Final project/predictive_controller/subject" + std::to_string(flag_subject_number.load()) +
                "mode"+std::to_string(flag_mode.load())+"mov"+std::to_string(flag_movement.load())+".csv";

        double initial_position_timer = 0.0;
        double goal_reached_timer = 0.0;

        while (1) // This is an example code
        {

            // KINOVA Feedback: Obtaining gen3 ACTUAL joint positions, velocities, torques & current
            for (int i = 0; i < ACTUATOR_COUNT; i++)
            {
                q[i] = (M_PI/180)*base_feedback.actuators(i).position();
                dq[i] = (M_PI/180)*base_feedback.actuators(i).velocity();
                torque[i] = base_feedback.actuators(i).torque();
                current[i] = base_feedback.actuators(i).current_motor();
            }
            //cout << "Feedback obtained" << endl;
            // GETTING END EFFECTOR EXTERNAL FORCE
            /*
            base_feedback = base_cyclic->RefreshFeedback();
            f_ext << base_feedback.base().tool_external_wrench_force_x(),
                 base_feedback.base().tool_external_wrench_force_y(),
                 base_feedback.base().tool_external_wrench_force_z(),
                 base_feedback.base().tool_external_wrench_torque_x(),
                 base_feedback.base().tool_external_wrench_torque_y(),
                 base_feedback.base().tool_external_wrench_torque_z();
            
            std::cout << "[INFO] External force " << f_ext << "\n";
            */

            // Apply the forward kinematics
            tie(p, T_B7) = forward(q);

            // initilize the controller and filter
            if(timer_count == 0){
                Vector<double, 3> pos = p.head(3);
                ini_controller(pos, T_B7, model_dir, max_obs);
                ini_butterworth();
            }

            // Compute the joint accelerations
            auto end_measure = chrono::high_resolution_clock::now();
            chrono::duration<double> measure_dur = end_measure - start_measure;
            dt = measure_dur.count();
            if (timer_count == 0){
                dt = 1.0 / CONTROL_FREQUENCY;
            }
            ddq = (dq - last_dq) / dt;
            start_measure = end_measure;
            last_dq = dq;
            //cout << "Computed acceleration" << endl;


            VectorXd acc_factor(6);
            /*
            // Impedance controller
            tie(u, dp, ddp, acc_factor) = impedance_controller(q, dq, ddq, T_B7, p_d, dp_d,
                                                               ddp_d, K_d_diag, K_n_diag,CONTROL_FREQUENCY, dt);
            */
            //cout << "Producing admittance matrix" << endl;
            VectorXd K_d_diag_admittance(6), D_d_diag(6), I_d_diag(6);
            D_d_diag << 20.0, 20.0, 20.0, 11.0,  11.0,  11.0;
            //cout << "D_d_diag produced" << endl;
            I_d_diag << 2,  2,  2,  3,  3,  3;
            //cout << "I_d_diag produced" << endl;
            //I_d_diag << 0.5,  0.5,  0.5,  0.3,  0.3,  0.3;
            K_d_diag_admittance << 0.0,  0.0,  0.0,  40.0,  40.0,  40.0;
            //cout << "K_d_diag_admittance produced" << endl;
            /*
            tie(u, dp, ddp, acc_factor) = admittance_controller(q, dq, ddq, T_B7, K_d_diag_admittance, D_d_diag, I_d_diag,
                                                                torque, CONTROL_FREQUENCY, dt);
            */
            //cout << "Starting predictive impedance control" << endl;
            //std::cout << "Starting predictive impedance control" << "\n";
            double sensor = read_sensor();
            //std::cout << "Sensor reading: " << sensor << "\n";

            int goal_reached = 0;
            int reset_happened = 0;
            int initial_position_reached = 0;

            if (flag_reset.load() == 1){
                filename = "C:/Imperial/Final project/predictive_controller/subject" + std::to_string(flag_subject_number.load()) +
                    "mode"+std::to_string(flag_mode.load())+"mov"+std::to_string(flag_movement.load())+".csv";
                //flag_reset.store(0);
            }
            tie(u, dp, ddp, acc_factor, goal_reached, reset_happened, initial_position_reached) = predictive_impedance_controller(q,dq,ddq,T_B7,p_d, dp_d, ddp_d,
                K_d_diag,K_n_diag,CONTROL_FREQUENCY,dt,override_duration,mode,
                record_dt_s,target_pos,K_d_diag_admittance,D_d_diag,I_d_diag,torque, current, sensor,
                flag_mode.load(), flag_reset.load(), flag_goal.load(), flag_running.load(), filename);

            if (reset_happened == 1){
                flag_reset.store(0);
            }

            if (goal_reached == 1)
            {
                if (goal_reached_timer < 0.01)
                {
                    std::wstring message = L"Hold goal position";
                    update_status_window(message, RGB(255, 255, 255), RGB(0, 0, 0));
                }
                goal_reached_timer = goal_reached_timer + record_dt_s;

                if (goal_reached_timer > 0.8)
                {
                    flag_running.store(0);
                }

            }else
            {
                if ((goal_reached_timer > 0.01) && (flag_running.load()==1))
                {
                    std::wstring message = L"Go to target position";
                    update_status_window(message, RGB(0, 255, 0), RGB(255, 255, 255));
                }
                goal_reached_timer = 0.0;
            }

            if (initial_position_reached == 1)
            {
                if ((initial_position_timer < 0.01) && (flag_initial_position.load() == 0))
                {
                    std::wstring message = L"Hold initial position";
                    update_status_window(message, RGB(255, 255, 255), RGB(0, 0, 0));
                }
                initial_position_timer = initial_position_reached + record_dt_s;
                if (initial_position_timer > 0.5)
                {
                    flag_initial_position.store(1);
                    initial_position_timer = 0.0;
                }
            }else
            {
                if (initial_position_timer > 0.01)
                {
                    std::wstring message = L"Go to initial position";
                    update_status_window(message, RGB(200, 0, 0), RGB(255, 255, 255));
                }
                initial_position_timer = 0.0;
            }
            //flag_initial_position.store(1);

            /*
            tie(u, dp, ddp, acc_factor) = position_admittance_controller(q, dq, ddq, T_B7, K_d_diag_admittance, D_d_diag, I_d_diag,
                                                                torque, CONTROL_FREQUENCY, dt);
            */

            // Set the control frequency
            auto end_control = chrono::high_resolution_clock::now();
            auto run_time = chrono::duration<double, milli>(end_control - start_control).count();
            int diff = 0;
            if (run_time < iteration_time) { // if real control loop time < desired iteration time
                auto start_delay = chrono::high_resolution_clock::now();
                diff = (iteration_time - run_time) * 1000;
                chrono::microseconds delay(diff);
                while (chrono::high_resolution_clock::now() - start_delay < delay); // delay the time difference
            }
            auto last_control_time = start_control;
            start_control = chrono::high_resolution_clock::now();
            auto control_duration = chrono::duration<double, milli>(start_control - last_control_time).count();

            for (int i = 0; i < ACTUATOR_COUNT; i++)
            {

                // -- Position
                base_command.mutable_actuators(i)->set_position(base_feedback.actuators(i).position());
                // -- Torque
                base_command.mutable_actuators(i)->set_torque_joint(u[i]);
                /*
                // -- Position
                base_command.mutable_actuators(i)->set_position(u[i]*180/M_PI);
                // -- Torque
                //base_command.mutable_actuators(i)->set_torque_joint(base_feedback.actuators(i).torque());

                std::cout
                    << "Joint " << i
                    << " feedback=" << base_feedback.actuators(i).position()
                    << " u(rad)=" << u[i]
                    << " u(deg)=" << u[i] * 180.0 / M_PI
                    << std::endl;
                */
            }


            // Incrementing identifier ensures actuators can reject out of time frames
            base_command.set_frame_id(base_command.frame_id() + 1);
            if (base_command.frame_id() > 65535)
                base_command.set_frame_id(0);

            for (int idx = 0; idx < ACTUATOR_COUNT; idx++)
            {
                base_command.mutable_actuators(idx)->set_command_id(base_command.frame_id());
            }
            /*
            auto information_servo_mode = base->GetServoingMode();
            std::cout << "Servoing mode = "
                << information_servo_mode.servoing_mode()
                << std::endl;
                */
            try
            {


                base_feedback = base_cyclic->Refresh(base_command, 0);
            }
            catch (k_api::KDetailedException& ex)
            {
                cout << "Kortex exception: " << ex.what() << endl;

                cout << "Error sub-code: " << k_api::SubErrorCodes_Name(k_api::SubErrorCodes((ex.getErrorInfo().getError().error_sub_code()))) << endl;
            }
            catch (runtime_error& ex2)
            {
                cout << "runtime error: " << ex2.what() << endl;
            }
            catch(...)
            {
                cout << "Unknown error." << endl;
            }

            timer_count++;
        }

        cout << "Torque control ^^ completed" << endl;

        // Set actuators back in position
        control_mode_message.set_control_mode(k_api::ActuatorConfig::ControlMode::POSITION);
        for (int id = 1; id < ACTUATOR_COUNT+1; id++)
        {
            actuator_config->SetControlMode(control_mode_message, id);
        }

        cout << "Torque control ^^ clean exit" << endl;

    }
    catch (k_api::KDetailedException& ex)
    {
        cout << "API error: " << ex.what() << endl;
        return_status = false;
    }
    catch (runtime_error& ex2)
    {
        cout << "Error: " << ex2.what() << endl;
        return_status = false;
    }

    // Set the servoing mode back to Single Level
    servoing_mode.set_servoing_mode(k_api::Base::ServoingMode::SINGLE_LEVEL_SERVOING);
    base->SetServoingMode(servoing_mode);

    // Wait for a bit
    this_thread::sleep_for(chrono::milliseconds(2000));

    return return_status;
}


//------------------------------------------
// Function of high-level movement
//-----------------------------------------
// 2 inputs:
// base: communication variables
// q_d: desired joint angular position
//-----------------------------------------
void Move_high_level(k_api::Base::BaseClient* base, VectorXd q_d)
{

    auto action = k_api::Base::Action();

    auto reach_joint_angles = action.mutable_reach_joint_angles();
    auto joint_angles = reach_joint_angles->mutable_joint_angles();

    auto actuator_count = base->GetActuatorCount();

    // Arm straight up
    for (size_t i = 0; i < actuator_count.count(); ++i)
    {
        auto joint_angle = joint_angles->add_joint_angles();
        joint_angle->set_joint_identifier(i);
        joint_angle->set_value(q_d[i]);
    }

    promise<k_api::Base::ActionEvent> finish_promise;
    auto finish_future = finish_promise.get_future();
    auto promise_notification_handle = base->OnNotificationActionTopic(
            create_event_listener_by_promise(finish_promise),
            k_api::Common::NotificationOptions()
    );

    base->ExecuteAction(action);

    const auto status = finish_future.wait_for(TIMEOUT_PROMISE_DURATION);
    base->Unsubscribe(promise_notification_handle);

    if(status != future_status::ready)
    {
        cout << "Timeout on action notification wait" << endl;
    }
    const auto promise_event = finish_future.get();
}

/**
 * @brief Return the value of a named CLI flag, or default_val if absent.
 */
static std::string get_arg(const std::vector<std::string>& args,
                            const std::string& flag,
                            const std::string& default_val = "")
{
    for (size_t i = 0; i + 1 < args.size(); ++i)
        if (args[i] == flag) return args[i + 1];
    return default_val;
}

/**
 * @brief Return true if @p flag appears anywhere in @p args.
 */
static bool has_flag(const std::vector<std::string>& args,
                     const std::string& flag)
{
    return std::find(args.begin(), args.end(), flag) != args.end();
}

// ============================================================
// Training mode
// ============================================================
static int run_train(const std::vector<std::string>& args)
{
    if (args.size() < 4 || has_flag(args, "--help")) {
        std::cerr <<
            "Usage: promp_rt train <demo_folder> <model_dir> [options]\n"
            "\n"
            "Demo CSV: time_s, pos0_rad, pos1_rad, ..., pos(n_dof-1)_rad\n"
            "\n"
            "Options:\n"
            "  --dof N            Number of joints (default 3)\n"
            "  --basis N          RBF basis functions (default 10)\n"
            "  --std S            Basis std-dev; <=0 = auto (default -1)\n"
            "  --steps N          Resampling grid size (default 100)\n"
            "  --max-duration D   Clip demos at D seconds (default: all)\n"
            "  --sg-window W      SG window length, odd (default 9)\n"
            "  --sg-poly P        SG polynomial order (default 4)\n"
            "  --obs-pos-noise V  Position observation noise var (default 1e-4)\n"
            "  --obs-vel-noise V  Velocity observation noise var (default 1e-2)\n"
            "  --goal-pos-noise V Goal position noise var        (default 1e-6)\n"
            "  --goal-vel-noise V Goal velocity noise var        (default 1e-2)\n";
        return 1;
    }

    promp_rt::TrainingConfig cfg;
    cfg.n_dof               = std::stoi(get_arg(args, "--dof",            "3"));
    cfg.num_basis_functions = std::stoi(get_arg(args, "--basis",          "15"));
    cfg.std_bf              = std::stod(get_arg(args, "--std",            "-1"));
    cfg.n_resample_steps    = std::stoi(get_arg(args, "--steps",          "300"));
    cfg.max_demo_duration   = std::stod(get_arg(args, "--max-duration",   "-1"));
    cfg.sg_window           = std::stoi(get_arg(args, "--sg-window",      "9"));
    cfg.sg_poly_order       = std::stoi(get_arg(args, "--sg-poly",        "4"));
    cfg.obs_pos_noise_var   = std::stod(get_arg(args, "--obs-pos-noise",  "0.0001"));
    cfg.obs_vel_noise_var   = std::stod(get_arg(args, "--obs-vel-noise",  "0.01"));
    cfg.goal_pos_noise_var  = std::stod(get_arg(args, "--goal-pos-noise", "0.000001"));
    cfg.goal_vel_noise_var  = std::stod(get_arg(args, "--goal-vel-noise", "0.01"));

    return promp_rt::train_from_folder(args[2], args[3], cfg) ? 0 : 1;
}

// ============================================================
// Prediction mode
// ============================================================
static void print_result(const promp_rt::PredictionResult& r)
{
    if (!r.valid) {
        std::cout << "PRED 0 " << r.current_phase << "\n";
        std::cout.flush();
        return;
    }

    std::cout << "PRED " << r.future_times_s.size()
              << " "     << r.current_phase << "\n";
    /*
    const int n     = static_cast<int>(r.future_times_s.size());
    const int ndims = r.n_dims;
    for (int i = 0; i < n; ++i) {
        std::cout << r.future_times_s[static_cast<size_t>(i)];
        for (int dim = 0; dim < ndims; ++dim) {
            // Interleave mean and std: mean[dim] std[dim]
            std::cout << " " << r.mean_traj[static_cast<size_t>(i * ndims + dim)];
            //          << " " << r.std_traj [static_cast<size_t>(i * ndims + dim)];
        }
        std::cout << "\n";
    }
    */
    //std::cout.flush();
}


static int run_predict(const std::vector<std::string>& args)
{
    if (args.size() < 3 || has_flag(args, "--help")) {
        std::cerr <<
            "Usage: promp_rt predict <model_dir> [options]\n"
            "\n"
            "Options:\n"
            "  --mode past_traj|past_traj_target   (default past_traj)\n"
            "  --target v0 v1 … v(n_dof-1)         target joint positions (rad)\n"
            "  --duration D                         override expected duration (s)\n"
            "  --record-dt DT                       prediction step in s; <=0 = auto from timestamps\n"
            "  --max-obs N                          obs per tick cap (default 0=all)\n"
            "  --input PATH                         read observations from file (default: stdin)\n"
            "\n"
            "stdin:  <time_s> <pos0> <pos1> ...  (one per line)\n"
            "stdout: PRED <n> <phase>  then n lines of predictions\n";
        return 1;
    }
    // Disable C++/C stdio sync for speed
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);


    const std::string model_dir    = args[2];
    const std::string mode_str     = get_arg(args, "--mode", "past_traj");
    const double override_duration = std::stod(get_arg(args, "--duration", "-1"));
    const double record_dt_s       = std::stod(get_arg(args, "--record-dt", "-1"));
    const int    max_obs           = std::stoi(get_arg(args, "--max-obs", "0"));

    // Empty default → read observations from stdin (see loop below).  A
    // non-empty default (e.g. "0") would make the code try to open a file by
    // that name and the documented stdin mode would be unreachable.
    const std::string input_path   = get_arg(args, "--input", "");

    // Construct predictor first so we can query n_dof for target parsing.
    promp_rt::OnlinePredictor predictor(model_dir, max_obs);
    const int n_dof = predictor.get_n_dof();

    // Parse --target  v0  v1  …  v(n_dof-1)
    // Must have exactly n_dof values following the flag.
    std::vector<double> target_pos(static_cast<size_t>(n_dof), 0.0);
    {
        auto it = std::find(args.begin(), args.end(), "--target");
        if (it != args.end()) {
            for (int d = 0; d < n_dof; ++d) {
                ++it;
                if (it == args.end()) {
                    std::cerr << "[predict] --target requires " << n_dof
                              << " values (n_dof=" << n_dof << ")\n";
                    return 1;
                }
                target_pos[static_cast<size_t>(d)] = std::stod(*it);
            }
        }
    }

    const promp_rt::ConditioningMode mode =
        (mode_str == "past_traj_target")
            ? promp_rt::ConditioningMode::PAST_TRAJ_TARGET
            : promp_rt::ConditioningMode::PAST_TRAJ;

    predictor.reset(override_duration);

    // ── Real-time streaming loop ──────────────────────────────────────────


    //while (std::getline(std::cin, line)) { // for manually input line
    std::ifstream input_file;
    std::istream* in = &std::cin;

    if (!input_path.empty()) {
        input_file.open(input_path);
        if (!input_file.is_open()) {
            std::cerr << "[predict] Could not open input file: " << input_path << "\n";
            return 1;
        }
        in = &input_file;
    }

    // Each stdin line produces exactly one PRED block on stdout.
    std::vector<double> pos_buf(static_cast<size_t>(n_dof));
    std::string line;

    while (std::getline(*in, line)) {

        if (line.empty()) break;

        std::istringstream ss(line);
        double time_s{};
        if (!(ss >> time_s)) {
            std::cerr << "[predict] Cannot parse line: " << line << "\n";
            continue;
        }
        bool ok = true;
        for (int d = 0; d < n_dof; ++d) {
            if (!(ss >> pos_buf[static_cast<size_t>(d)])) { ok = false; break; }
        }
        if (!ok) {
            std::cerr << "[predict] Expected " << n_dof
                      << " joint values on line: " << line << "\n";
            continue;
        }

        auto t_start = std::chrono::high_resolution_clock::now();
        promp_rt::PredictionResult result =
            predictor.update_and_predict(
                time_s, pos_buf, mode, record_dt_s, target_pos);
        auto t_end = std::chrono::high_resolution_clock::now();
        double us = std::chrono::duration<double, std::micro>(t_end - t_start).count();
        std::cout << result.n_obs_used << "," << us << "\n";
        //print_result(result);
    }

    return 0;
}



//----------------------------------------------------
// Main function of impedance control
//----------------------------------------------------

int main(int argc, char **argv)
{
    std::vector<std::string> args(argv, argv + argc);


    const std::string model_dir    = args[2];
    const std::string mode_str     = get_arg(args, "--mode", "past_traj");
    const double override_duration = std::stod(get_arg(args, "--duration", "-1"));
    const double record_dt_s       = std::stod(get_arg(args, "--record-dt", "-1"));
    const int    max_obs           = std::stoi(get_arg(args, "--max-obs", "0"));

    std::vector<double> target_pos(static_cast<size_t>(3), 0.0);
    {
        auto it = std::find(args.begin(), args.end(), "--target");
        if (it != args.end()) {
            for (int d = 0; d < 3; ++d) {
                ++it;
                if (it == args.end()) {
                    std::cerr << "[predict] --target requires " << 3
                              << " values (n_dof=" << 3 << ")\n";
                    return 1;
                }
                target_pos[static_cast<size_t>(d)] = std::stod(*it);
            }
        }
    }
    const promp_rt::ConditioningMode mode =
        (mode_str == "past_traj_target")
            ? promp_rt::ConditioningMode::PAST_TRAJ_TARGET
            : promp_rt::ConditioningMode::PAST_TRAJ;

    VectorXd p_d(6), dp_d(6), ddp_d(6), K_d_diag(6);
    // Define the desired position, velocity, acceleration and stiffness here.
     p_d << 0.442,-0.06,0.575,1.4,0,1.4;
    dp_d << 0,0,0,0,0,0;
     ddp_d << 0,0,0,0,0,0;
     K_d_diag << 500,500,500,10,10,10;

    // Experiment loop for different sets
    this_thread::sleep_for(chrono::milliseconds(1000));

    // Create API objects
    auto error_callback = [](k_api::KError err) { cout << "_________ callback error _________" << err.toString(); };

    cout << "Creating transport objects" << endl;
    auto transport = new k_api::TransportClientTcp();
    auto router = new k_api::RouterClient(transport, error_callback);
    transport->connect(IP_ADDRESS, PORT);

    cout << "Creating transport real time objects" << endl;
    auto transport_real_time = new k_api::TransportClientUdp();
    auto router_real_time = new k_api::RouterClient(transport_real_time, error_callback);
    transport_real_time->connect(IP_ADDRESS, PORT_REAL_TIME);

    // Set session data connection information
    auto create_session_info = k_api::Session::CreateSessionInfo();
    create_session_info.set_username("admin");
    create_session_info.set_password("admin");
    create_session_info.set_session_inactivity_timeout(60000);   // (milliseconds)
    create_session_info.set_connection_inactivity_timeout(2000); // (milliseconds)

    // Session manager service wrapper
    cout << "Creating sessions for communication" << endl;
    auto session_manager = new k_api::SessionManager(router);
    session_manager->CreateSession(create_session_info);
    auto session_manager_real_time = new k_api::SessionManager(router_real_time);
    session_manager_real_time->CreateSession(create_session_info);
    cout << "Sessions created" << endl;

    // Create services
    auto base = new k_api::Base::BaseClient(router);
    auto base_cyclic = new k_api::BaseCyclic::BaseCyclicClient(router_real_time);
    auto actuator_config = new k_api::ActuatorConfig::ActuatorConfigClient(router);
    cout << "Services created" << endl;

    // Move to the initial configuration
    VectorXd q_d(7);
    q_d << 0, 0.261887, 3.14162, 4.01422, 2.72361e-05, 0.959906, 1.57;
    q_d = q_d * (180 / M_PI);
    cout << "Moving to initial configuration" << endl;
    //Move_high_level(base, q_d);
    cout << "Moved to initial configuration" << endl;

    start_status_window();  // once, at the top of main()

    std::thread input_thread(console_input);
    input_thread.detach();


    // Atomic flag to control the loop termination for drawScene
    atomic<bool> is_impedance_control_running(true);

    std::thread sensor_thread(sensor_thread_func);

        // Thread for impedance_control function
        thread impedance_thread([&]() {
            auto isOk = impedance_control(base, base_cyclic, actuator_config, p_d, dp_d, ddp_d, K_d_diag,
             model_dir,  mode,  override_duration, record_dt_s, max_obs,target_pos);
            if (!isOk) {
                cout << "There has been an unexpected error in impedance_control() function." << endl;
            }
            // Signal that the impedance control has finished
            is_impedance_control_running = false;
        });

        // Join the threads back to the main thread
        sensor_thread.join();
        impedance_thread.join();

        // Close API session
        session_manager->CloseSession();
        session_manager_real_time->CloseSession();

        // Deactivate the router and cleanly disconnect from the transport object
        router->SetActivationStatus(false);
        transport->disconnect();
        router_real_time->SetActivationStatus(false);
        transport_real_time->disconnect();

        // Destroy the API
        delete base;
        delete base_cyclic;
        delete actuator_config;
        delete session_manager;
        delete session_manager_real_time;
        delete router;
        delete router_real_time;
        delete transport;
        delete transport_real_time;
}




// ============================================================
// Entry point
// ============================================================
/*
int main(int argc, char* argv[])
{
    std::vector<std::string> args(argv, argv + argc);

    if (args.size() < 2 || has_flag(args, "--help") || has_flag(args, "-h")) {
        std::cout <<
            "ProMP_1 real-time variable-DoF training and prediction tool\n"
            "Usage:\n"
            "  promp_rt train   <demo_folder> <model_dir> [--dof N] [options]\n"
            "  promp_rt predict <model_dir> --mode <mode> [options]\n"
            "Run 'promp_rt <subcommand> --help' for full option list.\n";
        return 0;
    }

    try {
        if (args[1] == "train")   return run_train(args);
        if (args[1] == "predict") return run_predict(args);

        std::cerr << "Unknown subcommand '" << args[1]
                  << "'.  Use 'train' or 'predict'.\n";
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "[error] " << e.what() << "\n";
        return 1;
    }
}
*/