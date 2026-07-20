cmake_minimum_required(VERSION 3.17)
project(promp_realtime LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---------------------------------------------------------------------------
# Compiler version check – require C++17-capable toolchain
# ---------------------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
        message(FATAL_ERROR
            "GCC ${CMAKE_CXX_COMPILER_VERSION} is too old.\n"
            "Minimum required: GCC 9.0 (for full C++17 support).\n"
            "Upgrade with: sudo apt install g++-11")
    endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "10.0")
        message(FATAL_ERROR
            "Clang ${CMAKE_CXX_COMPILER_VERSION} is too old.\n"
            "Minimum required: Clang 10.0 (for full C++17 support).\n"
            "Upgrade with: sudo apt install clang-12")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "19.20")
        message(FATAL_ERROR
            "MSVC ${CMAKE_CXX_COMPILER_VERSION} is too old.\n"
            "Minimum required: Visual Studio 2019 (MSVC 19.20).")
    endif()
endif()

message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

# ---------------------------------------------------------------------------
# Eigen3 – required external dependency
# ---------------------------------------------------------------------------
# Install Eigen3 before building:
#   Ubuntu / Debian : sudo apt install libeigen3-dev
#   macOS (Homebrew): brew install eigen
#   Windows (vcpkg) : vcpkg install eigen3:x64-windows
#   Windows (conda) : conda install -c conda-forge eigen
#
# To use a custom Eigen3 location pass: -DEIGEN3_ROOT=/path/to/eigen
# ---------------------------------------------------------------------------
set(EIGEN3_ROOT "" CACHE PATH
    "Optional path to a custom Eigen3 installation (parent of the Eigen/ folder).")

if(EIGEN3_ROOT)
    find_package(Eigen3 3.3 QUIET NO_MODULE HINTS "${EIGEN3_ROOT}")
    if(NOT Eigen3_FOUND)
        # Synthesise the target from the provided path so promp_1's
        # find_package(Eigen3) resolves correctly.
        if(NOT EXISTS "${EIGEN3_ROOT}/Eigen/Core")
            message(FATAL_ERROR
                "EIGEN3_ROOT is set to '${EIGEN3_ROOT}' but Eigen/Core was not found there.\n"
                "Provide the path to the directory that contains the 'Eigen/' folder.")
        endif()
        file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/eigen3_custom")
        file(WRITE "${CMAKE_BINARY_DIR}/eigen3_custom/Eigen3Config.cmake"
"if(NOT TARGET Eigen3::Eigen)
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES \"${EIGEN3_ROOT}\")
endif()
set(Eigen3_FOUND TRUE)
set(Eigen3_VERSION \"custom\")
")
        set(Eigen3_DIR "${CMAKE_BINARY_DIR}/eigen3_custom")
        find_package(Eigen3 CONFIG REQUIRED)
        message(STATUS "Eigen3: using custom path ${EIGEN3_ROOT}")
    endif()
else()

    find_package(Eigen3 3.3...6.0 NO_MODULE
        HINTS /usr/include/eigen3 /usr/local/include/eigen3)

    if(NOT Eigen3_FOUND)
        message(FATAL_ERROR
            "Eigen3 (>= 3.3) was not found on this system.\n"
            "\n"
            "Install it with one of the following commands, then re-run cmake:\n"
            "\n"
            "  Ubuntu / Debian : sudo apt install libeigen3-dev\n"
            "  macOS (Homebrew): brew install eigen\n"
            "  Windows (vcpkg) : vcpkg install eigen3:x64-windows\n"
            "  Windows (conda) : conda install -c conda-forge eigen\n"
            "\n"
            "Or point to a custom Eigen3 directory:\n"
            "  cmake .. -DEIGEN3_ROOT=/path/to/eigen-3.x.x")
    endif()
    message(STATUS "Eigen3: found system installation (version ${Eigen3_VERSION})")
endif()

# ---------------------------------------------------------------------------
# ProMP_1 library – bundled in lib/promp_1, no install needed
# ---------------------------------------------------------------------------
set(PROMP1_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/promp_1"
    CACHE PATH "Path to the bundled promp_1 source directory.")

if(NOT EXISTS "${PROMP1_SRC_DIR}/CMakeLists.txt")
    message(FATAL_ERROR
        "Bundled promp_1 source not found at: ${PROMP1_SRC_DIR}\n"
        "Expected: ${CMAKE_CURRENT_SOURCE_DIR}/lib/promp_1/CMakeLists.txt\n"
        "The lib/promp_1 directory must be present inside promp_realtime/.")
endif()

message(STATUS "ProMP_1: using bundled source at ${PROMP1_SRC_DIR}")

# Suppress extras we do not need from the promp_1 build.
set(BUILD_PYTHON_MODULE  OFF CACHE BOOL "" FORCE)
set(BUILD_DOCUMENTATIONS OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING        OFF CACHE BOOL "" FORCE)

# Forward Eigen3_DIR so promp_1's own find_package(Eigen3) resolves correctly.
if(DEFINED Eigen3_DIR)
    set(Eigen3_DIR "${Eigen3_DIR}" CACHE PATH "" FORCE)
endif()

add_subdirectory("${PROMP1_SRC_DIR}" "${CMAKE_BINARY_DIR}/promp_1_build")
# The add_subdirectory above defines the shared library target named 'promp'.

# ---------------------------------------------------------------------------
# promp_realtime executable
# ---------------------------------------------------------------------------
add_executable(promp_rt
    src/training.cpp
    src/prediction.cpp
    src/main.cpp
)

target_include_directories(promp_rt PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_libraries(promp_rt PRIVATE
    promp          # ProMP_1 shared library (brings Eigen3::Eigen transitively)
)

target_compile_features(promp_rt PRIVATE cxx_std_17)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
message(STATUS "")
message(STATUS "=== promp_realtime build summary ===")
message(STATUS "  Compiler    : ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  Build type  : ${CMAKE_BUILD_TYPE}")
message(STATUS "  Eigen3      : ${Eigen3_VERSION}")
message(STATUS "  ProMP_1 src : ${PROMP1_SRC_DIR}")
message(STATUS "  Output      : promp_rt")
message(STATUS "====================================")

# ---------------------------------------------------------------------------
# Install (optional)
# ---------------------------------------------------------------------------
install(TARGETS promp_rt RUNTIME DESTINATION bin)
