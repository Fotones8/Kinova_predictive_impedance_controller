cmake_minimum_required(VERSION 3.5)

project(ExperimentalCodes  LANGUAGES C CXX)

set(CMAKE_INCLUDE_CURRENT_DIR   ON)
set(CMAKE_VERBOSE_MAKEFILE      ON)
set(CMAKE_COLOR_MAKEFILE        ON)

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


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

#set_source_files_properties(
#        thirdParty/torque_sensor/Thread/mutex.c
#        thirdParty/torque_sensor/Thread/mutex.h
#        thirdParty/torque_sensor/Thread/rq_thread.c
#        thirdParty/torque_sensor/Thread/rq_thread.h
#        thirdParty/torque_sensor/rq_int.h
#        thirdParty/torque_sensor/rq_sensor_com.c
#        thirdParty/torque_sensor/rq_sensor_com.h
#        thirdParty/torque_sensor/rq_sensor_state.c
#        thirdParty/torque_sensor/rq_sensor_state.h
#        PROPERTIES LANGUAGE C
#)

option(USE_CONAN "Use the Conan package manager to automatically fetch the Kortex API" ON)
option(DOWNLOAD_API "Automatically download the API if conan is not used" ON)

# Activate C++ 11
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_C_STANDARD 11)  # Or C11 or whatever C version you are using

macro(configure_msvc_runtime)

    # Default to statically-linked runtime.
    if("${MSVC_RUNTIME}" STREQUAL "")
        set(MSVC_RUNTIME "static")
    endif()

    # Set compiler options.
    set(variables
            CMAKE_C_FLAGS_DEBUG
            CMAKE_C_FLAGS_MINSIZEREL
            CMAKE_C_FLAGS_RELEASE
            CMAKE_C_FLAGS_RELWITHDEBINFO
            CMAKE_CXX_FLAGS_DEBUG
            CMAKE_CXX_FLAGS_MINSIZEREL
            CMAKE_CXX_FLAGS_RELEASE
            CMAKE_CXX_FLAGS_RELWITHDEBINFO
    )

    if(${MSVC_RUNTIME} STREQUAL "static")
        message(STATUS
                "MSVC -> forcing use of statically-linked runtime."
        )
        foreach(variable ${variables})
            if(${variable} MATCHES "/MD")
                string(REGEX REPLACE "/MD" "/MT" ${variable} "${${variable}}")
            endif()
        endforeach()

    else()
        message(STATUS
                "MSVC -> forcing use of dynamically-linked runtime."
        )
        foreach(variable ${variables})
            if(${variable} MATCHES "/MT")
                string(REGEX REPLACE "/MT" "/MD" ${variable} "${${variable}}")
            endif()
        endforeach()

    endif()

endmacro()


include_directories(${PROJECT_SOURCE_DIR}/thirdParty/cxxopts/)

if(MSVC)
    configure_msvc_runtime()
else()
    add_compile_options(-Wall)
    add_compile_options(-Wno-reorder)

endif()

if(UNIX)
    add_definitions(-D_OS_UNIX)
elseif(WIN32)
    add_definitions(-D_OS_WINDOWS -DNOMINMAX)
    if(MSVC)
        add_compile_options(/bigobj)
    endif()
endif()

if(USE_CONAN)

    include(${PROJECT_SOURCE_DIR}/thirdParty/conan.cmake)

    conan_check(REQUIRED)
    conan_add_remote(
            NAME kinova_public
            URL https://artifactory.kinovaapps.com/artifactory/api/conan/conan-public)

    if(UNIX)
        conan_cmake_run(REQUIRES kortex_api_cpp/2.6.0-r.3@kortex/stable
                SETTINGS kortex_api_cpp:compiler=gcc
                SETTINGS kortex_api_cpp:compiler.version=5
                SETTINGS compiler.libcxx=libstdc++11
                PROFILE_AUTO build_type
                BASIC_SETUP
                UPDATE)
    elseif(WIN32)
        if(MSVC)

            _get_msvc_ide_version(_VISUAL_VERSION)

            # TODO Conan artifacts should not require to specify a target to download
            if (_VISUAL_VERSION EQUAL 14)
                set(kortex_api_cpp_target "msvc-2015")
            elseif(_VISUAL_VERSION EQUAL 15)
                set(kortex_api_cpp_target "msvc-2017")
            elseif(_VISUAL_VERSION EQUAL 16)
                set(kortex_api_cpp_target "msvc-2019")
            endif()

            conan_cmake_run(REQUIRES kortex_api_cpp/2.6.0-r.3@kortex/stable
                    PROFILE_AUTO build_type
                    BASIC_SETUP
                    UPDATE)

        else()
            conan_cmake_run(REQUIRES kortex_api_cpp/2.6.0-r.3@kortex/stable
                    SETTINGS kortex_api_cpp:compiler=gcc
                    SETTINGS kortex_api_cpp:compiler.version=5
                    SETTINGS compiler.libcxx=libstdc++11
                    PROFILE_AUTO build_type
                    BASIC_SETUP
                    UPDATE)
        endif()
    endif()

    link_libraries(${CONAN_LIBS})

else() # Not using Conan

    # Setup Kortex Api Path
    if(NOT KORTEX_SUB_DIR)
        set(KORTEX_SUB_DIR "")
    else()
        set(KORTEX_SUB_DIR "${KORTEX_SUB_DIR}/")
    endif()

    set(KORTEX_DIR "${PROJECT_SOURCE_DIR}/kortex_api/${KORTEX_SUB_DIR}")

    if(CMAKE_BUILD_TYPE EQUAL "Debug")
        set(KORTEX_LIB_SUBDIR "debug")
    else()
        set(KORTEX_LIB_SUBDIR "release")
    endif()

    # Download the API
    if(DOWNLOAD_API)
        if(UNIX)
            execute_process(COMMAND ./download_kortex_api.sh ${KORTEX_SUB_DIR}
                    WORKING_DIRECTORY ../scripts
                    RESULT_VARIABLE DOWNLOAD_API_RESULT
                    OUTPUT_VARIABLE DOWNLOAD_API_OUTPUT)
            if(NOT DOWNLOAD_API_RESULT EQUAL 0)
                message("Kortex API was not downloaded prior to running CMake.")
                message(FATAL_ERROR ${DOWNLOAD_API_OUTPUT})
            endif()
        elseif(WIN32)
            execute_process(COMMAND ./download_kortex_api.bat ${KORTEX_SUB_DIR}
                    WORKING_DIRECTORY ../scripts
                    RESULT_VARIABLE DOWNLOAD_API_RESULT
                    OUTPUT_VARIABLE DOWNLOAD_API_OUTPUT)
            if(NOT DOWNLOAD_API_RESULT EQUAL 0)
                message("Kortex API was not downloaded prior to running CMake.")
                message(FATAL_ERROR ${DOWNLOAD_API_OUTPUT})
            endif()
        endif()
    endif()
    if(UNIX)
        link_libraries(${KORTEX_DIR}lib/${KORTEX_LIB_SUBDIR}/libKortexApiCpp.a)
    elseif(WIN32)
        link_libraries(${KORTEX_DIR}lib/${KORTEX_LIB_SUBDIR}/KortexApiCpp.lib)
    endif()

    # Add Include Directories
    include_directories(${KORTEX_DIR}include)
    include_directories(${KORTEX_DIR}include/client)
    include_directories(${KORTEX_DIR}include/common)
    include_directories(${KORTEX_DIR}include/messages)
    include_directories(${KORTEX_DIR}include/client_stubs)

endif()

# link other libs
if(UNIX)
    link_libraries(pthread)
elseif(WIN32)
    link_libraries(winMM ws2_32)
else()
    MESSAGE(FATAL_ERROR "Unknown os! Not supported yet")
endif()

set(EIGEN_DIR "${PROJECT_SOURCE_DIR}/libs/eigen")
include_directories(${EIGEN_DIR})

# set(EasyX_ROOT "C:/Program Files (x86)/Microsoft Visual Studio/2019/Community/VC/Auxiliary/VS")
# set(EasyX_INC  ${EasyX_ROOT}/include)
# set(EasyX_LINK  ${EasyX_ROOT}/lib/x64)

# include_directories(
#        ${EasyX_INC}
# )

# link_directories(
#        ${EasyX_LINK}
# )

# include_directories(${PROJECT_SOURCE_DIR}/thirdParty/torque_sensor)
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



# Create executable for each example
add_executable(Experimental_Codes main.cpp
        Jacobian.cpp
        Jacobian.h
        Fwd_kinematics.cpp
        Fwd_kinematics.h
        Dynamics.cpp
        Dynamics.h
#        Interface.cpp
#        Interface.h
        Controller.cpp
        Controller.h
#        Trajectory.cpp
#        Trajectory.h
#        Saver.cpp
#        Saver.h
#        noise.cpp
#        noise.h
#        thirdParty/torque_sensor/Thread/mutex.c
#        thirdParty/torque_sensor/Thread/mutex.h
#        thirdParty/torque_sensor/Thread/rq_thread.c
#        thirdParty/torque_sensor/Thread/rq_thread.h
#        thirdParty/torque_sensor/rq_int.h
#        thirdParty/torque_sensor/rq_sensor_com.c
#        thirdParty/torque_sensor/rq_sensor_com.h
#        thirdParty/torque_sensor/rq_sensor_state.c
#        thirdParty/torque_sensor/rq_sensor_state.h
        Filter.cpp
        Filter.h
        src/training.cpp
        src/prediction.cpp

)

target_include_directories(Experimental_Codes PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_libraries(Experimental_Codes PRIVATE
        promp          # ProMP_1 shared library (brings Eigen3::Eigen transitively)
)


target_link_libraries(Experimental_Codes
#        EasyXa
#        EasyXw
)

target_compile_features(Experimental_Codes PRIVATE cxx_std_17)

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
