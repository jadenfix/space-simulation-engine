#ifndef SPACEWIND_VERSION_H
#define SPACEWIND_VERSION_H

#define SPACEWIND_VERSION_MAJOR 0
#define SPACEWIND_VERSION_MINOR 2
#define SPACEWIND_VERSION_PATCH 0
#define SPACEWIND_VERSION_STRING "0.2.0"

#if defined(__clang__)
#define SPACEWIND_COMPILER_ID "clang"
#define SPACEWIND_COMPILER_VERSION __clang_version__
#elif defined(__GNUC__)
#define SPACEWIND_COMPILER_ID "gcc"
#define SPACEWIND_COMPILER_VERSION __VERSION__
#elif defined(_MSC_VER)
#define SPACEWIND_COMPILER_ID "msvc"
#define SPACEWIND_COMPILER_VERSION "unknown"
#else
#define SPACEWIND_COMPILER_ID "unknown"
#define SPACEWIND_COMPILER_VERSION "unknown"
#endif

#if defined(__FAST_MATH__)
#define SPACEWIND_FAST_MATH_JSON "true"
#else
#define SPACEWIND_FAST_MATH_JSON "false"
#endif

#endif
