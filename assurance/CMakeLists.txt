cmake_minimum_required(VERSION 3.18)
project(spacewind_assurance VERSION 0.6.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

option(SPACEWIND_ASSURANCE_SANITIZERS "Enable address and undefined-behavior sanitizers" OFF)
option(SPACEWIND_ASSURANCE_WERROR "Treat warnings as errors" ON)
option(SPACEWIND_ASSURANCE_ENGINE_INTEGRATION "Build independent checks against the root engine" ON)

add_library(spacewind_assurance STATIC
    src/core.c
    src/physics.c
    src/invariants.c
    src/experiment.c
    src/evidence.c
    src/field_ledger.c
    src/charge_deposition.c
    src/spectral_poisson.c
    src/maxwell1d_reference.c
    src/covariant_em.c
)
target_include_directories(spacewind_assurance PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

function(spacewind_assurance_strict target)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes
        )
        if(SPACEWIND_ASSURANCE_WERROR)
            target_compile_options(${target} PRIVATE -Werror)
        endif()
    endif()
endfunction()

spacewind_assurance_strict(spacewind_assurance)

if(SPACEWIND_ASSURANCE_SANITIZERS AND
   CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(spacewind_assurance PUBLIC
        -fsanitize=address,undefined -fno-omit-frame-pointer
    )
    target_link_options(spacewind_assurance PUBLIC
        -fsanitize=address,undefined -fno-omit-frame-pointer
    )
endif()

if(NOT MSVC)
    target_link_libraries(spacewind_assurance PUBLIC m)
endif()

function(add_assurance_test target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} PRIVATE spacewind_assurance)
    spacewind_assurance_strict(${target})
    add_test(NAME ${target} COMMAND ${target})
    set_tests_properties(${target} PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endfunction()

enable_testing()
file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/output)

add_assurance_test(spacewind_assurance_core tests/test_assurance.c)
add_assurance_test(spacewind_assurance_invariants tests/test_invariants.c)
add_assurance_test(spacewind_assurance_experiment tests/test_experiment.c)
add_assurance_test(spacewind_assurance_evidence tests/test_evidence.c)
add_assurance_test(spacewind_assurance_field_ledger tests/test_field_ledger.c)
add_assurance_test(spacewind_assurance_charge_deposition tests/test_charge_deposition.c)
add_assurance_test(spacewind_assurance_spectral_poisson tests/test_spectral_poisson.c)
add_assurance_test(spacewind_assurance_maxwell1d_reference tests/test_maxwell1d_reference.c)
add_assurance_test(spacewind_assurance_covariant_em tests/test_covariant_em.c)
add_assurance_test(spacewind_assurance_maxwell_covariant_integration tests/test_maxwell_covariant_integration.c)

set(SPACEWIND_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)
if(SPACEWIND_ASSURANCE_ENGINE_INTEGRATION AND
   EXISTS ${SPACEWIND_ROOT}/src/plasma.c AND
   EXISTS ${SPACEWIND_ROOT}/src/math3.c AND
   EXISTS ${SPACEWIND_ROOT}/include/spacewind/plasma.h)
    function(add_pic_engine_integration_test target source)
        add_executable(${target}
            ${source}
            ${SPACEWIND_ROOT}/src/plasma.c
            ${SPACEWIND_ROOT}/src/math3.c
        )
        target_include_directories(${target} PRIVATE
            ${SPACEWIND_ROOT}/include
        )
        target_link_libraries(${target} PRIVATE spacewind_assurance)
        spacewind_assurance_strict(${target})
        add_test(NAME ${target} COMMAND ${target})
        set_tests_properties(${target} PROPERTIES
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
    endfunction()

    add_pic_engine_integration_test(
        spacewind_assurance_pic_transport_integration
        tests/test_pic_transport_integration.c
    )
    add_pic_engine_integration_test(
        spacewind_assurance_spectral_poisson_integration
        tests/test_spectral_poisson_integration.c
    )
    add_pic_engine_integration_test(
        spacewind_assurance_fdtd_spectral_integration
        tests/test_fdtd_spectral_integration.c
    )
elseif(SPACEWIND_ASSURANCE_ENGINE_INTEGRATION)
    message(STATUS "Root engine sources unavailable; skipping cross-layer PIC/FDTD integration")
endif()
