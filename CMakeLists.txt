cmake_minimum_required(VERSION 3.18)
project(spacewind VERSION 0.2.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

option(SPACEWIND_STRICT "Enable strict compiler warnings" ON)
option(SPACEWIND_SANITIZERS "Enable address and undefined behavior sanitizers" OFF)

set(SPACEWIND_SOURCES
    src/math3.c
    src/rng.c
    src/gravity_grid.c
    src/nbody.c
    src/metric.c
    src/fields.c
    src/spacecraft.c
    src/integrator.c
    src/plasma.c
    src/mhd1d.c
    src/simulation.c
)

add_library(spacewind STATIC ${SPACEWIND_SOURCES})
target_include_directories(spacewind PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(spacewind PUBLIC m)

if(SPACEWIND_STRICT)
    target_compile_options(spacewind PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes
    )
endif()

if(SPACEWIND_SANITIZERS)
    target_compile_options(spacewind PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(spacewind PUBLIC -fsanitize=address,undefined)
endif()

add_executable(spacewind_cli src/cli.c)
set_target_properties(spacewind_cli PROPERTIES OUTPUT_NAME spacewind)
target_link_libraries(spacewind_cli PRIVATE spacewind)

add_executable(spacewind_tests tests/test_main.c)
target_link_libraries(spacewind_tests PRIVATE spacewind)

enable_testing()
add_test(NAME spacewind_tests COMMAND spacewind_tests)
