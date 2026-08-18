cmake_minimum_required(VERSION 3.18)
project(spacewind_robustness_assurance VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
option(SPACEWIND_ROBUSTNESS_SANITIZERS "Enable ASan and UBSan" OFF)

add_library(spacewind_robustness_lib STATIC
    ../src/core.c
    ../src/physics.c
    ../src/invariants.c
    ../src/experiment.c
    ../src/evidence.c
    ../src/pde_checks.c
)
target_include_directories(spacewind_robustness_lib PUBLIC ../include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(spacewind_robustness_lib PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes -Werror
    )
endif()
if(SPACEWIND_ROBUSTNESS_SANITIZERS AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(spacewind_robustness_lib PUBLIC
        -fsanitize=address,undefined -fno-omit-frame-pointer
    )
    target_link_options(spacewind_robustness_lib PUBLIC
        -fsanitize=address,undefined -fno-omit-frame-pointer
    )
endif()
if(NOT MSVC)
    target_link_libraries(spacewind_robustness_lib PUBLIC m)
endif()

add_executable(spacewind_robustness_checks ../tests/test_robustness.c)
target_link_libraries(spacewind_robustness_checks PRIVATE spacewind_robustness_lib)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(spacewind_robustness_checks PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes -Werror
    )
endif()

enable_testing()
file(MAKE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/../output)
add_test(NAME spacewind_robustness_checks COMMAND spacewind_robustness_checks)
set_tests_properties(spacewind_robustness_checks PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/..
)
