cmake_minimum_required(VERSION 3.18)
project(spacewind_pde_assurance VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
option(SPACEWIND_PDE_SANITIZERS "Enable ASan and UBSan" OFF)

add_library(spacewind_pde_assurance STATIC
    ../src/core.c
    ../src/pde_checks.c
)
target_include_directories(spacewind_pde_assurance PUBLIC ../include)

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(spacewind_pde_assurance PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes -Werror
    )
endif()
if(SPACEWIND_PDE_SANITIZERS AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(spacewind_pde_assurance PUBLIC
        -fsanitize=address,undefined -fno-omit-frame-pointer
    )
    target_link_options(spacewind_pde_assurance PUBLIC
        -fsanitize=address,undefined -fno-omit-frame-pointer
    )
endif()
if(NOT MSVC)
    target_link_libraries(spacewind_pde_assurance PUBLIC m)
endif()

add_executable(spacewind_pde_checks ../tests/test_pde_checks.c)
target_link_libraries(spacewind_pde_checks PRIVATE spacewind_pde_assurance)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(spacewind_pde_checks PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes -Werror
    )
endif()

enable_testing()
file(MAKE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/../output)
add_test(NAME spacewind_pde_checks COMMAND spacewind_pde_checks)
set_tests_properties(spacewind_pde_checks PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/..
)
