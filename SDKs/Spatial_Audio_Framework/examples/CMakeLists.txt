
message(STATUS "Configuring SAF examples:")

if(NOT DEFINED saf_example_list)
    # List of examples which have the common structure (<example>.h/.c, <example>_internal.h/.c):
    set(saf_example_list
        calibration
endif()

# Find number of examples
list(LENGTH saf_example_list num_examples)
math(EXPR num_examples "${num_examples} - 1")

# Loop over all examples
foreach(current_example RANGE ${num_examples})
    # Get current example string
    list(GET saf_example_list ${current_example} current_string)

    # Set-up project
    message(STATUS "  ${current_string}")
    project("${example_prefix}${current_string}" LANGUAGES C)
    add_library(${PROJECT_NAME} STATIC)

    # Link with SAF
    target_link_libraries(${PROJECT_NAME} PRIVATE saf)

    # Source files
    target_sources(${PROJECT_NAME}
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/${current_string}/${current_string}_internal.c
        ${CMAKE_CURRENT_SOURCE_DIR}/src/${current_string}/${current_string}_internal.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/${current_string}/${current_string}.c
    )

    # enable compiler warnings
    if(UNIX)
        target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra)
    endif()

    # Include directory
    target_include_directories(${PROJECT_NAME}
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
endforeach()
