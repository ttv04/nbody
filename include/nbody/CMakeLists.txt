add_library(include INTERFACE)

target_sources(include
    INTERFACE 
        FILE_SET HEADERS
        BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/..
        FILES 
            core.hpp
            integrators.hpp
            solvers.hpp
            state.hpp
            generators.hpp
)