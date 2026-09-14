add_library(nbody STATIC)

target_compile_features(nbody PUBLIC cxx_std_20)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(nbody PUBLIC OpenMP::OpenMP_CXX)
endif()

add_subdirectory(state)
add_subdirectory(integrators)
add_subdirectory(solvers)
add_subdirectory(generators)
add_subdirectory(core)

target_link_libraries(nbody
    PUBLIC 
        include
    PRIVATE 
        state
        integrators
        solvers
        generators
        core
)

