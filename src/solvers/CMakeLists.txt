add_library(solvers OBJECT)

target_sources(solvers
    PRIVATE 
        naive.cpp
        barneshut.cpp
        gravity.cpp
)

target_link_libraries(solvers
    PRIVATE 
        include
)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(solvers PUBLIC OpenMP::OpenMP_CXX)
endif()