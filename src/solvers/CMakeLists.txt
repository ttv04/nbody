add_library(solvers OBJECT)

target_sources(solvers
    PRIVATE 
        naive.cpp
)

target_link_libraries(solvers
    PRIVATE 
        include
)