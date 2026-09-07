add_library(integrators OBJECT)

target_sources(integrators
    PRIVATE 
        euler.cpp
        eulercromer.cpp
        leapfrog.cpp
)

target_link_libraries(integrators
    PRIVATE 
        include
)