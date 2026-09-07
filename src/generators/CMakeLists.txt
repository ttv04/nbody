add_library(generators)

target_sources(generators
    PRIVATE
        uniform_galaxy_2d.cpp
)

target_link_libraries(generators
    PRIVATE 
        include
)