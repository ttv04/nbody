add_library(generators)

target_sources(generators
    PRIVATE
        uniform_2d.cpp
        uniform_3d.cpp
        uniform_galaxy_2d.cpp
        uniform_galaxy_3d.cpp
)

target_link_libraries(generators
    PRIVATE 
        include
)