add_executable(visualizer)
target_sources(visualizer
    PRIVATE
        main.cpp
)

target_link_libraries(visualizer
    PRIVATE
        nbody
        raylib
)
