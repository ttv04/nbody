add_library(core OBJECT)

target_sources(core
    PRIVATE 
        simulation.cpp
)

target_link_libraries(core
    PRIVATE 
        include
)