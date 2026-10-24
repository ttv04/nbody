add_library(system OBJECT)

target_sources(system
    PRIVATE
        system.cpp
)

target_link_libraries(system
    PRIVATE
        include
)
