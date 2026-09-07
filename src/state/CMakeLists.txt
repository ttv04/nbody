add_library(state OBJECT)

target_sources(state
    PRIVATE 
        state.cpp
)

target_link_libraries(state
    PRIVATE include
)