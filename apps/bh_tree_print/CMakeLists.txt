add_executable(bh_tree_print)
target_sources(bh_tree_print
    PRIVATE
        main.cpp
)

target_link_libraries(bh_tree_print
    PRIVATE
        nbody
)
