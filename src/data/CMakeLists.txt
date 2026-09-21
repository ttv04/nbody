add_library(data OBJECT)
target_sources(data 
    PRIVATE 
        data.cpp
)

target_link_libraries(data
    PRIVATE 
        include 
        nlohmann_json
)