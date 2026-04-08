add_executable(dtc.exe dtc.cpp)
target_link_libraries(dtc.exe bench_lib nlohmann_json::nlohmann_json nwgraph_hpx)
add_dependencies(bench dtc.exe)

add_executable(dbfs.exe dbfs.cpp)
target_link_libraries(dbfs.exe bench_lib nlohmann_json::nlohmann_json nwgraph_hpx)
add_dependencies(bench dbfs.exe)

add_executable(dpr.exe dpr.cpp)
target_link_libraries(dpr.exe bench_lib nlohmann_json::nlohmann_json nwgraph_hpx)
add_dependencies(bench dpr.exe)