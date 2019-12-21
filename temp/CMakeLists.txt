find_package(ROOT 6.16 CONFIG REQUIRED)
find_package(TBB REQUIRED)

add_executable(temp1 main.cpp)
target_link_libraries(temp1 SpaceCharge ROOT::Core ROOT::Gpad TBB::tbb)
