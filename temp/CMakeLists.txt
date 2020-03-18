find_package(ROOT 6.16 CONFIG REQUIRED)
find_package(TBB REQUIRED)
find_package(h5cpp REQUIRED)

add_executable(temp1 main.cpp field_h5.hpp vector.hpp)
target_link_libraries(temp1 SpaceCharge ROOT::Core ROOT::Gpad TBB::tbb h5cpp)
target_compile_features(temp1 PUBLIC cxx_std_11)
