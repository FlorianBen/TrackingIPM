find_package(ROOT 6.16 CONFIG REQUIRED)
find_package(TBB REQUIRED)
find_package(Boost REQUIRED)

add_executable(temp1 main.cpp field_h5.hpp vector.hpp)
target_link_libraries(temp1 SpaceCharge ROOT::Core ROOT::Gpad TBB::tbb h5cpp)
target_compile_features(temp1 PUBLIC cxx_std_11)

add_executable(temp2 main2.cpp field_h5.hpp vector.hpp)
target_link_libraries(temp2 SpaceCharge TBB::tbb h5cpp ROOT::Core ROOT::Gpad ROOT::Hist ROOT::Graf)
target_compile_features(temp2 PUBLIC cxx_std_11)