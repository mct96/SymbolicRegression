
add_library(load_data STATIC load_data.cpp)
add_library(stats STATIC statistics.cpp)
add_library(rnd STATIC random.cpp)

# Link agains sreg
add_executable(sreg main.cpp symbolic_regression.cpp)
target_include_directories(sreg PRIVATE /)
