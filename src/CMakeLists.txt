# Load rapid_csv Library
add_library(load_data STATIC load_data.cpp)

# Link agains sreg
add_executable(sreg main.cpp)
target_include_directories(sreg PRIVATE /)
