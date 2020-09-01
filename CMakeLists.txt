cmake_minimum_required(VERSION 3.15)

# Initialize project.
project(SymbolicRegression CXX)

# Define and link executable and librares.
add_subdirectory(lib/)
add_subdirectory(src/)

if (TARGET rapidcsv)
   message("rapidcsv FOUND")
else()
   message("rapidcsv NOT FOUND")
endif()

if (TARGET argh)
   message("argh FOUND")
else()
   message("argh NOT FOUND")
endif()

target_link_libraries(load_data PRIVATE rapidcsv)
target_link_libraries(sreg PRIVATE load_data)                    
set_property(TARGET sreg PROPERTY CXX_STANDARD 17)

# Export usage to docs directory.
configure_file(src/usage.txt docs/usage.txt COPYONLY)
