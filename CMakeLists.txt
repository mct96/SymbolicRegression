cmake_minimum_required(VERSION 3.15)

# Initialize project.
project(SymbolicRegression CXX)

# Define and link executable and librares.
add_subdirectory(lib/)
add_subdirectory(src/)

if (TARGET rapidcsv)
   message(INFO ": Library rapidcsv - FOUND")
else()
   message(ERROR ": Library rapidcsv - NOT FOUND")
endif()

if (TARGET argh)
   target_link_libraries(sreg PRIVATE argh)
   message(INFO ": Library argh - FOUND")
else()
   message(ERROR ": Library argh - NOT FOUND")
endif()

if (TARGET load_data)
   target_link_libraries(load_data PRIVATE rapidcsv)   
   target_link_libraries(sreg PRIVATE load_data)
   message(INFO ": Library load_data - FOUND")
else()
    message(ERROR ": Library load_data - NOT FOUND")
endif()

if (TARGET stats)
   target_link_libraries(sreg PRIVATE stats)
   message(INFO ": Library stats - FOUND")
else()
    message(ERROR ": Library stats - NOT FOUND")
endif()

if (TARGET rnd)
   target_link_libraries(sreg PRIVATE rnd)
   message(INFO ": Library rnd - FOUND")
else()
    message(ERROR ": Library rnd - NOT FOUND")
endif()

set_property(TARGET sreg rnd PROPERTY CXX_STANDARD 17)

# Export usage to docs directory.
configure_file(src/usage.txt docs/usage.txt COPYONLY)
