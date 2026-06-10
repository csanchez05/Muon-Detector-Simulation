cmake_minimum_required(VERSION 3.16...3.27)
project(muon_sim)

# Find Geant4 package, activating all available UI and Vis drivers by default
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

# Locate all source and header files
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Create the executable
add_executable(muon_sim muon_sim.cc ${sources} ${headers})
target_include_directories(muon_sim PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(muon_sim ${Geant4_LIBRARIES})

# Copy macros to the build directory (we will write init.mac later)
file(COPY init.mac DESTINATION ${PROJECT_BINARY_DIR} 
     PATTERN "*.mac")