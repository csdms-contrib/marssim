cmake_minimum_required (VERSION 2.6)

enable_language (Fortran)
set (MODEL marssim)
set (BUILD_SHARED_LIBS ON)

set (CMAKE_Fortran_FLAGS -ffree-line-length-none)

set (marssim_srcs
     program_global_variables.f90
     alternate_read_and_write.f90
     lake_flow_routing.f90
     boundary_conditions.f90
     lava_flows.f90
     default_flow_routing.f90
     determine_stream_network.f90
     erosion_iterations.f90
     sediment_routing.f90
     fluvial_slope_erosion.f90
     stream_network_properties.f90
     gradient_and_flow_directions.f90
     surface_erosion_deposition.f90
     groundwater_flow.f90
     weathering.f90
     impact_cratering.f90
     write_debug_information.f90
     initialize_variables.f90
     channel_properties.F90
     eolian_erosion_deposition.F90
     summary_statistics.F90
     mass_wasting.F90
     marssim_irf.f90
     main.f90)

add_library (${MODEL} ${marssim_srcs})
add_executable (run_${MODEL} ${marssim_srcs})
target_link_libraries (marssim m)

install (TARGETS run_${MODEL} RUNTIME DESTINATION bin)
install (TARGETS ${MODEL} DESTINATION lib)
install (FILES
	${PROJECT_BINARY_DIR}/src/accretion_globals.mod
        ${PROJECT_BINARY_DIR}/src/area_globals.mod
        ${PROJECT_BINARY_DIR}/src/crater_globals.mod
        ${PROJECT_BINARY_DIR}/src/eolian_globals.mod
        ${PROJECT_BINARY_DIR}/src/erode_globals.mod
        ${PROJECT_BINARY_DIR}/src/lake_globals.mod
        ${PROJECT_BINARY_DIR}/src/lava_globals.mod
        ${PROJECT_BINARY_DIR}/src/marssim_irf.mod
        ${PROJECT_BINARY_DIR}/src/ocean_globals.mod
        ${PROJECT_BINARY_DIR}/src/seddebug_globals.mod
        ${PROJECT_BINARY_DIR}/src/sedroute_globals.mod
        DESTINATION include)

