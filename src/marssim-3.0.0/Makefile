.SUFFIXES: .f90

PROG=marssim

FC=/usr/local/gcc/bin/gfortran

FFLAGS=-ffree-line-length-none

SOURCES=  \
         program_global_variables.f90 \
         alternate_read_and_write.f90 \
         lake_flow_routing.f90 \
         boundary_conditions.f90 \
         lava_flows.f90 \
         default_flow_routing.f90 \
         determine_stream_network.f90 \
         erosion_iterations.f90 \
         sediment_routing.f90 \
         fluvial_slope_erosion.f90 \
         stream_network_properties.f90 \
         gradient_and_flow_directions.f90 \
         surface_erosion_deposition.f90 \
         groundwater_flow.f90 \
         weathering.f90 \
         impact_cratering.f90 \
         write_debug_information.f90 \
         initialize_variables.f90 \
         channel_properties.f90 \
         eolian_erosion_deposition.f90 \
         summary_statistics.f90 \
         mass_wasting.f90

#         read_and_write_data_files.f90

MAIN_SOURCE= \
         main_program.f90

OBJS=${SOURCES:.f90=.o}

print:
	echo ${OBJS}

${PROG}: ${OBJS}
	${FC} -o ${PROG} ${MAIN_SOURCE} ${OBJS}

clean:
	rm -f core ${PROG} ${OBJS} *mod

.f90.o:
	${FC} -c ${FFLAGS} $<

