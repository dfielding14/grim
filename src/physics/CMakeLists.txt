add_library(physics physics.c physics.h conduction.c viscosity.c)
target_link_libraries(physics geometry problem)
