cmake_minimum_required(VERSION 3.16)
project(cs2_monte C)

set(CMAKE_C_STANDARD 11)

#include_directories(.)

add_executable(cs2_monte
        adjust.c
        energy.c
        lattice.c
        mc_move.c
        mc_main.c
        ran_uniform.c
        ran_uniform.h
        readinput.c
        sample.c
        system.c
        system.h
        writepdb.c
        ewald.c
        ewald.h
        mc_move.h
        energy.h)


target_link_libraries(cs2_monte m)
