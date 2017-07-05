GWGEN_F90: The FORTRAN code base for a global weather generator
===============================================================

Welcome! This repository contains the code base for the **G**lobal **W**eather
**Gen**erator. This files can be used to run the weather generator.

To compile the software, you need a FORTRAN 95 compiler, e.g. via::

    sudo apt-get install gfortran

on debian, or

    brew install gcc

on Mac OSX using homebrew.

To create the executable, just modify the Makefile_ with the correct fortran
compiler and run::

    make

This will create the `weathergen` executable in the same directory.

The weather generator requires a comma separated file with 11 columns, that are

station id (character string of 11 characters)
    a unique identifier of the weather station
lon (float)
    the longitude of the weather station (only necessary if the
    ``use_geohash``  namelist parameter in the main namelist of
    ``weathergen.nml`` is True (the default))
lat (float)
    the latitude of the weather station (see ``lon``)
year (int)
    the year of the month
month (int)
    the month
min. temperature (float)
    the minimum temperature degrees Celsius
max. temperature (float)
    the maximum temperature in degrees Celsius
cloud fraction (float)
    the mean cloud fraction during the month between 0 and 1
wind speed (float)
    the mean wind speed during the month in m/s
precipitation (float)
    the total precipitation in the month in mm/day
wet (int)
    the number of wet days in the month

and can be run via::

    ./weathergen monthly_input.csv daily_output.csv

.. _Makefile: https://github.com/ARVE-Research/gwgen_f90/blob/master/Makefile
