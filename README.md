# DISPLAY (Detection and Identification of SPectroscopical Lines for Atacama large milimeter-submilimeter arraY) #

The objective of this project is to identify spectroscopical lines from ALMA-like cubes of spectral data.

Authors:
* Andrés Riveros
* Karim Pichara
* Pavlos Protopapas
* Diego Mardones

The syntethic data was generated using the ASYDO (Astronomical SYnthetic Data Observatory) Project:
 * Mauricio Araya
 * Teodoro Hochfärber


## Introduction ##

### The Problem ###

### Our Proposal ###

## Downloading and Installing ##

> git clone https://github.com/ChileanVirtualObservatory/DISPLAY.git

### Dependencies ###
 * scipy
 * numpy
 * matplotlib
 * astropy
 * pysqlite
 * ASYDO (package included)


###  Installing ###
The DISPLAY package and database creation files should be copy to your working directory. For example, in a unix-based system:

>  user@machine:~/yourproject$ cp -r ../path\_to\_DISPLAY/src/display .
>  user@machine:~/yourproject$ cp ../path\_to\_ASYDO/src/db\* .

### ASYDO Data ###

To create synthetic data, instructions to use the DISPLAY Project are available in (https://github.com/ChileanVirtualObservatory/ASYDO)


