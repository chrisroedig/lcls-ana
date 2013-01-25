# LCLS AMO Analysis Code

This is code for use with the LCLS myana analysis framework. The focus is to provide simple functionality in order to preform fast, preliminary dat areduction and analysis, during experimental runs.

The intention is to drop these files, including the compile script into a directory that contains the myana framework.

## Note:

The code has not been recently tested using actual XTC data. The code builds without errors, but it should be test run on some smaple data.

## How to use

* place these files in the same dir as `myana.cc`
* you will be replacing `comp` (compiler script) so rename the old one first
* make sure comp is executable `chmod +x comp`
* build the code `./comp`
* make sure the directory outdata exists
* it is very helfull to symlink your xtc directories into the working dir `ln -s /path/to/xtc xtcdata`
* build the programs `./comp`
* run the desired routine ex. `./etof_extra -f mydata.xtc`
* use the n flag to test on the first few shots ex. `./etof_extra -f mydata.xtc -n 1000`

## Included analysis routines

* `etof_extra.cc` 
  * averages all each ETOF
  * converts time to energy
  * corrects for ETOF transmission efficiency
  * smooths data with moving average
  * outputs raw, corrected and smoothed data as ASCII data files
   
* `etof_scan.cc` 
  * sorts traces by FEL / IR timing delay
  * helpfull for finding t-zero during experiments
  * produces 2D ASCII matrix for each ETOF
  * currently implemented via laser AMO delay stage and PCAV signals
   
* `etof_ptrac.cc`
  * intended to sort data by streaking energy displacement
  * currently the code averages electron detections in one of the ETOFs and uses this as a sorting variable
  * outputs a 2D ASCII matrix for each ETOF
  * outputs histogram of streaking energies
   


