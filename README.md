# LCLS AMO Analysis Code

This is code for use with the LCLS myana analysis framework. The focus is to provide simple functionality in order to preform fast, preliminary dat areduction and analysis, during experimental runs.

The intention is to drop these files, including the compile script into a directory that contains the myana framework.

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


