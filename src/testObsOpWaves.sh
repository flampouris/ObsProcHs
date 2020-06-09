#!/bin/bash
#========================================================================
#
#          FILE:   testObsOpWaves.sh
# 
#       PURPOSE:   Driver script for ObsOpWaves.exe and tester for 
#                  correct compilation
# 
#         USAGE:   1. Update the Content of the section "Modify Me"
#                  2. Update the access permission to the file, e.g.
#                  >  chmod 755 testObsOpWaves.sh
#                  3. ./testObsOpWaves.sh
# 
#   DESCRIPTION:
# 
#       OPTIONS: All the options are internal in the present file and 
#                have to be modified by the user
#      LOCATION: The file is included to the ObsOpWaves package at the 
#                path [...]/ObsOpWaves/ObsOpExamples/
#  REQUIREMENTS: The script requires the data required for ObsOpWaves.exe
#                WW3 output in grib format, altimeter data in prepbufr 
#                format and the actual executable ObsOpWaves.exe
#        OUTPUT: Check the documentation at the obsophs.f90
#                For this example:
#                 3 outputfiles are produced:
#                    a binary output $ALTFILE.dat 
#                    a text output (activated byt the -textout .true.) $ALTFILE.txt and 
#                    a stats file ErrorStats_HTSGW.txt
#                 The binary and text $ALTFILE output: 
#                 alt_ID              : Sensor ID
#                 YYYY MM DD HH mm ss : time of data acquisition 
#                 lat lon             : lat lon of observation 
#                 altim_hs            : observed hs  
#                 mod_hs              : model hs
#                 altim_ws            : observed wind speed
#                 mod_ws              : model ws 
#                 indRealX,Y          : decimal index in X- and Y- direction for the obs at the model grid  
#                 observation_typical_error  : error variance of observation. 
#                 Format of ErrorStats_HTSGW.txt: 
#                 GuessFileName, FHR, Window Length(TMWND), nobs (number of observations), Bias, RMSE, SI, Percentile 90, P95, P99
#   LIMITATIONS: This script works only at WCOSS or THEIA in conjuction 
#                with Stelios' wgrib2api paths
#          BUGS: ---
#        AUTHOR: Jessica Meixner, jessica.meixner@noaa.gov
#       COMPANY: 
#       VERSION: 1.2
#       CREATED: 2017-Apr-19 20:30 zulu
#      REVISION: Stelios Flampouris, stylianos.flampouris@noaa.gov
#                Documentation, minor changes, Inclusion to ObsOpWaves 
#                package
#
#========================================================================
#
#========================== Modify Me Starts ============================
# 1.svn host path, in case that you want to checkout the latest version
svn_address='https://svnemc.ncep.noaa.gov/projects/ww3_utils/ObsOpWaves'
# 2. Check out the latest version, define: yes or no
bool_co='no'
# 3. Define Server: THEIA or WCOSS
server_name='THEIA'
# 4. Path of the example (assumption this script is located at /ObsOpWaves/
example_path='ObsOpExamples/'
# 5.
# GRIBFILE, wave model output
GRIBFILE=multi_1.glo_30m.t00z.f012.grib2
# FHR, f time of prediction
FHR=12
# TMWND, the time window that the algorithm is looking for data
#  -TMWND<=z+f<=TMWND
TMWND=1.5
#SATFILE, the satelite files with the obs, tip: if user prefers can CAT them
SATFILE='xx115 xx122 xx123 xx124'
#SATFILE='xx122'
#ALTFILE, the name of the output file 
ALTFILE=outputaltfile
# 
#=========================== Modify Me Ends =============================
#
# 1.
#if [ "$bool_co" = "yes" ]
#   then
#   svn co ${svn_address}
#fi
#
# 2. 
module purge
if [ "$server_name" = "WCOSS" ]
   then
   module load ibmpe
   module load w3nco
   module load ics/16.0.3
   module load bufr/v11.0.2
fi
if [ "$server_name" = "THEIA" ]
   then
   module load intel/16.1.150
   module load nco/4.4.5 #it is compiled even without this...
fi
#
# 3.
build_script='build.obsop_hs.sh'
rm ${build_script}
ln -s build.obsop_hs_${server_name}.sh ${build_script}
chmod 700 ${build_script}
./${build_script}
# 4. 
#cd ${example_path}
#5.
#cp -rf ../ObsOpWaves.exe ./

#chmod 700 ObsOpWaves.exe
#./ObsOpWaves.exe -gues $GRIBFILE -davar HTSGW WIND -fhours $FHR -timewnd $TMWND -obssatin $SATFILE -textout .true. -despike .true. 2.3 -obsout $ALTFILE -movavg .true. 7 
