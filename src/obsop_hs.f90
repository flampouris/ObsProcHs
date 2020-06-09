!>@file obsop_hs.f90
!!@brief The \c obsop_hs.f90 includes the main program of ObsOpWaves and the approapriate documentation.
!!@mainpage ObsOpWaves
!!
!!@section Preface_doc Preface
!!
!!@subsection Purpose_of_doc Purpose 
!!This document/webpage is the extended manual of ObsOpWaves program. It provides guidance and template 
!!material which is intended to assist the users and the potential developers to use the capabilities of
!!the software at their full extend.
!!
!!@subsection Use_of_doc Use of the Manual
!!Trivial reminder but some times useful: <b>Just read it and keep it for further reference.</b>
!!
!!@section Introduction_man Introduction.
!!The ObsOpWaves is the forward operator for wave data assimilation systems. The ObsOpWaves 
!!prepares the input files for any Wave Data Assimilation System developed for the NOAA/NWS/EMC 
!!by Stelios Flampouris. In addition, the intermediate output of the program is the collocated 
!!data of observations and predictions.\n
!!This document describes:
!! - The basics about ObsOpWaves
!! - Where the ObsOpWaves can be retireved
!! - How to be installed and 
!! - How to be used.
!!
!!This is the user manual and system documentation of version 1.1. 
!!
!!@section scope Scope of ObsOpWaves
!!The objective is the creation of an operational modular interface combining the wave
!!measurements and any wave model, which follow the standards of the WMO: 
!!
!!  - BUFR  : Binary Universal Form for the Representation of meteorological data and\n
!!  - GRIB2 : Representation and Exchange of Regularly Spaced Data In Binary Form for the model output\n
!!
!!The outputs of the program are:
!!  - The input for the LETKF-Waves, 
!!  - The collocated data and 
!!  - The basic error statistics of the model.
!!
!!The code is inspired by the example for sst for LETKF by S. Penny and T. Miyoshi,
!!but the program has been redesign for 1. focusing on the waves and 2. optimization. \n
!!
!!The source code ObsOpWaves has been written exclusively in FORTRAN 95 according to 
!!the NOAA/NCO programing standards and it is machine independent. \n
!!
!!@subsection contact_info General Information
!!
!!@author Dr. Rer. Nat. Stylianos Flampouris \n
!!IMSG at EMC/NCEP/NOAA\n
!!NOAA Centers for Weather and Climate Prediction\n
!!5830 University Research Court\n
!!College Park, Maryland 20740, USA\n
!!phone: +1-301-683-3789 // +1-857-472-0097\n
!!email: stylianos.flampouris@noaa.gov \n
!!
!!@version 1.02
!!The current version of the program is focused on the significant wave
!!height.  It inputs observations (yo) and a single forecast (xf) and 
!!computes the innovations (yo-H(xf)) associated with that member, after
!!applying the appropriate operators.
!! 
!!The ObsOpWaves is work in progress and it develops as the Wave Data Assimilation Systems develop. \n
!!
!!In case of interest, anybody is welcomed to use, add, modify the code according to 
!!their needs and desires. But they are kindly asked to follow the existing code as example 
!!(a good read: <http://www.fortran90.org/src/best-practices.html>),
!!document all the modifications and additions with doxygen, add their names and affiliations to the authors lists. \n
!!@date October 2016 The Creation
!!@date 06-Nov-2016 Multiple files with satellite observations can be used
!!@date 16-Nov-2016 The wind obs are also imported and exported at the required
!!format for validation 
!!@date 17-Nov-2016 Jason3 (xx124) data are imported.
!!@date 20-Nov-2016 The collocated data are sorted
!!@date 18-Nov-2016 Debugging option has been extended, see below
!!@date 20-Jan-2017 The collocation was optimized for cartesian grids.
!!@date 24-Jan-2017 New input argument is added: "textout", when .true. the collocated data are exported, at file named "obsout".txt (see obsout argument)
!!@date 24-Jan-2017 Error Statistics module is added and tables (text files) named: ErrorStats_DAQuantity.txt are created.
!!@date 12-Feb-2017 The format of ErrorStats_DAQuantity.txt has been updated to include number of observations.
!!@date 02-Nov-2017 Despiking was added.
!!@date 05-Nov-2017 Moving Average smoothing was added.
!!
!!@copyright 
!!The ObsOpWaves is free software: you can redistribute it and/or modify
!!it under the terms of the GNU General Public License as published by
!!the Free Software Foundation, either version 3 of the License, or
!!(at your option) any later version./n
!!This program is distributed in the hope that it will be useful,
!!but WITHOUT ANY WARRANTY; without even the implied warranty of
!!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!GNU General Public License for more details./n
!!You should have received a copy of the GNU General Public License
!!along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!and all the related material are  other GNU./n
!!
!!Libraries that are used and they have not been developed by the author(s) of 
!!ObsOpWaves, are used according to their licenses./n
!!
!!@section source_code_access Access to the Source Code
!! The code can be retrieved from the NOAA internal svn
!!> https://svnemc.ncep.noaa.gov/trac/ww3_utils
!! OR
!!> email: stylianos.flampouris@noaa.gov
!!
!!@section Fast Fast Start
!!The code comes with an installation script (<c>testObsOpWaves.sh</c>) for WCOSS and THEIA and an example for checking the installation.
!!
!!Follow the instructions in the <c>testObsOpWaves.sh</c>, modify it and it should work.
!!
!!@section Compilation How To Compile ObsOpWaves
!!The source code can be compiled with gfortan and intel fortran. Three external libraries are required: 
!!1. The grib2api by Wesley Ebisuzaki (wesley.ebisuzaki@noaa.gov). The library is currently available by request only. Users of THEIA can use Stelios'
!! copy of grib2api, the exact locations are given below. If you prefer your local copy
!!of the API, just copy it! Important Note: the Makefiles have been modified to work at THEIA. /n
!!2. The bufr library, already installed at NOAA's computational systems or available at <http://www.nco.ncep.noaa.gov/sib/decoders/BUFRLIB/> /n
!!3. The production utilities by NOAA/NCO <http://www.nco.ncep.noaa.gov/pmb/docs/libs/w3lib/ncep_w3lib.shtml/> and <http://www.nco.ncep.noaa.gov/pmb/codes/nwprod/>. Both websites have useful tools and cleverly implemented code. \n
!!The extra: The doxygen <http://www.doxygen.nl/> is not used for the compilation but for the documentation of the code./n
!!
!!@subsection Fast_Mode Fast Mode
!!> ifort -c -no-wrap-margin kinds.f90 constants.f90 variables_obsop_hs.f90 general_functions.f90 grib2_ww3_io.f90 qsort_mod.f90 despike_mod.f90 read_obs.f90 error_statistics_mod.f90 obsop_hs.f90 -I /scratch4/NCEPDEV/ocean/save/Stylianos.Flampouris/grib2_working/li
!!
!!@subsection Debugging_Mode Debugging Mode 
!!> ifort -c -no-wrap-margin -D_REAL8_ -DWRF -qopenmp -O3 -traceback -check all -fp-stack-check -g -fp-model source -convert big_endian -assume byterecl kinds.f90 constants.f90 variables_obsop_hs.f90 general_functions.f90 grib2_ww3_io.f90 qsort_mod.f90 despike_mod.f90 read_obs.f90 error_statistics_mod.f90 obsop_hs.f90 -I /scratch4/NCEPDEV/ocean/save/Stylianos.Flampouris/grib2_working/lib
!!
!!@subsection Create_Executable Executable! 
!!> ifort -o ObsOpWaves.exe -mkl -Wl,-Map,loadmap.txt *.o -L/scratch3/NCEPDEV/nwprod/lib/bufr/v10.2.5/ -lbufr_v10.2.5_d_64 -L/scratch3/NCEPDEV/nwprod/lib/w3nco/v2.0.6/ -lw3nco_4 -L/scratch4/NCEPDEV/ocean/save/Stylianos.Flampouris/grib2_working/lib -lwgrib2 -lgfortran -lz -lm
!!
!!@subsection MPI MPI
!!The code is NOT currently paralalized, but it is designed to be parallelized in the space of the observations. /n
!! Still the program can be compiled as following:
!!> mpiifort -c -no-wrap-margin  -D_REAL8_ -DWRF -openmp -O3 -traceback -g -fp-model source -convert big_endian -assume byterecl kinds.f90 constants.f90 variables_obsop_hs.f90 grib2_ww3_io.f90 read_obs.f90 general_functions.f90 obsop_hs.f90 -I/scratch4/NCEPDEV/ocean/save/Stylianos.Flampouris/ensemble/Ocean-LETKF/src/model_specific/ww3/grib2/lib/
!!
!!> mpiifort -o ObsOpWaves.exe -mkl -Wl,-Map,loadmap.txt *.o -L/scratch3/NCEPDEV/nceplibs/gwv/nwprod.0927.2013/lib.0927.2013/ -lbufr_v10.2.2_8_64 -L/scratch3/NCEPDEV/nwprod/lib/g2/v2.5.0/libg2_v2.5.0_4.a -L/scratch3/NCEPDEV/nwprod/lib/w3nco/v2.0.6/ -lw3nco_4 -L/scratch4/NCEPDEV/ocean/save/Stylianos.Flampouris/ensemble/Ocean-LETKF/src/model_specific/ww3/grib2/lib -lwgrib2 -lgfortran -lz -lm -openmp
!!
!!@section Execution_of_ObsOpWaves How To Run ObsOpWaves
!!
!!The ObsOpWaves.exe is user friendly, for terminal friends (<-this is a joke), it runs with specific flags which determine the values of the necessary variables.
!!
!!@subsection ToRun To Run it
!!The following line is a typical call of the ObsOpWaves 
!!> ObsOpWaves.exe -gues [guess_filename].grib2 -fhours [time of prediction in hours] -davar [Variable Name] -obssatin [filename of satellite_observations].bufr -obsinsin [filename of in situ obervations].bufr -superob [.true.] [degrees in x and y] -timewnd [time in hours for synchronization]  
!!
!!@subsection Description_flags Description of the flags 
!!All the flags are case sensitive Details about the flags and the default values are given below: \n
!!\b -gues 	[guess_filename].grib2 					: [guess_filename].grib2 is the name of the grib2 file with the model output. Default value \c gues.grb2 \n
!!\b -davar	[Variable Name_1,..., [Variable Name_n]	: Variable Name is the grib2 name for the quantity of interest, the user can define as they exist in the file. Default value \c HTSGW \n
!!\b -fhours 	[hours]									: Time of prediction in hours. It is an integer define as [dT X No of time steps]. Default value \0, which corresponds to the analysis time \n
!!\b -obssatin N x [filename of satellite_observations].bufr	: N Filenames of satellite observations in (prep)BUFR format. Default value \c obssatin.prepbufr \n
!!\b -obsinsin [filename of in situ obervations].bufr 	: Filename of in situ observations in (prep)BUFR format. Default value \c obsinsin.prepbufr \n
!!\b -obsout [filename for ouput]                        : Filename of for saving in binary format the collated data. Default value \c obsout.dat \n
!!\b -textout [flag] :  Logical flag for having formatted output of the collocated data, the file has the same name as the "obsout" with suffix ".txt". Default value \c .false. \n 
!!\b -superob [.flag.] [degrees in x and y]	: Logical flag to create or not Super Observations. Default value \c .false. \n
!!In case of .true. the distance in degrees in X-direction and optionally in Y-direction can be given, in order to define the area of each Super Observation. Default value \c 0.5 \n
!!\b -timewnd [hours]	: Real value for defining the time window for assimilation. Default value \c 1.0 \n
!!\b -filtersat  [.flag.] : Logical flag to filter or not the satellite data.  Default value \c .true. \n
!!\b -despike  [.flag.] [despikeStdFac]: Logical flag to despike spatially the data according to miu-despikeStdFac*sigma <= value <=miu+despikeStdFac*sigma. \n
!!Default values  <c> .true.</c> and <c>2.0</c>\n
!!\b -movavg  [.flag.] [Length of MA, No of obs]: Logical flag to apply moving average for period equal to the defined No of obs. \n
!!Default values <c> .true.</c> and <c>5</c> \n
!!\b -debug [flag]	: Logical flag for debugging, mainly screen outputs of each step; limited use. Default value <c>.false.</c> \n
!!\n There are more flags but currently not in use. \n
!!
!!@subsection Output Output
!!The collocated (in temporal and physical space) information is stored at the structure \b collocated. \n
!!The error statistics output has the following format:
!!GuessFileName, F number, Window Length, nobs, Bias, RMSE, SI, Percentile 90, P95, P99 
 
!!
!!@subsection Debugging Debugging
!!The option debugging has been extended. The modules are given access to the
!!variable dodebug1 from variables_obsop_hs module, therefore it is used for
!!visualizing intermediate steps.\n
!! The main debugging outputs are the following: \n
!! 1. Many comments appear on the screen, starting with the name of the
!! subroutine \n
!! 2. When active more files are created:\n
!!    model_out_[grib2 variable].txt -- The gridded model data for the specific variable
!!has been exported.
!! 
!!@section Further_Development Further Development
!!The ObsOpWaves is work in progress and it develops as the Wave Data Assimilation Systems are developped. \n
!!In case of interest, anybody is welcomed to use, add, modify the code according to their needs and desires. But they are kindly asked to follow the existing code as example (a good read: <http://www.fortran90.org/src/best-practices.html>),
!!document all the modifications and additions with doxygen, add their names and affiliations to the authors lists.
!!
!>@bug The current version of the ObsOpWaves collocates the wind from model and
!!observations. But it does NOT match the quantites from the 2 data sources
!!according to the name automatically. So in case that the user wants wind
!!validation, she/he has to put it second at the calling line "-davar WHATEVER WIND" \n
!!A dictionary that translates the \b grib2 keywords to \b bufr keywords has to be
!!implemented. Maybe there is a problem with the dimensions of the arrays for the observations. \n
!<
!>@bug FIXED : The collocated superobs are not sorted. \n
!!
!!@todo The following list is ideas in random order that they should and will be implemented: \n
!!\Support Unstructured grid. \n
!!\> Create Makefile \n
!!\> All the variables to be case insensitive; use the string module of George Benthien. \n
!!\> WORK in Progress : Extensive use of debugging mode. \n
!!\> Implementation of MPI. \n
!!\> Clean up the code and extend the documentation. \n 
!!\> The current version reads in only the significant wave height from the *bufr files, even though that satellites provide K- wave spectra and the in-situ sensors at dcom have also Wave Peak Period and Direction. 
!!\> Define all the attributes for exported structure according to the needs of users. 
!<
!!*/
program obsop_hs
   use variables_obsop_hs
   use constants_hs, only: i1,i2,i25
   use grib2_ww3_io, only : setup_read_member, read_member, define_domain
   use read_obs_mod, only : read_obs
   use general_func, only : collocation, write_obs, writeMatrix,write_error_stat
   use error_statistics, only : AmIWrong
   implicit none
   character(i25) :: myname =' obsop_hs '
   integer :: ll
!
   call process_command_line
!  
   call setup_read_member(guesfile, len(DAquantity),DAquantity(i1),f_number,nlon,nlat) 
!   
   call allocate_member(nlon,nlat,nv2d) 
!
   call read_member(guesfile,  nv2d, len(DAquantity), DAquantity  &
                   , f_number, nlon, nlat                         &
                   , latg, long, mbr_data, gstime                 )
!
   call define_domain(nlon, nlat, latg, long, size(mbr_domain), mbr_domain)
!
   call read_obs (nsatflnm, obsSatinfile, obsInSinfile, gstime    &
                 , winlen, mbr_domain, do_superobs, do_despiking  &
                 , despikeStdFac, do_ma, MaLngth, SprObsDim, nobs ) 
!
   call allocate_collo(nobs)
!
   call collocation(nlon, nlat, nv2d, latg, long, mbr_data        &
                   , size(all_obs), all_obs, collocated,nobs      )
!  
   call write_obs(collocated,nobs,obsoutfile,dotextout)
!
   call AmIWrong(nobs, collocated%mbr_value(i1)                   &
                , collocated%obs_value, bias, rmse, si, prcntErr  &
                , prcntl                                          )     !howv
!
   call write_error_stat(DAquantity(1), guesfile, f_number, winlen&
                        , nobs, bias, rmse, si, prcntErr          )
!
   call AmIWrong(nobs, collocated%mbr_value(i2)                   &
                , collocated%wspa, bias, rmse, si, prcntErr       &
                , prcntl                                          )     !wspa
!
   call write_error_stat(DAquantity(2),guesfile,f_number,winlen,nobs,bias, rmse, si, prcntErr)
!   
   call destroy_arrays
!
end program obsop_hs
