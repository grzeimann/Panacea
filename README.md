(NOTE: This is still in progress and will finish 3/01/2019.  Beware of absolute fluxing issues.)


# Panacea v1.0 (Automatic LRS2 Pipeline)
## Table of Contents
[Overview](https://github.com/grzeimann/Panacea/blob/master/README.md#Overview)

[TACC](https://github.com/grzeimann/Panacea/blob/master/README.md#Working-on-TACC)

[Data Access](https://github.com/grzeimann/Panacea/blob/master/README.md#How-to-get-your-data)

[Data Products](https://github.com/grzeimann/Panacea/blob/master/README.md#Data-Products)

[Running Panacea](https://github.com/grzeimann/Panacea/blob/master/README.md#Running-the-reductions-yourself)

[Panacea: the code](https://github.com/grzeimann/Panacea/blob/master/README.md#Code-Description)

[FAQ](https://github.com/grzeimann/Panacea/blob/master/README.md#frequently-asked-questions)


## Overview
This package is the reduction pipeline for LRS2 observations at the Hobby Eberly Telescope. Every morning the pipeline reduces data taken 
the previous night.  Below we discuss the algorithms and products of Panacea, how to access your data reductions, and how to run the 
pipeline yourself with varying options. All of the data reduction products live on the Texas Advanced Computing Center (TACC).  We start 
with the instructions to log on to TACC, and where you reductions are placed.

## Working on TACC 
The reductions are designed to be run on TACC where a copy of the raw data lives.  We will describe how to get started on TACC,  where the automatic reduction products live, how to run the code yourself, and the products that are produced.

### Signing up for an account
https://portal.tacc.utexas.edu/
<p align="center">
  <img src="images/tacc_create_account.png" width="650"/>
</p>

After creating an accounting using the link above, please send Greg Zeimann <gregz@astro.as.utexas.edu> your TACC username and he will add you to the HET group.  When that step is complete, you can ssh into TACC using:
```
ssh -Y USERNAME@wrangler.tacc.utexas.edu
```

## How to get your data
The reduction pipeline run each morning puts your data products in the following path:
```
/work/03946/hetdex/maverick/LRS2/PROGRAM-ID
```
where PROGRAM-ID, is your program number, for example HET19-1-999.  To get all of the current reductions for your program, simply:
```
scp -r username@wrangler.tacc.utexas.edu:/work/03946/hetdex/maverick/LRS2/PROGRAM-ID .
```
You merely have to use your "username" and your "PROGRAM-ID" and you can copy over your products.  Now, the data reduction products are
extensive, that is to say they for every Mb of raw data there is 24 Mb of reduced data.  Without going into the data products yet,
you may just a single product or a single night.  Below is an example, which grabs all spectra within your program for a given data:
```
scp username@wrangler.tacc.utexas.edu:/work/03946/hetdex/maverick/LRS2/PROGRAM-ID/spec*20190105*.fits .
```

## Data Products
There are three main data products: spectrum*.fits, multi*.fits, and *cube*.fits.  The first product, spectrum*.fits, 
is produced for all exposures and all channels.  Within the fits image, lie rows corresponding to different attributes. 
```
row1: wavelength (air)
row2: extracted object spectrum (f_lambda: ergs/s/cm^2/A)
row3: extracted sky spectrum from same aperture and weighting as object (s_lambda: ergs/s/cm^2/A)
row4: error for extracted object spectrum (e_f_lambda: ergs/s/cm^2/A)
row5: error for extracted sky spectrum (e_s_lambda: ergs/s/cm^2/A)
row6: response function (ergs / e-)
```

The multi*{uv,orange,red,farred}.fits are multi-extension fits files and contain the following attributes:

```
Rectified Spectra: flux calibrated spectrum (object + sky) for each fiber
Rectified Sky Model:flux calibrated sky spectrum for each fiber
Rectified Sky Subtracted Spectra: flux calibrated sky subtracted spectrum for each fiber
Rectified Error Frame: flux calibrated error spectrum for each fiber
Collapsed image: a collapsed frame for visualization of the source(s)
Positions (IFU, Focal, Sky): ifu x and y positions, focal x and y position, and ra and dec
Extracted Spectra and Response: This is identical to the spectrum*.fits extension above
ADR: The atmospheric differential refraction as a function of wavelength.  The columns are wavelength, x_adr, y_adr
CCD Wavelength: The wavelength of each pixel in the 2d frame
Image: the initial reduction of the 2d raw frame.
Flat Fielded image: same as the image frame above but divided by the flat field (fiber profile and fiber to fiber normalization)
Central Trace Pixels: location of the pixels for each fiber (central two pixels)
Cosmics: identified cosmics in the central four pixels of the trace
Unrectified Spectra: Unrectified, uncalibrated spectra for each fiber
```

## Running the reductions yourself
This section covers how to run your own reductions with modifications to achieve specific science objectives.

### Setting up your Python environment
To begin on TACC, point to the common python environment. In your home "~/.bashrc" file, add the following line at the bottom:
```
export PATH=/home/00115/gebhardt/anaconda2/bin:/work/03946/hetdex/maverick/bin:$PATH
```

### Running Panacea in the command line
To run in the command line, TACC wants users to create an interactive development environment which basically gets you a single CPU to 
yourself.  Just type the following:
```
idev
```

Then you can check out the Panacea reduction options:
```
python /work/03730/gregz/maverick/Panacea/full_lrs2_reduction.py -h

usage: full_lrs2_reduction.py [-h] [-d DATE] [-s SIDES] [-o OBJECT] [-uf]
                              [-cf] [-cw CENTRAL_WAVE] [-wb WAVELENGTH_BIN]
                              [-sx SOURCE_X] [-sy SOURCE_Y]
                              [-ssd STANDARD_STAR_DATE]
                              [-sso STANDARD_STAR_OBSID]

optional arguments:
  -h, --help            show this help message and exit
  -d DATE, --date DATE  Date for reduction
  -s SIDES, --sides SIDES
                        "uv,orange,red,farred"
  -o OBJECT, --object OBJECT
                        Object name, no input reduces all objects
  -uf, --use_flat       Use FLT instead of Twi
  -cf, --correct_ftf    Correct fiber to fiber
  -cw CENTRAL_WAVE, --central_wave CENTRAL_WAVE
                        Central Wavelength for collapsed Frame
  -wb WAVELENGTH_BIN, --wavelength_bin WAVELENGTH_BIN
                        Wavelength Bin to collapse over (+/- bin size)
  -sx SOURCE_X, --source_x SOURCE_X
                        Source's x position at the central_wave
  -sy SOURCE_Y, --source_y SOURCE_Y
                        Source's y position at the central_wave
  -ssd STANDARD_STAR_DATE, --standard_star_date STANDARD_STAR_DATE
                        Standard Star Date for response function, example:
                        20181101
  -sso STANDARD_STAR_OBSID, --standard_star_obsid STANDARD_STAR_OBSID
                        Standard Star ObsID for response function, example:
                        0000012
```

If you want to reduce a given object on a given night you can use the following options:

```
python /work/03730/gregz/maverick/Panacea/full_lrs2_reduction.py -d DATE -o TARGET_NAME -s "uv"
```

You can reduce any side you want, above I choose the "uv" channel, and the TARGET_NAME only has to be in the full name of the target
(e.g., HD which is in HD_19445_056_E).

### Running Panacea in batch
To run a reduction of a given target on a given date for all four channels simply:
```
cdw
cp /work/03946/hetdex/maverick/run_lrs2/runlrs2general .
runlrs2general DATE TARGET_NAME
```

You will see an immediate output like:
```
----------------------------------------------------------------
          Welcome to the Maverick Supercomputer                 
----------------------------------------------------------------

No reservation for this job
--> Verifying valid submit host (login2)...OK
--> Verifying valid jobname...OK
--> Enforcing max jobs per user...OK
--> Verifying availability of your home dir (/home/03730/gregz)...OK
--> Verifying availability of your work dir (/work/03730/gregz/maverick)...OK
--> Verifying valid ssh keys...OK
--> Verifying access to desired queue (gpu)...OK
--> Verifying job request is within current queue limits...OK
--> Checking available allocation (Hobby-Eberly-Telesco)...OK
Submitted batch job 900134
```

This means you successfully submitted your job to the supercomputer and the reductions are in progress.  You can see the log of the 
reductions in the file "reductionlrs2daily.oXXXXXX" where the XXXXXX is the job number as shown above in the line 
"Submitted batch job 900134".  The reductions should finish in 20 minutes or so depending on computer availability 
and number of exposures of the target.  The simplest way to see the effectiveness of the reduction is look at the source extraction
information in the log.

```
cat reductionlrs2daily.oXXXXXX | grep source
```

If you would like more flexibility in your batch processing, you can always edit "rlrs2_daily" to run any four reduction 
call you may want and submit the job manually with:
```
sbatch rgeneral_lrs2.slurm
```

The reductions will be in "LRS2/ORPHANS" for reductions before 2018/07/01 and in "LRS2/PROGRAM-ID" for reductions after this date.
The standard stars will be in "LRS2/STANDARDS" and the calibrations used are in "LRS2/CALS".  

## Code Description
Panacea is a general integral field unit (IFU) spectroscopic reduction tool tailored specifically for the Hobby Eberly Telescope (HET).
The code is primarily used for reducing science data from the LRS2 and VIRUS spectrographs.  
Before we dive into the details of the reduction code it is useful to visualize the mapping between 
raw data and on sky fiber layout for both instruments. 

### LRS2 Layout
The LRS2 instrument has two spectrographs each with two arms (LRS2-B: UV and Orange, LRS2-R: Red and FarRed).  
<p align="center">
  <img src="images/lrs2_mapping_visual.png" width="850"/>
</p>

### Bias Subtraction
The first step in Panacea's reduction is to measure the bias pedestal in the overscan region for each amplifier.  We subtract the pedestal row by row, excluding the first column in the overscan region using the remaining 31 or 63 pixels (2x1 or 1x1 binning, respectively) in a given row.  After the bias pedestal is removed, there remains a bias structure and excess charge that still needs to be removed. Each night, twenty bias exposures are taken, which are insufficient to accurately measure the low-level structure, so we use 100 consecutive frames over 5-6 nights, which is a compromise between sufficient statitics and minimal evolution of the structure from the passage of time.  The evolution of the bias structure depends on temperature control and the intricacies of each amplifier's controller. 

The quality of the bias subtraction is demonstrated in two ways: subtracting consecutive master bias frames from one another and looking at "blank" areas of science frames for residual counts and structure.  

<p align="center">
  <img src="images/bias_subtract.png" width="850"/>
</p>

<p align="center">
  <img src="images/trace_fibmodel.png" width="850"/>
</p>



## Frequently Asked Questions

Q: Are the wavelength units in vacuum or air?

A: Air

## Authors

* Greg Zeimann, UT Austin

## NOTE
* COPYRIGHTS from astropy, free software foundation were used
* cosmics.py is a copy from Malte Tewes and Pieter van Dokkum's code available online: http://obswww.unige.ch/~tewes/cosmics_dot_py/
