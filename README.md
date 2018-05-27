# Panacea 
This package is intended to be the base reduction pipeline for VIRUS and LRS2 at the Hobby Eberly Telescope.  Instructions for installation and use at the Texas Advanced Computing Center (TACC) are below.

## Getting Started with LRS2
To begin on TACC, point to the common python environment. In your home ".bashrc" file, add the following line at the bottom:
```
export PATH=”/home/00115/gebhardt/anaconda2/bin:/work/03946/hetdex/maverick/bin:$PATH”
```

Then move to your work directory and clone Panacea: 
```
cdw
git clone https://github.com/grzeimann/Panacea.git
```

The next step is to generate the necessary set of scripts for your target:
```
python Panacea/build_panacea_call.py --start_date 20180515 --date_length 1 --rootdir /work/03946/hetdex/maverick --instrument lrs2 --side blue --target sdss
```

The following scripts are generated from that call and printed to screen:
```
sbatch rtwi_blue_1.slurm
sbatch rsci_blue_1.slurm
sbatch rstd_blue_1.slurm
sbatch rresponse_blue_1.slurm
```

### Authors

* Greg Zeimann, UT Austin
* Karl Gebhardt, UT Austin
* Alex Hagen, Penn State

#### NOTE
* COPYRIGHTS from astropy, free software foundation were used
* cosmics.py is a copy from Malte Tewes and Pieter van Dokkum's code available online: http://obswww.unige.ch/~tewes/cosmics_dot_py/
