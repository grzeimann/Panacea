"""
EXAMPLE BLOCK CODE FOR PANACEA PYSPARK ASSIMILATION
"""

# Import Pyspark
import pyspark

# The class amplifier loads the fits files (astropy.io.fits), picks out
# info from the header, performs operations on the fits as well as other
# calculated products, and is the heart of a reduction system called "panacea"
# which is a play on the instrument VIRUS.  This code also works with
# the instruments of LRS2 (on the HET) and VIRUS-W (on 2.7m at McD).
from amplifier import Amplifier

# Import dictionary of arguments for Amplifier specific to Virus
# These get set to variables in the Amplifier class and get carried around
# through functions as needed.  Thus, there is some redundancy here,
# which can be dealt with later if necessary.
# This makes it easy for a user to specify new parameters in the reduction
# without touching the main script.
from config_amplifier_virus import kwargs, date_start, date_end

# Call pyspark's SparkContext() with stand-in 'local[*]'
# This block should get replaced with a context specific load and
# will help maximize the parallelization.  More help needed here.
sc = pyspark.SparkContext('local[*]')


def getfiles(start_date, end_date, exposure_type='twi'):
    # GET FILES
    pass


def getAmplifier(filename_full_path):
    return Amplifier(filename_full_path, **kwargs)


def reduceTwilight(amp):
    # The Amplifier class has many functions that either manipulate
    # the original raw fits file, or calculate a new array/list
    # Thus the "amp" here grows in memory with a ceiling
    # and a step can be added after these operations that only keeps the
    # final desired values.

    operations = ['prepare_image', 'get_trace', 'get_fibermodel',
                  'get_wavelength_solution', 'get_fiber_to_fiber']

    for op in operations:
        getattr(amp, op)()

    return amp


def reduceScience(amp):
    # The Amplifier class has many functions that either manipulate
    # the original raw fits file, or calculate a new array/list
    # Thus the "amp" here grows in memory with a ceiling
    # and a step can be added after these operations that only keeps the
    # final desired values.

    operations = ['sky_subtraction', 'clean_cosmics', 'sky_subtraction',
                  'calculate_significance']

    for op in operations:
        getattr(amp, op)()

    return amp


def map_twi_to_sci(meta_twi, meta_sci, raw_sci):
    # TO DO
    pass


def main():
    # ##################### TWI REDUCTIONS ####################################
    # Return all file names (full path) of a certain type between dates
    # For example: 'twi', 'sci', 'zro', 'flt', or 'cmp'
    filelist_twi = getfiles(date_start, date_end, exposure_type='twi')

    # Initialize the amplifier class and make rdd with file names as key
    # and Amplifier class object as value
    raw_twi = sc.parallelize(filelist_twi).mapValues(getAmplifier)

    # Reduce the amplifiers for a twilight frame
    reduced_twi = raw_twi.mapValues(reduceTwilight).collect()
    # Collect the metadata about each Amplifier object for later pairing with
    # science frame
    metadata_twi = reduced_twi.map(lambda x: (x[1].date, x[1].amp,
                                              x[1].ifuslot, x[1].specid,
                                              x[1].ifuid))

    # ##################### SCI REDUCTIONS ####################################
    # The current mode is to write out certain twi products and the filename
    # connects to a given science frame.  All the science frame needs is the
    # path and it will load the file, extract the info, and assimilate with
    # the Amplifier class.  We can change that and keep it all in memory,
    # but I will have to enhance Amplifier for this new mode.
    filelist_sci = getfiles(date_start, date_end, exposure_type='sci')

    # Initialize the amplifier class and make rdd with file names as key
    # and Amplifier class object as value
    raw_sci = sc.parallelize(filelist_sci).mapValues(getAmplifier)
    metadata_sci = raw_sci.map(lambda x: (x[1].date, x[1].amp,
                                          x[1].ifuslot, x[1].specid,
                                          x[1].ifuid))

    # Need to "copy/apply" products from twi reduction to sci Amplifier
    # based on closest date for a match same amp, ifuslot, specid, ifuid
    intermediate_sci = map_twi_to_sci(metadata_twi, metadata_sci, raw_sci)
    reduced_sci = intermediate_sci.mapValues(reduceScience).collect()

    # This concludes the basic reduction
    # Information should then be stored into a database or pass on to more
    # advanced steps.

if __name__ == '__main__':
    main()
