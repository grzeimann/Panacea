# -*- coding: utf-8 -*-
"""
Argument Parser
------------------
Built for the VIRUS instrument as well as LRS2 on HET

@author: gregz
"""

import argparse as ap
import textwrap
import glob
import os.path as op
import config
from amplifier import Amplifier

def parse_args(argv=None):
    """Parse the command line arguments

    Parameters
    ----------
    argv : list of string
        arguments to parse; if ``None``, ``sys.argv`` is used

    Returns
    -------
    Namespace
        parsed arguments
    """
    description = textwrap.dedent('''Panacea ''')
                     
    parser = ap.ArgumentParser(description=description,
                            formatter_class=ap.RawTextHelpFormatter)
                            
    parser.add_argument("-rt","--reduce_twi", 
                        help='''Reduce Twighlight frames for calibration''',
                        action="count", default=0)

    parser.add_argument("-rs","--reduce_sci", 
                        help='''Reduce Science frames''',
                        action="count", default=0)

    parser.add_argument("-c","--combine_reductions", 
                        help='''Produce Fe, CuFe, S products''',
                        action="count", default=0)
                                                                                                
    parser.add_argument("--ifuslot", nargs='?', type=str, 
                        help='''Single IFUSLOT for processing. [REQUIRED]
                        Ex: "075".''', default = None)

    parser.add_argument("--instr", nargs='?', type=str, 
                        help='''Instrument to process. 
                        Default: "virus"
                        Ex: "camra" for lab data,
                            "lrs2" for lrs2.''', default = "virus")
                            
    parser.add_argument("--instr_side", nargs='?', type=str, 
                        help='''Instrument side to process. 
                        Default: "blue"
                        Ex: "blue" for LRS2B,
                            "red" for LRS2R.''', default = "blue")
                                                
    parser.add_argument("-sd","--scidir_date", nargs='?', type=str,
                        help='''Science Directory Date.     [REQUIRED, if --reduce_sci]
                        Ex: \"20160412\"''', default=None)

    parser.add_argument("-so","--scidir_obsid", nargs='?', type=str,
                        help='''Science Directory ObsID.    [REQUIRED, if --reduce_sci]
                        Ex: \"3\" or \"102\"''', default=None)
                        
    parser.add_argument("-se","--scidir_expnum", nargs='?', type=str,
                        help='''Science Directory exposure number.
                        Ex: \"1\" or \"05\"''', default=None) 

    parser.add_argument("-td","--twidir_date", nargs='?', type=str,
                        help='''Twi Directory Date.     [REQUIRED, if --reduce_twi]
                        Ex: \"20160412\"''', default=None)

    parser.add_argument("-to","--twidir_obsid", nargs='?', type=str,
                        help='''Twi Directory ObsID.    [REQUIRED, if --reduce_twi]
                        Ex: \"3\" or \"102\"''', default=None)
                        
    parser.add_argument("-te","--twidir_expnum", nargs='?', type=str,
                        help='''Twi Directory exposure number.
                        Ex: \"1\" or \"05\"''', default=None) 

    parser.add_argument("-skd","--skydir_date", nargs='?', type=str,
                        help='''Sky Directory Date.     [REQUIRED, --use_other_sky]
                        Ex: \"20160412\"''', default=None)

    parser.add_argument("-sko","--skydir_obsid", nargs='?', type=str,
                        help='''Sky Directory ObsID.    [REQUIRED, --use_other_sky]
                        Ex: \"3\" or \"102\"''', default=None)
                        
    parser.add_argument("-ske","--skydir_expnum", nargs='?', type=str,
                        help='''Sky Directory exposure number.
                        Ex: \"1\" or \"05\"''', default=None) 

    parser.add_argument("-q","--quiet", help='''Turn off logging.''',
                        action="count", default=0)
                          
    args = parser.parse_args(args=argv)

    args = read_in_raw(args, parser)
    
    return args
    
def read_in_raw(args, parser):
    # Check that the arguments are filled
    if args.ifuslot:
        args.ifuslot = "%03d" %int(args.ifuslot)
    else:
        msg = 'No IFUslot was provided.'
        parser.error(msg)              
    
    if args.instr.lower() == 'virus':
        instr = 'virus_'
    if args.instr.lower() == 'lrs2':
        if args.instr_side.lower() == 'blue':
            instr = 'lrs2b_'
        if args.instr_side.lower() == 'red':
            instr = 'lrs2r_'
    if args.instr.lower() == 'virusw':
        instr = 'virusw_' 
    
    args.kwargs = {}
    for param in config.param_dict:
        args.kwargs[param] = getattr(config, instr+config.param_dict[param])
    for param in config.param_amp_dict:
        setattr(args, param, getattr(config, 
                                     instr+config.param_amp_dict[param])) 
    for con in config.config_dict:
        args.kwargs[con] = getattr(config, config.config_dict[con]) 
    
    args.use_other_sky = config.use_other_sky
    args.refit_fiber_to_fiber = config.refit_fiber_to_fiber
    args.ifucen_fn = getattr(config,instr+'fn')
    args.side_dict = config.Side_dict
    args.disp = getattr(config,instr+'di')
    args.scale = getattr(config,instr+'cs')
    args.nfibers = getattr(config, instr+'nf')
    
    if args.quiet:
        args.kwargs['verbose'] = False

    labels = ['dir_date', 'dir_obsid', 'dir_expnum']
    observations=[]
    if args.combine_reductions and not args.reduce_sci:
        observations.append('sci')
        args.sci_operations = {}
        att_list = [attr for attr in dir(config) if not callable(getattr(config, attr)) 
                                                    and not attr.startswith("__") 
                                                 if attr[:4] == 'sci_']
        for att in att_list:
            args.sci_operations[att[4:]] = getattr(config, att)
    if args.reduce_sci:
        observations.append('sci')
        if not args.reduce_twi:
            observations.append('twi')
        args.sci_operations = {}
        att_list = [attr for attr in dir(config) if not callable(getattr(config, attr)) 
                                                    and not attr.startswith("__") 
                                                 if attr[:4] == 'sci_']
        for att in att_list:
            args.sci_operations[att[4:]] = getattr(config, att)
    if args.reduce_twi:
        observations.append('twi')
        args.twi_operations = {}
        att_list = [attr for attr in dir(config) if not callable(getattr(config, attr)) 
                                                    and not attr.startswith("__") 
                                                 if attr[:4] == 'twi_']
        for att in att_list:
            args.twi_operations[att[4:]] = getattr(config, att) 
    if args.use_other_sky:
        observations.append('sky')
    for obs in observations:
        amp_list = []
        for label in labels[:2]:
            getattr(args, obs+label)
            if getattr(args, obs+label) is None:
                msg = '%s%s was not provided' %(obs, label)
                parser.error(msg) 
            else:
                setattr(args, obs+label, 
                        getattr(args, obs+label).replace(" ", "").split(','))
        if getattr(args, obs+labels[2]) is not None:
            setattr(args, obs+labels[2], 
                    getattr(args, obs+labels[2]).replace(" ", "").split(','))
        for date in getattr(args, obs+labels[0]):
            for obsid in getattr(args, obs+labels[1]):
                if getattr(args, obs+labels[2]) is not None:   
                    for expnum in getattr(args, obs+labels[2]):
                        folder = op.join(date, 
                                         args.instr, 
                                         "{:s}{:07d}".format(args.instr,
                                         int(obsid)), 
                                         "exp{:02d}".format(int(expnum)),
                                         args.instr)
                        files = sorted(glob.glob(op.join(config.rootdir, folder, 
                                                         '*_%s*' %args.ifuslot)))
                        path = op.join(config.output, folder)
                        for fn in files:
                            amp_list.append(Amplifier(fn, path, **args.kwargs))
                            amp_list[-1].type = obs
                            for Amp in config.Amps: 
                                if (amp_list[-1].amp == Amp) or (amp_list[-1].amp ==config.Amp_dict[Amp][0]):
                                    for param in config.param_amp_dict:
                                        setattr(amp_list[-1], param, 
                                                getattr(args, param)[Amp])


                else:
                    folder = op.join(date, args.instr,
                                     "{:s}{:07d}".format(args.instr, 
                                                         int(obsid)))
                    files = sorted(glob.glob(op.join(config.rootdir, folder, '*', 
                                                     args.instr, '*_%s*' %args.ifuslot)))
                    for fn in files:
                        expn = op.basename(op.dirname(op.dirname(fn)))
                        path = op.join(config.output, folder, expn, args.instr)
                        amp_list.append(Amplifier(fn, path, **args.kwargs))
                        amp_list[-1].type = obs
                        for Amp in config.Amps: 
                                if (amp_list[-1].amp == Amp) or (amp_list[-1].amp ==config.Amp_dict[Amp][0]):
                                    for param in config.param_amp_dict:
                                        setattr(amp_list[-1], param, 
                                                getattr(args, param)[Amp])
        setattr(args,obs+'_list', amp_list)                    

    return args