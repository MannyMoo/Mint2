'''General utils for MINT.'''

import os, subprocess
from Mint2.ConfigFile import ConfigFile
from ROOT import TVector3, TLorentzVector

def run_job(exe, workingdir, configs = [], parameters = {}, stdout = 'stdout', stderr = 'stderr') :
    '''Run a MINT executable with the given config files/parameters.'''

    if configs and isinstance(configs, (str, ConfigFile)) :
        configs = [configs]

    pwd = os.getcwd()
    os.chdir(workingdir)

    conf = ConfigFile(*configs, **parameters)
    conf.write_file('config.txt')
    fstdout = open(stdout, 'w')
    if stdout == stderr :
        fstderr = fstdout
    else :
        fstderr = open(stderr, 'w')
    retval = subprocess.call(exe + ' < config.txt', shell = True, stdout = fstdout, stderr = fstderr)
    fstdout.close()
    fstderr.close()
    os.chdir(pwd)
    return retval

def gen_time_dependent(name, integratorsdir, mintdatadir, configs = [], **parameters) :
    '''Generate time dependent MINT MC with genTimeDependent.exe.'''

    if not configs and not parameters :
        raise ValueError('Must give at least one config file or pass the config via **parameters!')

    if isinstance(configs, str) :
        configs = [configs]
    config = ConfigFile(*configs, **parameters)
    outputdir = os.path.join(integratorsdir, name)
    if not os.path.exists(outputdir) :
        os.makedirs(outputdir)
    integsdir = os.path.join(outputdir, 'integrators')
    config['integratorsDirectory'] = [integsdir]

    datadir = os.path.join(mintdatadir, name)
    if not os.path.exists(datadir) :
        os.makedirs(datadir)
    config['outputFileName'] = [os.path.join(datadir, os.path.split(config['outputFileName'][0])[1])]

    return run_job(exe = 'genTimeDependent.exe', workingdir = outputdir, configs = config,
                   stdout = 'stdout', stderr = 'stdout')

def gen_time_dependent_main(defaultconfigs, defaultintegratorsdir, defaultdatadir) :
    '''Main function to generate time dependent MINT MC.'''
    from argparse import ArgumentParser
    import sys

    parser = ArgumentParser()
    parser.add_argument('name', help = 'Name of the dataset to generate')
    parser.add_argument('--configs', nargs = '*', help = 'Config files to use', default = defaultconfigs)
    parser.add_argument('--integratorsdir', help = 'Integrators directory', default = defaultintegratorsdir)
    parser.add_argument('--mintdatadir', help = 'Directory to save generated data', default = defaultdatadir)
    
    args, remainder = parser.parse_known_args()
    variableslists = {}
    for arg in remainder :
        if arg.startswith('--') :
            varargs = []
            variableslists[arg[2:]] = varargs
            continue
        varargs.append(arg)

    sys.exit(gen_time_dependent(name = args.name, configs = args.configs, integratorsdir = args.integratorsdir,
                                mintdatadir = args.mintdatadir, **variableslists))

def three_body_event(pattern, s13, s23) :
    '''To be implemented.'''

    m0 = pattern[0].mass()
    m1 = pattern[1].mass()
    m2 = pattern[2].mass()
    m3 = pattern[3].mass()

    '''
    s13 = (e1 + e3)**2 - (px1 + px3)**2 - (py1 + py3)**2 - (pz1 + pz3)**2
    s23 = (e2 + e3)**2 - (px2 + px3)**2 - (py2 + py3)**2 - (pz2 + pz3)**2
    
    # set px3 = py3 = 0, py2 = 0
    s13 = ((m1**2 + px1**2 + py1**2 )**.5 + (m3**2 + pz3**2))**2 - px1**2 - py1**2 - (pz1 + pz3)**2
    s23 = ((e2 + e3)**2 - (px2 + px3)**2 - (py2 + py3)**2 - (pz2 + pz3)**2
    '''
    
