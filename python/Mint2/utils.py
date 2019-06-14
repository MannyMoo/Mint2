'''General utils for MINT.'''

import os, subprocess
from Mint2.ConfigFile import ConfigFile
from ROOT import TVector3, DalitzEvent
import ROOT

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

def gen_time_dependent(name, integratorsdir, mintdatadir, configs = [], number = None, zfill = 3, **parameters) :
    '''Generate time dependent MINT MC with genTimeDependent.exe.'''

    if not configs and not parameters :
        raise ValueError('Must give at least one config file or pass the config via **parameters!')

    if isinstance(configs, str) :
        configs = [configs]
    config = ConfigFile(*configs, **parameters)

    topdir = os.path.join(integratorsdir, name)
    if number == None :
        number = 0
        while os.path.exists(os.path.join(topdir, str(number).zfill(zfill))) :
            number += 1

    outputdir = os.path.join(topdir, str(number).zfill(zfill))
    if not os.path.exists(outputdir) :
        os.makedirs(outputdir)
    integsdir = os.path.join(topdir, 'integrators')
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
    parser.add_argument('--number', help = 'Number of the job (used as the name of the working directory).', default = None)
    parser.add_argument('--zfill', help = 'How many zeros to pad the job number with', default = 3)

    args, remainder = parser.parse_known_args()
    variableslists = {}
    for arg in remainder :
        if arg.startswith('--') :
            varargs = []
            variableslists[arg[2:]] = varargs
            continue
        varargs.append(arg)

    sys.exit(gen_time_dependent(name = args.name, configs = args.configs, integratorsdir = args.integratorsdir,
                                mintdatadir = args.mintdatadir, number = args.number, zfill = args.zfill,
                                **variableslists))

def three_body_event(pattern, s13, s23) :
    '''Create a 3-body DalitzEvent with the given DalitzEventPattern, s13 and s23.'''

    if pattern.size() != 4 :
        raise ValueError('The given DalitzEventPattern is for a '
                         '{0}-body mode, not 3-body!'.format(pattern.size()-1))

    m0 = pattern[0].mass()
    m1 = pattern[1].mass()
    m2 = pattern[2].mass()
    m3 = pattern[3].mass()

    px1 = 0.
    py1 = 0.
    pz1 = 0.
    
    px2 = 0.
    py2 = 0.
    pz2 = 0.

    px3 = 0.
    py3 = 0.
    pz3 = 0.

    #Calculate energies of pi+, pi- from known parameters
    E1 = ( m0**2 + m1**2 - s23 ) / (2 * m0)
    E2 = ( m0**2 + m2**2 - s13 ) / (2 * m0)

    #Get E3 from conservation of energy
    E3 = m0 - E1 - E2

    #Reject events where kinematics are impossible
    pz1Sq = E1**2 - m1**2
    if pz1Sq < 0 :
        return None

    #Get pz1 from mass/energy relation (co-ordinate system defined so that px1 = py1 = 0)
    pz1 = pz1Sq**0.5

    #Calculate pz3 from conservation of momentum 
    pz3 = ( E2**2 - m2**2 - E3**2 + m3**2 - pz1**2 ) / ( 2*pz1 )

    #Use pz3 to get py3 (co-ordinate system designed so that px3 = 0)
    py3Sq =  E3**2 - m3**2 - pz3**2

    #Reject events where kinematics are impossible
    if py3Sq < 0 :
        return None

    py3 = py3Sq**0.5
    #Conservation of momentum => py2 = -py3 and pz2 = -(pz1 + pz3)
    py2 = -1 * py3
    pz2 = -1 * (pz1 + pz3)

    momenta = ROOT.vector('TVector3')()
    momenta.push_back(TVector3(0., 0., 0.)) # mother momentum
    momenta.push_back(TVector3(px1, py1, pz1))
    momenta.push_back(TVector3(px2, py2, pz2))
    momenta.push_back(TVector3(px3, py3, pz3))

    return DalitzEvent(pattern, momenta)
