'''General utils for MINT.'''

import os, subprocess
from Mint2.ConfigFile import ConfigFile, set_default_config
from ROOT import TVector3, DalitzEvent, FlexiFastAmplitudeIntegrator, FitAmpSum, DalitzEventPattern, \
    DalitzEventList, DiskResidentEventList
import ROOT
from array import array

def run_job(exe, workingdir, configs = [], parameters = {}, stdout = 'stdout', stderr = 'stderr') :
    '''Run a MINT executable with the given config files/parameters.'''

    if configs and isinstance(configs, (str, ConfigFile)) :
        configs = [configs]

    try :
        os.makedirs(workingdir)
    except OSError :
        if not os.path.exists(workingdir) :
            raise

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

def run_job_main() :
    '''Run a MINT job from a main script.'''

    from argparse import ArgumentParser
    import sys

    parser = ArgumentParser()
    parser.add_argument('exe', help = 'Executable to call.')
    parser.add_argument('workingdir', help = 'Name of the dataset to generate')
    parser.add_argument('--configs', nargs = '+', help = 'Config files to use')

    args, remainder = parser.parse_known_args()
    variableslists = {}
    for arg in remainder :
        if arg.startswith('--') :
            varargs = []
            variableslists[arg[2:]] = varargs
            continue
        varargs.append(arg)

    sys.exit(run_job(exe = args.exe, workingdir = args.workingdir, configs = args.configs, parameters = variableslists))

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
    try :
        os.makedirs(outputdir)
    except OSError :
        if not os.path.exists(outputdir) :
            raise

    integsdir = os.path.join(topdir, 'integrators')
    config['integratorsDirectory'] = [integsdir]

    datadir = os.path.join(mintdatadir, name)
    try :
        os.makedirs(datadir)
    except OSError :
        if not os.path.exists(datadir) :
            raise

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

def get_fit_fractions(fname, precision = 1e-2):
    '''Get the fit fractions for each resonance from the model in the given file.'''
    
    set_default_config(fname)
    config = ConfigFile(fname)
    pattern = DalitzEventPattern(*map(int, config['Event Pattern']))
    amps = FitAmpSum(pattern)
    integ = FlexiFastAmplitudeIntegrator(pattern, amps)
    integ.setPrecision(precision)
    integ.doFractions()
    return {frac.name() : frac.frac() for frac in integ.getFractions()}

def get_amp_scales(fname, targetfracs, normchannel, fout = None, precision = 1e-2):
    '''Get the scale factors for the amplitudes from the fit fractions, given the target fractions.
    Optionally save updated config to a new file.'''
    
    fracs = get_fit_fractions(fname, precision)
    scales = {normchannel : 1.}
    for name, frac in fracs.items():
        if name == normchannel:
            continue
        rexp = targetfracs[name]/targetfracs[normchannel]
        robs = frac/fracs[normchannel]
        scale = (rexp/robs)**.5
        scales[name] = scale
    if not fout:
        return scales
    config = ConfigFile(fname)
    for name, scale in scales.items():
        for suff in 'Re', 'Im':
            val = float(config[name + '_' + suff][1])
            val *= scale
            config[name + '_' + suff][1] = str(val)
    config.write_file(fout)
    return scales

def smear_momenta(evt, preslist):
    pmother = ROOT.TLorentzVector(0, 0, 0, 0)
    momenta = ROOT.vector('TLorentzVector')()
    for i in xrange(1, evt.eventPattern().size()):
        p = evt.p(i)
        m = p.M()
        ptot = p.P()
        psmeared = ptot + preslist[i-1](ptot)
        pscale = psmeared/ptot
        psmeared = ROOT.TLorentzVector(p.Px()*pscale, p.Py()*pscale, p.Pz()*pscale, (m**2+psmeared**2)**.5)
        pmother += psmeared
        momenta.push_back(psmeared)
    momenta.insert(momenta.begin(), pmother)
    return DalitzEvent(evt.eventPattern(), momenta)

class MomentumResolution(object):
    '''Callable to generate momentum smearing values.'''
    def __init__(self, c0 = 0.005, c1 = 0.01/100e3, seed = 0):
        '''Resolution generated is (c0 + c1*p)*p. Defaults are for charged tracks.'''
        self.c0 = c0
        self.c1 = c1
        self.rndm = ROOT.TRandom3(seed)

    def resolution(self, p):
        return (self.c0 + self.c1*p)*p        

    def __call__(self, p):
        return self.rndm.Gaus(0, self.resolution(p))

def add_resolutions_and_efficiencies(fname, fout, plab = (lambda : 0.), timeeff = (lambda t : True), 
                                     phasespaceeff = (lambda evt : True), timeres = (lambda t : 0.),
                                     pres = (lambda p : 0.), pres0 = (lambda p : 0.)):
    '''Add efficiencies and resolutions to an existing DalitzEventList.'''
    fin = ROOT.TFile.Open(fname)
    tdalitz = fin.Get('DalitzEventList')
    
    newlist = DiskResidentEventList(fout, 'recreate')
    ttime = ROOT.TTree('TimeEventList', 'TimeEventList')
    t = array('f', [0])
    tsmeared = array('f', [0])
    tag = array('i', [0])
    ttime.Branch('decaytime', t, 'decaytime/F')
    ttime.Branch('smeareddecaytime', tsmeared, 'smeareddecaytime/F')
    ttime.Branch('tag', tag, 'tag/I')
    
    tdalitz.GetEntry(0)
    evt = DalitzEvent()
    evt.fromNtuple(tdalitz)

    pattern = DalitzEventPattern(evt.eventPattern())
    preslist = []
    for i in xrange(1, pattern.size()):
        if pattern[i].charge() == '0':
            preslist.append(pres0)
        else:
            preslist.append(pres)

    for i in xrange(tdalitz.GetEntries()):
        tdalitz.GetEntry(i)
        evt = DalitzEvent()
        evt.fromNtuple(tdalitz)
        t[0] = tdalitz.decaytime
        tag[0] = tdalitz.tag
        if not (ROOT.gRandom.Rndm() < timeeff(t[0]) and ROOT.gRandom.Rndm() < phasespaceeff(evt)):
            continue
        tsmeared[0] = t[0] + timeres(t[0])
        p = ROOT.TVector3(0., 0., plab())
        evt.setMothers3Momentum(p)
        smearedevt = smear_momenta(evt, preslist)
        # ismear = 1
        # while not smearedevt.kinematicallyAllowed() and ismear < 1000:
        #     smearedevt = smear_momenta(evt, preslist)
        #     ismear += 1
        # if ismear == 1000:
        #     continue
        newlist.Add(smearedevt)
        ttime.Fill()
    print 'selected', ttime.GetEntries(), '/', tdalitz.GetEntries()
    ROOT.gDirectory.Get('DalitzEventList').AddFriend(ttime)
    ttime.Write()
    newlist.save()
    newlist.Close()
    fin.Close()

def event_loop(tree):
    '''Loop over events in the tree and build DalitzEvents from them.'''
    for i in xrange(tree.GetEntries()):
        tree.GetEntry(i)
        evt = DalitzEvent()
        evt.fromTree(tree)
        yield evt
