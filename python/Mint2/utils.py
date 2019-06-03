'''General utils for MINT.'''

import os, subprocess
from Mint2.ConfigFile import ConfigFile

def run_job(exe, workingdir, *configs, **parameters) :
    '''Run a MINT executable with the given config files/parameters.'''

    pwd = os.getcwd()
    conf = ConfigFile(*configs, **parameters)
    conf.write_file('config.txt')
    stdout = open('stdout', 'w')
    stderr = open('stderr', 'w')
    retval = subprocess.call(app + ' < config.txt', shell = True, stdout = stdout, stderr = stderr)
    stdout.close()
    stderr.close()
    return retval

