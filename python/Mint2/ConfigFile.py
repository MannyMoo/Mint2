'''A class to access parameters in a MINT config file.'''

import os

def parse_line(line) :
    '''Parse a line from a config file, taking into account quotes.'''
    line = line.strip()
    if not line :
        return []
    if line.startswith('"') :
        istart = 1
        iend = line.find('"', 1)
        val = line[1:iend]
        remainder = line[iend+1:]
    else :
        splitline = line.split()
        val = splitline[0]
        remainder = ' '.join(splitline[1:])
    if remainder :
        return [val] + parse_line(remainder)
    return [val]

def quotes(val, width = 20) :
    '''Return a string version of the value enclosed in double quotes, left adjusted to
    the given width.'''
    return ('"' + str(val) + '"').ljust(width)

class ConfigFile(dict) :
    '''A class to access parameters in a MINT config file.'''
    
    def __init__(self, *fnames, **parameters) :
        '''Constructor, takes names of config files to use and/or the keyword arguments of parameters
        to use.'''

        self.fnames = []

        for fname in fnames :
            self.add_config_file(fname)
        self.update(parameters)

    def add_config_file(self, fname) :
        '''Add the parameters in the given config file.'''

        if isinstance(fname, ConfigFile) :
            self.add(fname)
            return

        fname = os.path.abspath(os.path.expandvars(fname))
        if not os.path.exists(fname) :
            raise OSError('ConfigFile.add_config_file: Could not find file "{0}"!'.format(fname))

        self.fnames.append(fname)
        with open(fname) as configFile :
            for line in configFile :
                line = line.strip()
                if not line or line.startswith('*') :
                    continue 
                line = parse_line(line)
                self[line[0]] = line[1:]
                    
    def add(self, other) :
        '''Combine this config with another.'''
        self.update(other)
        self.fnames += other.fnames

    def write_file(self, fname, namewidth = 30, valwidth = 6) :
        '''Write the configuration to a file.'''

        with open(fname, 'w') as f :
            for par, vals in sorted(self.items()) :
                line = quotes(par, namewidth) + ' '
                line += ' '.join(quotes(v, valwidth) for v in vals)
                f.write(line + '\n')

    
    def float(self, name) :
        '''Get a (the first) parameter value as a float.'''
        return float(self[name][0])

    def floats(self, name) :
        '''Get the parameter values as a list of floats.'''
        return map(float, self[name])
