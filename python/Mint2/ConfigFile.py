'''A class to access parameters in a MINT config file.'''

class ConfigFile(object) :
    '''A class to access parameters in a MINT config file.'''
    
    def __init__(self, *fnames, **parameters) :
        '''Constructor, takes names of config files to use and/or the keyword arguments of parameters
        to use.'''

        self.fnames = []
        self.parameters = {}

        for fname in fnames :
            self.add_config_file(fname)
        self.parameters.update(parameters)

    def add_config_file(self, fname) :
        '''Add the parameters in the given config file.'''

        if isinstance(fname, ConfigFile) :
            self.add(fname)
            return

        if not os.path.exists(fname) :
            raise OSError('ConfigFile.add_config_file: Could not find file "{0}"!'.format(fname))

        self.fnames.append(fname)
        with open(fname) as configFile :
            for line in configFile :
                line = line.split()
                if len(line) == 0 or line[0].startswith('*') :
                    continue 
                self.parameters[line[0]] = line[1:]
                    
    def add(self, other) :
        '''Combine this config with another.'''
        self.parameters.update(other.parameters)
        self.fnames += other.fnames

    def write_file(self, fname, namewidth = 20, valwidth = 5) :
        '''Write the configuration to a file.'''

        with open(fname, 'w') as f :
            for par, vals in sorted(self.parameters.items()) :
                line = par.ljust(namewidth)
                line += ' '.join(repr(v).ljust(valwidth) for v in vals)
                f.write(line + '\n')
