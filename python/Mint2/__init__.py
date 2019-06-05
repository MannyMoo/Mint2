
def load_cpp_lib() :
    '''Load the C++ libary for this package into ROOT.'''
    import ROOT
    ROOT.gSystem.Load('libMint2Lib.so')
    ROOT.gSystem.Load('libMint2Dict.so')

# Uncomment the below if you want to automatically load your C++ libraries
load_cpp_lib()
