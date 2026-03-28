import numpy as np
import matplotlib.pyplot as plt
import os


def parse_infile(infile, verbose=True):
    """
    Take a string (folder/file), or 2 element sequence ([folder, file]), and
    return (folder, file).
    """

    if type(infile) == str:
        folder, file = os.path.split(infile)
    elif type(infile) in [list, tuple]:
        if len(infile) != 2:
            raise Exception('The infile must be a string or a length 2 sequence')
        folder, file = infile
    else:
        raise Exception('The infile must be a string or a length 2 sequence')

    return folder, file

def drop_extension(fullpath):

    path, file = os.path.split(fullpath)

    filestub = os.path.splitext(file)[0]

    return os.path.join(path, filestub)

def swap_extension(fullpath, new_ext):

    fullpath = drop_extension(fullpath)

    return fullpath + new_ext

def write_dict_as_str(filename, mydictionary, security='low'):
    """
    Writes a dict as a text file. Import requires eval so not particularly safe. 
    """
    
    with open(filename, 'w') as f:
        print(mydictionary, file=f)
    
def read_dict_from_str(filename, security='low'):
    """
    Reads text file produced by write_dict_as_str. Use eval so not particularly safe operation.
    """
    
    with open(filename, 'r') as f:
        line = f.readline().strip()
        mydictionary = eval(line)
    
    return mydictionary

