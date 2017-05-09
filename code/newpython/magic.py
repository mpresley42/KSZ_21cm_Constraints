import fnmatch as fnm
import re
import os 

def match_file(location,pattern):
    """Returns a file from the given location that matches 
       the given pattern."""
    for file in os.listdir(location):
        if fnm.fnmatch(file, pattern):
            return file
    return None

def match_files(location,pattern):
    """Returns a list of all the files from the given location 
       that match the given pattern."""
    files = []
    for f in os.listdir(location):
        if fnm.fnmatch(f, pattern):
            files.append(f)
    return files

def nums_from_string(f):
    """Uses regex magic to extract all of the numbers from a string
       and returns them as a list."""
    return [float(i) for i in re.findall("[-+]?\d+[\.]?\d*",f)] 
