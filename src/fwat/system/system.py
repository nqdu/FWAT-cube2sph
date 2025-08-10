import os 
import shutil 
import glob 

def sys_remove(pattern,verbose = False):
    """
    Mimics 'rm -rf' with wildcard support.
    Removes files, directories, and symlinks matching the pattern.
    """
    for path in glob.glob(pattern, recursive=True):
        if verbose:
            print(f"removing {path}")
        try:
            if os.path.islink(path) or os.path.isfile(path):
                os.remove(path)
            elif os.path.isdir(path):
                shutil.rmtree(path)
            elif not os.path.exists(path):
                pass
        except Exception as e:
            print(f"Failed to remove {path}: {e}")