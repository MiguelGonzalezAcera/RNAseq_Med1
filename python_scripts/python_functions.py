import os
import logging
import glob

def run_command(command, logger):
    """Runs a command. Raises error if it misses
    """
    try:
        os.system(command)
    except ErrorCommandRun as errcomm:
        print(errcomm)
        logger.error(f"Command:\n\t{command}\ncould not run.")

def create_logger(logpath, name='logger'):
    """Create a logger

    args: name (str): name of logger

    returns: logger (obj): logging.Logger instance
    """
    logger = logging.getLogger(name)
    fmt = logging.Formatter('%(asctime)s - %(name)s -'
                            ' %(levelname)s -%(message)s')
    hdl = logging.FileHandler(logpath)
    hdl.setFormatter(fmt)

    logger.addHandler(hdl)

    return logger

def list_files_dir(dirpath, ext=""):
    """List files in a directory. Can specify an extension.
    """

    filelist = glob.glob(f"{dirpath}/{ext}")
    if len(filelist) == 0:
        raise ValueError(f'Directory {filepath} is empty.')
    return filelist
