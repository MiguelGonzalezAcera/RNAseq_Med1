import subprocess
import logging
import glob

def run_command(command):
    """Runs a command. Raises error if it misses
    """
    logging.info(f'Running command: {command}')
    for cmd in command.split('; '):
        output = subprocess.run(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stout = output.stdout.decode('utf-8')
        error = output.stderr.decode('utf-8')
        if output.returncode == 1:
            logging.error(f'{cmd}\n\n{error}')
            raise ValueError(f'Error in command: {cmd}\n\n{error}')
        elif output.returncode == 0:
            logging.info(stout)
            logging.info(error)

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
        raise ValueError(f'Directory {dirpath} is empty.')
    return filelist
