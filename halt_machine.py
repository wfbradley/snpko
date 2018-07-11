import utils_snpko as utils
from version_snpko import __version__
import subprocess

logger = utils.logger

def possibly_halt(args):
    if args.halt:
        logger.info("Halting machine.")
        subprocess.call('sudo shutdown -h now',shell=True)
