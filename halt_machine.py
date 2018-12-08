import utils_snpko as utils
import subprocess

logger = utils.logger


def possibly_halt(args):
    if args.halt:
        logger.info("Halting machine.")
        subprocess.call('sudo shutdown -h now', shell=True)
