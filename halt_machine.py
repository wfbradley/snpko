import utils_snpko as utils
import subprocess
import logging

logger = utils.logger


def possibly_halt(args):
    if args.halt:
        logger.info("Halting machine.")
        # Flush any log buffers
        logging.shutdown()
        subprocess.call('sudo shutdown -h now', shell=True)
