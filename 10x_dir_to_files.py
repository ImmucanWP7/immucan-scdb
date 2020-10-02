
import logging
import glob 
import os 
import re 
import shutil
import sys

logging.basicConfig()
log = logging.getLogger()
log.setLevel(logging.DEBUG)

directories_submitted_as_commandline_arguments = sys.argv[1:]
for directory_to_act_on in directories_submitted_as_commandline_arguments:
    log.debug(f"Acting on directory {directory_to_act_on}...")
    assert(not "/" in directory_to_act_on)

    all_filenames = glob.glob(f"{directory_to_act_on}/*")
    log.debug(f"Found these files: {all_filenames}")

    dir_name = os.path.split(os.path.abspath(directory_to_act_on))[1]
    separator = "_"

    for filename in all_filenames:
        log.debug(f"Moving {filename}...")
        cropped_filename = filename[len(directory_to_act_on)+1:]
        shutil.move(filename, dir_name+separator+cropped_filename)

    os.rmdir(directory_to_act_on)
