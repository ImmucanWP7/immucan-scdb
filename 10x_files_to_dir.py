
import logging
import glob 
import os 
import re 
import shutil

logging.basicConfig()
log = logging.getLogger()
log.setLevel(logging.DEBUG)

log.debug("Searching for matrix files...")

all_matrix_filenames = glob.glob("*matrix.mtx")
log.debug(f"Found these matrix files: {all_matrix_filenames}")

for matrix_filename in all_matrix_filenames:
    m = re.search("^(.*)matrix.mtx$", matrix_filename)
    assert(m is not None) ## should find a hit

    hits = m.groups()
    assert(len(hits)==1) ## should be exactly one hit

    first_part = hits[0]

    log.debug(f"Creating directory {first_part}")
    os.makedirs(first_part)

    all_files_like_this = glob.glob(f"{first_part}*.*")
    log.debug(f"Will move the following files: {all_files_like_this}")

    for filename in all_files_like_this:
        cropped_filename = filename[len(first_part):]
        shutil.move(filename, first_part + "/" + cropped_filename)
