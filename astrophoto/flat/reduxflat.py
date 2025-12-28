import glob

import numpy as np
from astropy.io import fits

from astrophoto.io import printf, println, prompt, COMBINE_MEDIAN, prompt_choices, COMBINE_MEAN
from astrophoto.io.reduxio import saveandquit
from astrophoto.viz import plot_frame


def calib_flat(root_dir):

    file_pattern = root_dir + "/" + prompt("File pattern (glob pattern)", "flat_fpattern")
    file_list = glob.glob(file_pattern)

    band = prompt("Instrument Filter (str)", "flat_band", str)

    data_list = []
    texp_list = []

    for f in file_list:

        hdr = fits.getheader(f)
        band_f = hdr["FILTER"]

        if band_f == band:
            data_list.append(fits.getdata(f).astype(np.float64))
            texp_list.append(hdr["EXPTIME"])

    if len(data_list) == 0:
        printf("No files loaded.")
        saveandquit(20)

    flats = np.array(data_list)
    texps = np.array(texp_list)

    bias_fname = prompt(
        "Master bias file (path) ", "flat_bias")
    bias_file = root_dir + "/" + bias_fname
    bias = fits.getdata(bias_file).astype(np.float64)

    flats -= bias

    dark_fname = prompt(
        "Master dark current file (path) ", "flat_dark")
    if dark_fname != "":
        dark_file = root_dir + "/" + dark_fname
        dark = fits.getdata(dark_file).astype(np.float64)
        darks = np.array([t * dark for t in texps])
        flats -= darks

    # try this... see if combining is better after normalizing flats
    # flats /= np.mean(flats, axis=0)

    combine_mode = prompt_choices(
        "How to combine frames (num)?",
        """
        (1) Median
        (2) Mean
        """, "flat_combine", int)

    if combine_mode == COMBINE_MEDIAN:
        flat = np.median(flats, axis=0)
    elif combine_mode == COMBINE_MEAN:
        flat = np.mean(flats, axis=0)
    else:
        printf("Invalid frame combination mode")
        saveandquit(10)

    # normalize combined flat
    flat /= np.mean(flat)

    plot_frame(flat)

    out_file = root_dir + "/" + prompt(
        f"Output file for master {band} flat frame (path)?", "flat_fout")
    fits.writeto(out_file, flat, overwrite=True)

    printf(f"Master {band} flat frame saved.")
    println()

    return flat