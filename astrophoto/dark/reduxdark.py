import glob

import numpy as np
from astropy.io import fits

from astrophoto.io import printf, println, prompt, COMBINE_MEDIAN, prompt_choices, COMBINE_MEAN
from astrophoto.io.reduxio import saveandquit
from astrophoto.viz import plot_frame


def calib_dark(root_dir):

    file_pattern = root_dir + "/" + prompt("File pattern (glob pattern)", "dark_fpattern")
    file_list = glob.glob(file_pattern)

    texp = prompt("Exposure Time (sec) (float)", "dark_texp", float)

    data_list = []
    texp_list = []
    for f in file_list:

        hdr = fits.getheader(f)
        texp_f = hdr["EXPTIME"]

        if texp == 0.0 or texp_f == texp:
            data_list.append(fits.getdata(f).astype(np.float64))
            texp_list.append(texp_f)

    if len(data_list) == 0:
        printf("No files loaded.")
        saveandquit(20)

    darks = np.array(data_list)
    texps = np.array(texp_list)

    is_current = prompt("Calculate dark frame (0) or dark current (1)?", "dark_is_current", int)
    if is_current:
        bias_fname = prompt(
            "Master bias file for calculating dark current (path)?", "dark_bias")
        bias_file = root_dir + "/" + bias_fname
        bias = fits.getdata(bias_file).astype(np.float64)

        darks -= bias
        darks /= texps[..., None, None]

    combine_mode = prompt_choices(
        "How to combine frames (num)?",
        """
        (1) Median
        (2) Mean
        """, "dark_combine", int)

    if combine_mode == COMBINE_MEDIAN:
        dark = np.median(darks, axis=0)
    elif combine_mode == COMBINE_MEAN:
        dark = np.mean(darks, axis=0)
    else:
        printf("Invalid frame combination mode")
        saveandquit(10)

    plot_frame(dark)

    out_file = root_dir + "/" + prompt(
        "Output file for master dark frame (path)?", "dark_fout")
    fits.writeto(out_file, dark, overwrite=True)

    printf("Master dark frame saved.")
    println()

    return dark