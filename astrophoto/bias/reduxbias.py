import glob

import numpy as np
from astropy.io import fits

from astrophoto.io.reduxio import printf, println, prompt, COMBINE_MEDIAN, prompt_choices, COMBINE_MEAN, save_and_quit
from astrophoto.viz.plot import plot_frame


def calib_bias(root_dir):
    # save inputs for later
    file_pattern = root_dir + "/" + prompt("File pattern (glob pattern)", "bias_fpattern")
    file_list = glob.glob(file_pattern)

    data_list = []
    for f in file_list:
        data_list.append(fits.getdata(f).astype(np.float64))

    if len(data_list) == 0:
        printf("No files loaded.")
        save_and_quit(20)

    biases = np.array(data_list)

    combine_mode = prompt_choices(
        "How to combine frames (num)?",
        """
        (1) Median
        (2) Mean
        """, "bias_combine", int)

    if combine_mode == COMBINE_MEDIAN:
        bias = np.median(biases, axis=0)
    elif combine_mode == COMBINE_MEAN:
        bias = np.mean(biases, axis=0)
    else:
        printf("Invalid frame combination mode")
        save_and_quit(10)

    plot_frame(bias)

    out_file = root_dir + "/" + prompt(
        "Output file for master bias frame (path)?", "bias_fout")
    fits.writeto(out_file, bias, overwrite=True)

    printf("Master bias frame saved.")
    println()

    return bias