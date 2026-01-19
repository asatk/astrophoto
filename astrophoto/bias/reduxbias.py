from astropy.io import fits
import numpy as np

from astrophoto.core import combine
from astrophoto.io.reduxio import printf, println, prompt, load_files
from astrophoto.viz.plot import plot_frame



def calib_bias(root_dir):

    biases = load_files(root_dir)

    bias = combine(biases)

    plot_frame(bias, show=True)

    save_file(root_dir, bias)

    return bias

def save_file(root_dir: str, data: np.ndarray):

    out_file = root_dir + "/" + prompt(
        "Output file for frame (path)?", "fout")
    fits.writeto(out_file, data, overwrite=True)

    printf("Frame saved.")
    println()