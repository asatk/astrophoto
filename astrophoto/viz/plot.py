from matplotlib import pyplot as plt
import numpy as np

def plot_frame(frame: np.ndarray, title: str=None, sigma:float=1.0, show: bool=False):
    med = np.median(frame)
    mad = np.median(np.abs(frame - med))
    vmin = med - sigma * mad
    vmax = med + sigma * mad

    fig, ax = plt.subplots()

    im = ax.imshow(frame, origin="lower", vmin=vmin, vmax=vmax, cmap="grey")
    cb = plt.colorbar(im)

    if title is not None:
        ax.set_title(title)

    if show:
        plt.show()

    return ax