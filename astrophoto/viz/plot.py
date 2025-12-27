from matplotlib import pyplot as plt
import numpy as np

def plot_frame(frame, title: str=None, sigma:float=1.0):
    med = np.median(frame)
    mad = np.median(np.abs(frame - med))
    vmin = med - sigma*mad
    vmax = med + sigma*mad

    fig, ax = plt.subplots()

    im = ax.imshow(frame, origin="lower", vmin=vmin, vmax=vmax, cmap="grey")
    cb = plt.colorbar(im)

    if title is not None:
        ax.set_title(title)

    # plt.show(block=False)
    return ax