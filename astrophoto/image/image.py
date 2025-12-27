import photutils.aperture
from astropy.io import ascii, fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats, SigmaClip
import glob
from matplotlib import pyplot as plt
import numpy as np
import os

from matplotlib.patches import Circle
from photutils.aperture import CircularAperture
from photutils.background import MedianBackground, Background2D, SExtractorBackground
from photutils.centroids import centroid_2dg, centroid_sources
from photutils.detection import DAOStarFinder
from photutils.psf import CircularGaussianPSF, PSFPhotometry, CircularGaussianPRF
from photutils.segmentation import detect_threshold
from scipy.ndimage import shift

from astrophoto.io import prompt, prompt_choices, printf, COMBINE_MEDIAN, COMBINE_MEAN, println
from astrophoto.viz import plot_frame

# shift and scale data WITHOUT combining -- frames are returned in the same order but aligned
def coadd(data_cube, xsh, ysh, skyval=None, starval=None):

    # reference image index -- doesn't matter since the subtractions/division are all relative
    ind_ref = 0

    # ds9 x coordinate offsets
    xoff = np.subtract.outer(xsh, xsh)[ind_ref]
    # xoff = np.round(xoff).astype(np.int64)
    xpadlo = np.max(np.clip(xoff, a_min=0, a_max=None))
    xpadhi = np.min(np.clip(xoff, a_min=None, a_max=0))

    # ds9 y coordinate offsets
    yoff = np.subtract.outer(ysh, ysh)[ind_ref]
    # yoff = np.round(yoff).astype(np.int64)
    ypadlo = np.max(np.clip(yoff, a_min=0, a_max=None))
    ypadhi = np.min(np.clip(yoff, a_min=None, a_max=0))

    # no padding for time/frame axis
    # python and ds9 x/y flipped
    pad_arr = np.array([(0, 0), (xpadlo, abs(xpadhi)), (ypadlo, abs(ypadhi))])
    pad = np.astype(pad_arr, np.int64)

    # pad all frames in data cube
    data_pad = np.pad(data_cube, pad, constant_values=np.nan)

    # shift each frame in cube by each's offset
    # data_roll = np.array([
    #     np.roll(data_pad[i], (yoff[i], xoff[i]), axis=(0,1))
    #     for i in range(data_pad.shape[0])
    # ])

    data_roll = np.array([
        shift(data_pad[i], (xoff[i], yoff[i]), mode="constant", cval=np.nan)
        for i in range(data_pad.shape[0])
    ])

    # trim padding/nan/bad data out of shifted frames
    xtrimlo = None if xpadlo == 0 else 2*xpadlo
    xtrimhi = None if xpadhi == 0 else 2*xpadhi
    ytrimlo = None if ypadlo == 0 else 2*ypadlo
    ytrimhi = None if ypadhi == 0 else 2*ypadhi

    data_trim = data_roll[:,ytrimlo:ytrimhi,xtrimlo:xtrimhi]

    # shift and scale counts in each frame by bkgd sky value and star counts
    if skyval is not None and starval is not None:
        data_trim -= skyval[:,None,None]
        starval -= skyval

        starval_norm = np.divide.outer(starval, starval)[ind_ref]
        data_trim *= starval_norm[:,None,None]

    return data_trim



def calib_image(root_dir):

    file_pattern = root_dir + "/" + prompt("File pattern (glob pattern)", "fpattern")
    all_file_list = glob.glob(file_pattern)

    band = prompt("Instrument Filter (str)", "flat_band", str)

    data_list = []
    texp_list = []
    file_list = []

    for f in all_file_list:

        hdr = fits.getheader(f)
        band_f = hdr["FILTER"]

        if band_f == band:
            data_list.append(fits.getdata(f).astype(np.float64))
            texp_list.append(hdr["EXPTIME"])
            file_list.append(f)

    if len(data_list) == 0:
        printf("No files loaded.")
        exit(20)

    imgs = np.array(data_list)
    texps = np.array(texp_list)

    dark_mode = prompt("Subtract none (0), dark frame (1), or dark current (2)?", "dark_mode", int)

    if dark_mode == 2:
        bias_fname = prompt(
            "Master bias file (path) ", "bias")
        bias_file = root_dir + "/" + bias_fname
        bias = fits.getdata(bias_file).astype(np.float64)

        imgs -= bias

        dark_fname = prompt("Master dark current file (path) ", "dark")
        dark_file = root_dir + "/" + dark_fname
        dark = fits.getdata(dark_file).astype(np.float64)

        darks = np.array([t * dark for t in texps])
        imgs -= darks

    elif dark_mode == 1:
        dark_fname = prompt("Master dark file (path) ", "dark")
        dark_file = root_dir + "/" + dark_fname
        dark = fits.getdata(dark_file).astype(np.float64)

        imgs -= dark

    # add check for same band
    # add info to calib file data
    flat_mode = prompt("Divide by none (0) or master flat (1)?", "flat_mode", int)

    if flat_mode:
        flat_fname = prompt(f"Master {band} flat (path) ", "flat")
        if flat_fname != "":
            flat_file = root_dir + "/" + flat_fname
            flat = fits.getdata(flat_file).astype(np.float64)
            imgs /= flat

    sigma_clip = SigmaClip(sigma=3.0, maxiters=10)

    final_frames = []

    for img in imgs:
        # estimate background and background RMS
        _, bmedian, bstd = sigma_clipped_stats(img, sigma=3.0)

        # determine detection threshold for sources
        threshold = detect_threshold(img, nsigma=2.0, sigma_clip=sigma_clip)

        #
        is_good = False

        while is_good:

            fwhm = prompt("FWHM of sources (pixels, float): ", "fwhm", float)
            roundlo = prompt("Minimum roundness for sources ([-1.0..1.0]): ", "roundlo", float)
            roundhi = prompt("Maximum roundness for sources ([-1.0..1.0]): ", "roundhi", float)

            # initialize instance of star finder object
            finder = DAOStarFinder(threshold=threshold, fwhm=fwhm, roundlo=roundlo, roundhi=roundhi)\

            # identify sources in image
            srcs = finder.find_stars(img - bmedian)

            if len(srcs) != 0:

                xcen = srcs["xcentroid"]
                ycen = srcs["ycentroid"]
                cens = np.array([xcen, ycen]).T

                ax = plot_frame(img, sigma=5.0)
                CircularAperture(cens, r=50.0).plot(ax=ax)

            else:
                printf("No sources found with the provided finder parameters.")

            search_option = int(input("Exclude image (0), continue to next frame (1), or refine source search (2)?"))

            if search_option == 0:
                break
            elif search_option == 1:
                final_frames.append(img)
                break





    centroids_in_fname = prompt("Use centroids file (path)? ", "centroids_in")
    if centroids_in_fname != "":
        centroids_in_file = root_dir + "/" + centroids_in_fname
        centroids_table = ascii.read(centroids_in_file, guess="csv")

        # read shift data for coadd
        centroids_table["path"] = [os.path.abspath(f) for f in centroids_table["file"]]
        ind_files = [list(centroids_table["path"]).index(os.path.abspath(f))
                     for f in file_list
                     if os.path.abspath(f) in centroids_table["path"]]
        centroids_table = centroids_table[ind_files]

    else:

        println()
        printf("Roughly identify the same star in each frame for alignment")
        centroids = []
        for img_ in imgs:
            plot_frame(img_)
            xi = float(input("Star x (pixel): "))
            yi = float(input("Star y (pixel): "))
            x, y = centroid_sources(img_, xi, yi, box_size=101,
                                    centroid_func=centroid_2dg)
            print(f"Star centroid: {x}, {y}")
            centroids.append((x[0], y[0]))
            plt.close()

        centroids = np.array(centroids)
        centroids_table = Table({
            "file": file_list,
            "xsh": centroids[:, 0],
            "ysh": centroids[:, 1],
        })


        centroids_out_fname = prompt("Save centroids to file (path)? ", "centroids_out")
        if centroids_out_fname != "":
            centroids_out_file = root_dir + "/" + centroids_out_fname
            ascii.write(centroids_table, centroids_out_file, format="csv")

    imgs_align = coadd(imgs, centroids_table["xsh"], centroids_table["ysh"])

    combine_mode = prompt_choices(
        "How to combine frames (num)?",
        """
        (1) Median
        (2) Mean
        """, "flat_combine", int)

    if combine_mode == COMBINE_MEDIAN:
        img = np.median(imgs_align, axis=0)
    elif combine_mode == COMBINE_MEAN:
        img = np.mean(imgs_align, axis=0)
    else:
        printf("Invalid frame combination mode")
        exit(10)

    plot_frame(img)
    plt.waitforbuttonpress()

    out_fname = prompt(f"Output file for combined {band} image (path)?", "fout")
    if out_fname != "":
        out_file = root_dir + "/" + out_fname
        fits.writeto(out_file, img, overwrite=True)

    printf(f"Combined {band} image saved.")

    return img



def stack_image(root_dir):

    file_pattern = root_dir + "/" + prompt("File pattern (glob pattern)", "stack_fpattern")
    file_list = glob.glob(file_pattern)

    data_list = [fits.getdata(f).astype(np.float64) for f in file_list]

    if len(data_list) == 0:
        printf("No files loaded.")
        exit(20)

    # imgs = np.array(data_list)
    imgs = data_list

    centroids_in_fname = prompt("Use centroids file (path)? ", "stack_centroids_in")
    if centroids_in_fname != "":
        centroids_in_file = root_dir + "/" + centroids_in_fname
        centroids_table = ascii.read(centroids_in_file, guess="csv")

        # read shift data for coadd
        centroids_table["path"] = [os.path.abspath(f) for f in centroids_table["file"]]
        ind_files = [list(centroids_table["path"]).index(os.path.abspath(f))
                     for f in file_list
                     if os.path.abspath(f) in centroids_table["path"]]
        centroids_table = centroids_table[ind_files]

    else:

        println()
        printf("Subtract sky background and fit 2D Gaussian PSF to a star common to all frames.")
        centroids = []
        bkgds = []
        scales = []

        bkg_size = prompt("Bkg mesh size (pixel) (larger than star; small enough for variations): ", "mesh_size", int)
        threshold = prompt("PSF threshold (count): ", "threshold", float)

        sclip = SigmaClip(sigma=3.0, maxiters=10)
        bkg_estimator = SExtractorBackground()
        psf_model = CircularGaussianPRF()
        fit_shape = (5,5)
        finder = DAOStarFinder(threshold, 10.0, roundlo=0.3)

        plt.ion()

        for img in imgs:

            img_ = img.copy()

            # Calculate 2D background (correct for gradients possibly present)
            bkg = Background2D(img_, (bkg_size, bkg_size), filter_size=(3, 3),
                               sigma_clip=sclip, bkg_estimator=bkg_estimator)
            plot_frame(bkg.background, title="background pattern")

            img_ -= bkg.background

            ax_ = plot_frame(img_, title="background subtracted frame")

            plt.waitforbuttonpress()

            phot = psf_phot(img_)

            res = psf_phot.finder_results

            printf(f"Plotting identified sources ({len(res)}).")

            for src in res:
                circle = Circle((src["xcentroid"], src["ycentroid"]), 5*aperture,
                                edgecolor="red", facecolor="none", linewidth=2, alpha=0.5)
                ax_.add_patch(circle)


            plt.show()

            # xi = float(input("Star x (pixel): "))
            # yi = float(input("Star y (pixel): "))
            # x, y = centroid_sources(img_, xi, yi, box_size=101,
            #                         centroid_func=centroid_2dg)
            # printf(f"Star centroid: {x}, {y}")
            # centroids.append((x[0], y[0]))

            # aperture = CircularAperture([[x], [y]], r=fwhm)

            input("Refine source extraction? [press any key to refine]")

            plt.close()

        plt.ioff()

        centroids = np.array(centroids)
        bkgds = np.array(bkgds)
        centroids_table = Table({
            "file": file_list,
            "xsh": centroids[:, 0],
            "ysh": centroids[:, 1],
            "sky": bkgds,
            "scale": scales,
        })

        centroids_out_fname = prompt("Save centroids to file (path)? ", "centroids_out")
        if centroids_out_fname != "":
            centroids_out_file = root_dir + "/" + centroids_out_fname
            ascii.write(centroids_table, centroids_out_file, format="csv")

    # imgs_align = coadd(imgs, centroids_table["xsh"], centroids_table["ysh"],
    #                    centroids_table["sky"], centroids_table["scale"])

    combine_mode = prompt_choices(
        "How to combine frames (num)?",
        """
        (1) Median
        (2) Mean
        """, "flat_combine", int)

    if combine_mode == COMBINE_MEDIAN:
        img = np.median(imgs_align, axis=0)
    elif combine_mode == COMBINE_MEAN:
        img = np.mean(imgs_align, axis=0)
    else:
        printf("Invalid frame combination mode")
        exit(10)

    plot_frame(img)
    plt.waitforbuttonpress()

    out_fname = prompt(f"Output file for combined {band} image (path)?", "fout")
    if out_fname != "":
        out_file = root_dir + "/" + out_fname
        fits.writeto(out_file, img, overwrite=True)

    printf(f"Combined {band} image saved.")

    return img
