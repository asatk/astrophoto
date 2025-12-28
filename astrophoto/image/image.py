import astroalign
from astropy.io import ascii, fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats, SigmaClip
import glob
from matplotlib import pyplot as plt
import numpy as np
import os

from matplotlib.patches import Circle
from photutils.aperture import CircularAperture
from photutils.background import Background2D, SExtractorBackground
from photutils.detection import DAOStarFinder
from photutils.psf import PSFPhotometry, CircularGaussianPRF
from scipy.ndimage import shift

from astrophoto.io import prompt, prompt_choices, printf, COMBINE_MEDIAN, COMBINE_MEAN, println
from astrophoto.io.reduxio import update_default, get_default, saveandquit
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
        saveandquit(20)

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

        # todo check that exposure times are equal to dark exp time if in master dark metadata

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

    println()
    printf("Identifying sources for frame alignment.")

    # large-scale background subtraction
    sclip_sigma = 3.0
    sigma_clip = SigmaClip(sigma=sclip_sigma, maxiters=10)
    bkg_estimator = SExtractorBackground()
    bkg_size = (int(imgs[0].shape[0] / 4), int(imgs[0].shape[1] / 4))

    # finder parameters
    thresh_nsigma = get_default("thresh_nsigma")
    fwhm = get_default("fwhm")
    roundlim = get_default("roundlim")

    final_frames = []
    ref_frame = None

    for img in imgs:

        # subtract any large-scale background gradients to make source extraction easier

        bkg = Background2D(img, bkg_size, filter_size=(3, 3),
                                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

        # background-subtracted image
        img_bsub = img - bkg.background

        # estimate background RMS for threshold val
        _, _, bstd = sigma_clipped_stats(img_bsub, sigma=sclip_sigma)


        keep_searching = True
        while keep_searching:

            # initialize instance of star finder object
            finder = DAOStarFinder(threshold=thresh_nsigma * bstd, fwhm=fwhm, roundlo=-roundlim, roundhi=roundlim)

            # identify sources in image
            srcs = finder.find_stars(img_bsub)

            if srcs is not None and len(srcs) >= 3:

                printf(f"{len(srcs)} source(s) identified.")

                # plot each identified source as small circle
                xcen = srcs["xcentroid"]
                ycen = srcs["ycentroid"]
                cens = np.array([xcen, ycen]).T

                ax = plot_frame(img_bsub, sigma=5.0)
                CircularAperture(cens, r=50.0).plot(ax=ax, color="red")
                plt.show()

            else:
                printf("Insufficient number of sources (< 3) found with the provided finder parameters.")

            search_option_str = input(
                """Verdict on current frame's source extraction [0]:
                [0] Continue to next frame
                [1] Exclude current frame
                [2] Refine source extraction
                """)
            search_option_str = search_option_str.strip()
            if search_option_str == "":
                search_option = 1
            else:
                search_option = int(search_option_str)

            if search_option == 0:

                # set the reference frame
                if len(final_frames) == 0:
                    ref_frame = img
                # align subsequent frames to reference frames
                else:
                    img, _ = astroalign.register(img, ref_frame)

                final_frames.append(img)
                keep_searching = False

            elif search_option == 1:
                keep_searching = False

            elif search_option == 2:

                # update
                thresh_nsigma = prompt(f"Number of sigma (={bstd:.3f}) above background level for source detection",
                                       "thresh_nsigma", float)
                fwhm = prompt("FWHM of sources (pixels, float): ", "fwhm", float)
                roundlim = abs(
                    prompt("Upper bound on absolute roundness for sources ([0.0..1.0]): ", "roundlim", float))

                update_default("thresh_nsigma", thresh_nsigma)
                update_default("fwhm", fwhm)
                update_default("roundlim", roundlim)

    if len(final_frames) == 0:
        printf("No frames selected to combine.")
        saveandquit(30)

    printf(f"Combining {len(final_frames)} frames.")

    combine_mode = prompt_choices(
        "How to combine frames (num)?",
        """
        (1) Median
        (2) Mean
        """, "flat_combine", int)

    if combine_mode == COMBINE_MEDIAN:
        img = np.median(final_frames, axis=0)
    elif combine_mode == COMBINE_MEAN:
        img = np.mean(final_frames, axis=0)
    else:
        printf("Invalid frame combination mode")
        saveandquit(10)

    plot_frame(img, sigma=5.0)
    plt.show()

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
