import astroalign as aa
import skimage.transform
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

from astrophoto.core import combine
from astrophoto.io import prompt, prompt_choices, printf, COMBINE_MEDIAN, COMBINE_MEAN, println
from astrophoto.io.reduxio import update_default, get_default, save_and_quit, load_files
from astrophoto.viz import plot_frame



def calib_image(root_dir):

    file_options = prompt("Use a stored file list (0) or match files to a pattern (1)", "calib_file_opt", int)
    if file_options == 0:
        flist_fname_in = prompt("File containing list of files: ", "calib_flist_in", str)
        all_file_list = []
        with open(flist_fname_in, "r") as flist_in:
            temp_list = flist_in.readlines()

        all_file_list = [s.strip() for s in temp_list]

    elif file_options == 1:
        file_pattern = root_dir + "/" + prompt("File pattern (glob pattern)", "calib_fpattern")
        all_file_list = glob.glob(file_pattern)
    else:
        printf("Invalid selection.")
        save_and_quit(20)



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
        save_and_quit(20)

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
    sharplo = get_default("sharplo")

    final_frames = []
    shifts = []
    final_files = []
    img_ref = None

    check_quality = prompt("Check each frame for quality (0=No, 1=Yes)?", "calib_quality", int)

    # todo extract sources in each image in parallel
    for i, img in enumerate(imgs):

        printf(f"Aligning frame [{i+1}/{len(imgs)}]: {file_list[i]}")

        # subtract any large-scale background gradients to make source extraction easier

        bkg = Background2D(img, bkg_size, filter_size=(3, 3),
                           sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

        # background-subtracted image
        img_bsub = img - bkg.background

        # estimate background RMS for threshold val
        _, _, bstd = sigma_clipped_stats(img_bsub, sigma=sclip_sigma)

        if check_quality:

            while True:

                # initialize instance of star finder object
                finder = DAOStarFinder(threshold=thresh_nsigma * bstd, fwhm=fwhm,
                                       roundlo=-roundlim, roundhi=roundlim,
                                       sharplo=sharplo)

                # identify sources in image
                srcs = finder.find_stars(img_bsub)

                if srcs is not None and len(srcs) >= 3:

                    printf(f"{len(srcs)} sources identified.")

                    # plot each identified source as small circle
                    xcen = srcs["xcentroid"]
                    ycen = srcs["ycentroid"]
                    cens = np.array([xcen, ycen]).T

                    ax = plot_frame(img_bsub, sigma=5.0)
                    CircularAperture(cens, r=50.0).plot(ax=ax, color="red")
                    plt.show()

                else:
                    printf("Insufficient number of sources (< 3) found with the provided finder parameters.")

                search_option = prompt_choices(
                    "Verdict on current frame's source extraction",
                    """
                    (0) Continue to next frame
                    (1) Exclude current frame
                    (2) Refine source extraction
                    """, "src_ext_decision", int)

                if search_option == 0 and srcs is not None and len(srcs) >= 3:

                    # set the reference frame
                    if len(final_frames) == 0:
                        img_ref = img_bsub
                        cens_ref = np.array([srcs["xcentroid"], srcs["ycentroid"]]).T
                    # align subsequent frames to reference frames
                    else:

                        try:
                            xform, _ = aa.find_transform(cens, cens_ref, max_control_points=100)
                        except:
                            printf("No matches between transformed frame and reference frame.\n" + \
                                   "Something went wrong estimating the transform; most likely, more sources are needed.")
                            continue


                        img, _ = aa.apply_transform(xform, img, img_ref)

                        shift_xy = xform.params[1::-1, 2]
                        shift_xy = np.astype(np.sign(shift_xy) * np.ceil(np.abs(shift_xy)), int)
                        shifts.append(shift_xy)

                    final_frames.append(img)
                    final_files.append(file_list[i])

                    printf(f"Aligned frame included in final data cube [{len(final_frames)} frame(s)].")

                    break

                elif search_option == 1:
                    break

                elif search_option == 2:

                    # update
                    thresh_nsigma = prompt(f"Number of sigma (={bstd:.3f}) above background level for source detection",
                                           "thresh_nsigma", float)
                    fwhm = prompt("FWHM of sources (pixels, float): ", "fwhm", float)
                    roundlim = abs(
                        prompt("Upper bound on absolute roundness for sources (float): ", "roundlim", float))
                    sharplo = abs(
                        prompt("Lower bound on sharpness for sources (float): ", "sharplo", float))

                    update_default("thresh_nsigma", thresh_nsigma)
                    update_default("fwhm", fwhm)
                    update_default("roundlim", roundlim)

        else:
            if i == 0:
                img_ref = img_bsub
                shifts.append([0, 0])
            else:
                # xform, _ = aa.find_transform(img_bsub, img_ref)
                xform, _ = aa.find_transform(img_bsub, img_ref, detection_sigma=thresh_nsigma)
                img, _ = aa.apply_transform(xform, img, img_ref)
                shift_xy = xform.params[1::-1, 2]
                shift_xy = np.astype(np.sign(shift_xy) * np.ceil(np.abs(shift_xy)), int)
                shifts.append(shift_xy)

            final_frames.append(img)
            final_files.append(file_list[i])

    if len(final_frames) == 0:
        printf("No frames selected to combine.")
        save_and_quit(30)

    printf(f"Combining {len(final_frames)} frames.")

    combine_mode = prompt_choices(
        "How to combine frames (num)?",
        """
        (1) Median
        (2) Mean
        """, "calib_combine", int)

    if combine_mode == COMBINE_MEDIAN:
        img = np.median(final_frames, axis=0)
    elif combine_mode == COMBINE_MEAN:
        img = np.mean(final_frames, axis=0)
    else:
        printf("Invalid frame combination mode")
        save_and_quit(10)

    # trim combined frame based on shifts
    shifts = np.transpose(shifts).astype(int)
    shiftsx = shifts[0]
    shiftsy = shifts[1]

    xpadlo = np.max(np.clip(shiftsx, a_min=0, a_max=None))
    xpadhi = np.min(np.clip(shiftsx, a_min=None, a_max=0))

    ypadlo = np.max(np.clip(shiftsy, a_min=0, a_max=None))
    ypadhi = np.min(np.clip(shiftsy, a_min=None, a_max=0))

    # trim padding/nan/bad data out of shifted frames
    xtrimlo = None if xpadlo == 0 else xpadlo
    xtrimhi = None if xpadhi == 0 else xpadhi
    ytrimlo = None if ypadlo == 0 else ypadlo
    ytrimhi = None if ypadhi == 0 else ypadhi

    img = img[xtrimlo:xtrimhi, ytrimlo:ytrimhi]

    plot_frame(img, sigma=5.0)
    plt.show()

    flist_fname_out = prompt("Output file for list of files selected in the final stack: ", "calib_flist_out", str)
    if flist_fname_out != "":
        with open(flist_fname_out, "w") as flist_out:
            for f in final_files:
                flist_out.write(f + "\n")
        printf("File list saved.")

    out_fname = prompt(f"Output file for combined {band} image (path)?", "calib_fout")
    if out_fname != "":
        out_file = root_dir + "/" + out_fname
        fits.writeto(out_file, img, overwrite=True)
        printf(f"Combined {band} image saved.")

    return img



def stack_image(root_dir):

    imgs = load_files(root_dir)

    println()
    printf("Identifying sources for frame alignment.")

    # background subtraction
    sclip_sigma = 3.0
    sigma_clip = SigmaClip(sigma=sclip_sigma, maxiters=10)
    bkg_estimator = SExtractorBackground(sigma_clip=sigma_clip)

    # bkg_size = int(np.sqrt(np.prod(imgs[0].shape)) / 32)
    # update_default("stack_bkg_size", bkg_size)

    # finder parameters
    bkg_size = get_default("stack_bkg_size")
    thresh_nsigma = get_default("thresh_nsigma")
    fwhm = get_default("fwhm")
    roundlim = get_default("roundlim")
    sharplo = get_default("sharplo")

    final_frames = []
    shifts = []
    img_ref = None

    check_quality = prompt("Check each frame for quality (0=No, 1=Yes)?", "stack_quality", int)

    i = 0
    while i < len(imgs):

        img = imgs[i].copy()

        printf(f"Aligning frame [{i + 1}/{len(imgs)}]: {file_list[i]}")

        plot_frame(img, sigma=5.0)
        plt.show(block=False)

        has_mask = prompt("Mask a region (0=No, 1=Yes)?", "stack_use_mask", int)
        mask = np.zeros_like(img, dtype=int) == 1
        if has_mask:
            while has_mask:
                mask_xlo = prompt("Mask rectangle lower x coordinate (pixel):", "stack_mask_xlo", int)
                mask_ylo = prompt("Mask rectangle lower y coordinate (pixel):", "stack_mask_ylo", int)
                mask_xhi = prompt("Mask rectangle upper x coordinate (pixel):", "stack_mask_xhi", int)
                mask_yhi = prompt("Mask rectangle upper y coordinate (pixel):", "stack_mask_yhi", int)

                mask_ones = np.zeros_like(img, dtype=int)
                mask_ones[mask_xlo:mask_xhi,mask_ylo:mask_yhi] = 1
                mask |= (mask_ones == 1)

                has_mask = prompt("Mask a region (0=No, 1=Yes"")?", "stack_use_mask", int)

        # subtract any large-scale background gradients to make source extraction easier

        bkg = Background2D(img, (bkg_size, bkg_size), filter_size=(3, 3),
                           sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
                           mask=mask)

        # background-subtracted image
        img_bsub = img - bkg.background

        if check_quality:

            plot_frame(bkg.background, "background pattern")
            plot_frame(img_bsub, "background-subtracted frame", sigma=5.0)
            plt.show()

            background_option = prompt_choices(
                "Verdict on background estimate",
                """
                (0) Continue to alignment; background pattern sufficient
                (1) Re-estimate background pattern
                """, "background_option", int)

            if background_option == 1:
                bkg_size = prompt(f"Background mesh size (img shape: {img.shape})", "stack_bkg_size", int)
                continue
            elif background_option != 0:
                printf("Invalid selection")
                save_and_quit(40)

            # estimate background RMS for threshold val
            _, _, bstd = sigma_clipped_stats(img_bsub, sigma=sclip_sigma)

            while True:

                # initialize instance of star finder object
                finder = DAOStarFinder(threshold=thresh_nsigma * bstd, fwhm=fwhm,
                                       roundlo=-roundlim, roundhi=roundlim,
                                       sharplo=sharplo)

                # identify sources in image
                srcs = finder.find_stars(img_bsub)

                if srcs is not None and len(srcs) >= 3:

                    printf(f"{len(srcs)} sources identified.")

                    # plot each identified source as small circle
                    xcen = srcs["xcentroid"]
                    ycen = srcs["ycentroid"]
                    cens = np.array([xcen, ycen]).T

                    ax = plot_frame(img_bsub, sigma=5.0)
                    CircularAperture(cens, r=50.0).plot(ax=ax, color="red")
                    plt.show()

                else:
                    printf("Insufficient number of sources (< 3) found with the provided finder parameters.")

                search_option = prompt_choices(
                    "Verdict on current frame's source extraction",
                    """
                    (0) Continue to next frame
                    (1) Exclude current frame
                    (2) Refine source extraction
                    """, "src_ext_decision", int)

                if search_option == 0 and srcs is not None and len(srcs) >= 3:

                    # set the reference frame
                    if len(final_frames) == 0:
                        img_ref = img_bsub
                        cens_ref = np.array([srcs["xcentroid"], srcs["ycentroid"]]).T
                        flux_ref = np.array(srcs["flux"])

                        printf("Flux scaling: 1.0 (reference frame)")

                    # align subsequent frames to reference frames
                    else:

                        cens_bsub = np.array([srcs["xcentroid"], srcs["ycentroid"]]).T
                        flux_bsub = np.array(srcs["flux"])

                        xform, _ = aa.find_transform(cens_bsub, cens_ref, max_control_points=100)
                        cens_xform = skimage.transform.matrix_transform(cens_bsub, xform)

                        shift_xy = np.astype(xform.params[:2, 2], int)
                        shifts.append(shift_xy)

                        # half-pixel uncertainty should be good enough to confirm that two pairs of coordinates
                        # refer to the same soure
                        unc_px = 0.5
                        cens_xmatch = np.abs(np.subtract.outer(cens_xform[:,0], cens_ref[:,0])) < unc_px
                        cens_ymatch = np.abs(np.subtract.outer(cens_xform[:,1], cens_ref[:,1])) < unc_px
                        cens_match = cens_xmatch & cens_ymatch
                        nmatches = np.sum(cens_match)

                        if nmatches < 1:
                            printf("No matches between transformed frame and reference frame.\n" + \
                                   "Something went wrong estimating the transform; most likely, more sources are needed.")
                            continue

                        ind_bsub = np.nonzero(cens_match @ np.ones(len(cens_ref)))
                        ind_ref = np.nonzero(np.ones(len(cens_bsub)) @ cens_match)

                        flux_ratios = flux_bsub[ind_bsub] / flux_ref[ind_ref]

                        scaling = np.median(flux_ratios)
                        scaling_mad = np.mean(np.abs(flux_ratios - np.mean(scaling)))

                        printf(f"Flux scaling: {scaling} +/- {scaling_mad}")

                        img_bsub, _ = aa.apply_transform(xform, img_bsub, img_ref)
                        #todo re-run ngc2264 b/c had scaling backwards
                        img_bsub /= scaling

                    printf("Aligned frame included in final data cube.")
                    final_frames.append(img_bsub)

                    fname_old = file_list[i]
                    fname = fname_old[:fname_old.rfind(".")] + "_aligned.fit"

                    fname_output = input(f"Save aligned and background-subtracted band frame to (path) [{fname}]: ")
                    fname_output = fname_output.strip()
                    if fname_output != "":
                        fname = fname_output

                    fits.writeto(fname, img_bsub, overwrite=True)

                    break

                elif search_option == 1:
                    break

                elif search_option == 2:

                    # update
                    thresh_nsigma = prompt(f"Number of sigma (={bstd:.3f}) above background level for source detection",
                                           "thresh_nsigma", float)
                    fwhm = prompt("FWHM of sources (pixels, float): ", "fwhm", float)
                    roundlim = abs(
                        prompt("Upper bound on absolute roundness for sources (float): ", "roundlim", float))
                    sharplo = abs(
                        prompt("Lower bound on sharpness for sources (float): ", "sharplo", float))

                    update_default("thresh_nsigma", thresh_nsigma)
                    update_default("fwhm", fwhm)
                    update_default("roundlim", roundlim)

                else:
                    printf("Invalid selection")
                    save_and_quit(40)

        else:
            if i == 0:
                img_ref = img_bsub
                shifts.append([0, 0])

            else:
                xform, _ = aa.find_transform(img_bsub, img_ref)

                shift_xy = np.astype(xform.params[:2, 2], int)
                shifts.append(shift_xy)

                img_bsub, _ = aa.apply_transform(xform, img_bsub, img_ref)

            printf("Aligned frame included in final data cube.")
            final_frames.append(img_bsub)

            fname_old = file_list[i]
            fname = fname_old[:fname_old.rfind(".")] + "_aligned.fit"

            fname_output = input(f"Save aligned and background-subtracted band frame to (path) [{fname}]: ")
            fname_output = fname_output.strip()
            if fname_output != "":
                fname = fname_output

            fits.writeto(fname, img_bsub, overwrite=True)

        i += 1

    img = combine(final_frames)

    # trim combined frame based on shifts
    shifts = np.transpose(shifts).astype(int)
    shiftsx = shifts[0]
    shiftsy = shifts[1]

    xpadlo = np.max(np.clip(shiftsx, a_min=0, a_max=None))
    xpadhi = np.min(np.clip(shiftsx, a_min=None, a_max=0))

    ypadlo = np.max(np.clip(shiftsy, a_min=0, a_max=None))
    ypadhi = np.min(np.clip(shiftsy, a_min=None, a_max=0))

    # trim padding/nan/bad data out of shifted frames
    xtrimlo = None if xpadlo == 0 else xpadlo
    xtrimhi = None if xpadhi == 0 else xpadhi
    ytrimlo = None if ypadlo == 0 else ypadlo
    ytrimhi = None if ypadhi == 0 else ypadhi

    # crop still not right
    img = img[xtrimlo:xtrimhi, ytrimlo:ytrimhi]

    plot_frame(img, sigma=5.0)
    plt.show()

    out_fname = prompt(f"Output file for total intensity image (path)?", "align_fout")
    if out_fname != "":
        out_file = root_dir + "/" + out_fname
        fits.writeto(out_file, img, overwrite=True)
        printf(f"Total intensity image saved.")

    return img
