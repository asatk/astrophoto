import argparse as ap
from astropy.io import ascii, fits
import aplpy
from matplotlib import pyplot as plt
import numpy as np
import os
import re

#parser = ap.ArgumentParser()
#parser.add_argument()

def prompt(msg: str, default: str):
    resp = input(f"{'='*80}\n{msg} [{default}] ")
    if resp == "":
        resp = str(default)
    return resp

root_dir = os.path.abspath(prompt("Project's root directory", "."))


# define vars for useful data locations
root_dir = "data/M82/"
# MANDATORY!!! directory structure for data set
bias_dir = root_dir + "/bias/"
dark_dir = root_dir + "/dark/"
flat_dir = root_dir + "/flat/"
data_dir = root_dir + "/raw/"

calib_dir = root_dir + "/calib/"
if not os.path.isdir(calib_dir):
    os.mkdir(calib_dir)

output_dir = root_dir + "/reduced/"
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)




def get_data_cube(d: str, p: str=r"", **kwargs):
    """
    Retrieve FITS image data from a directory `d` with file names matching a
    string regex pattern `p` with a header matching specific values in kwargs.
    """

    file_list = []
    data_list = []

    for f in os.listdir(d):
        # test filename matching regex pattern
        if re.search(p, f) is None:
            continue

        is_match = True
        hdr = fits.getheader(f)
        # test fits header matches header conditions
        for k, v in kwargs.items():
            if hdr.get(k, None) != v:
                is_match = False

        if not is_match:
            continue

        # append file's name to list of file names
        file_list.append(f)

        # get raw data from fits file
        data - fits.getdata(f)

        # append file's data to data cube
        data_list.append(data)

    # make numpy data cube from list of numpy data
    data_cube = np.array(data_list)

    return data_cube, file_list



def trim(data, xlo: int=None, xhi: int=None, ylo: int=None, yhi: int=None):
    """
    Trims the last two dimensions of >2D numpy array data. The x- and y-
    coordinates are those are read in DS9 which has x/y axes flipped relative
    to Python.
    """
    return data[...,ylo:yhi,xlo:xhi]



def qplot(data, title: str=None):
    """
    'Quick' plotting function. Does not have many features but is meant to be
    used in a pinch for debugging.
    """
    med = np.nanmedian(data)
    mad = np.nanmean(np.abs(data - med))
    print(f"Median: {med:.0f}\nMAD: {mad:.3f}")

    im = plt.imshow(data, vmin=med-mad, vmax=med+mad, cmap="grey_r", origin="lower")
    cb = plt.colorbar(im)

    if title is not None:
        plt.title(title)

    plt.show()



# load all calib data
biases, bias_files = get_data_cube(bias_dir)
darks, darks_files = get_data_cube(dark_dir)
flats_B, flats_B_files = get_data_cube(flat_dir, INSFILTE="B")
flats_V, flats_V_files = get_data_cube(flat_dir, INSFILTE="V")
flats_R, flats_R_files = get_data_cube(flat_dir, INSFILTE="R")

# trim all calib data
biases = trim(biases)
darks = trim(darks)
flats_B = trim(flats_B)
flats_V = trim(flats_V)
flats_R = trim(flats_R)


# combine calib files for master bias and master dark frames
bias = np.median(biases, axis=0)
dark = np.median(darks - bias, axis=0)
texp_dark = float(fits.getheader(darks_files[0])["EXPTIME"])
print(f"Median Dark Current: {np.median(dark) / texp_dark}")

texp

# write master frames to disk
fits.writeto(calib_dir + "bias.fit", bias, fits.getheader(bias_files[0]), overwrite=True)
fits.writeto(calib_dir + "dark.fit", dark, fits.getheader(dark_files[0]), overwrite=True)

# plot master frames
qplot(bias, "Master Bias")
qplot(dark, "Master Bias-Subtracted Dark Counts")


# calculate master normalized flat field response frames for each band
flat_B = np.median(flats_B - bias, axis=0)
flat_B_norm = flat_B / np.mean(flat_B)
flat_V = np.median(flats_V - dark - bias, axis=0)
flat_V_norm = flat_V / np.mean(flat_V)
flat_R = np.median(flats_R - dark - bias, axis=0)
flat_R_norm = flat_R / np.mean(flat_R)

# write master flats to disk
fits.writeto(flat_dir + "flat_B.fit", flat_B, fits.getheader(flat_dir + "flat_B01.fit"), overwrite=True)
fits.writeto(flat_dir + "flat_V.fit", flat_V, fits.getheader(flat_dir + "flat_V01.fit"), overwrite=True)
fits.writeto(flat_dir + "flat_R.fit", flat_R, fits.getheader(flat_dir + "flat_R01.fit"), overwrite=True)

# plot master flats
qplot(flat_B_norm, "Normalized B Flat Field")
qplot(flat_V_norm, "Normalized V Flat Field")
qplot(flat_R_norm, "Normalized R Flat Field")


# In[10]:


# shift and scale data WITHOUT combining -- frames are returned in the same order but aligned
def coadd(data_cube, xsh, ysh, skyval=None, starval=None):

    # reference image index -- doesn't matter since the subtractions/division are all relative
    ind_ref = 0

    # ds9 x coordinate offsets
    xoff = np.subtract.outer(xsh, xsh)[ind_ref]
    xpadlo = np.max(np.clip(xoff, a_min=0, a_max=None))
    xpadhi = np.min(np.clip(xoff, a_min=None, a_max=0))

    # ds9 y coordinate offsets
    yoff = np.subtract.outer(ysh, ysh)[ind_ref]
    ypadlo = np.max(np.clip(yoff, a_min=0, a_max=None))
    ypadhi = np.min(np.clip(yoff, a_min=None, a_max=0))

    # no padding for time/frame axis
    # python and ds9 x/y flipped
    pad = ((0, 0), (ypadlo, abs(ypadhi)), (xpadlo, abs(xpadhi)))

    # pad all frames in data cube
    data_pad = np.pad(data_cube, pad, constant_values=np.nan)

    # shift each frame in cube by each's offset
    data_roll = np.array([
        np.roll(data_pad[i], (yoff[i], xoff[i]), axis=(0,1))
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


# In[11]:


# reduce raw data frames for a given band
def reduce_data(band: str, flat: np.ndarray, shift_file: str):

    # only B, V, and R images for this dataset
    if band not in ["B", "V", "R"]:
        print("`band` must be one of ['B', 'V', 'R']")
        return None

    # clean raw science files (trim, bias, dark, flat)
    raw, fname = get_data_cube(data_dir + f"M57_{band}??.fit")
    raw = trim(raw)
    data = (raw - dark - bias) / flat

    # read shift data for coadd
    shift = ascii.read(shift_file)
    shift["path"] = [os.path.abspath(data_dir + f) for f in shift["name"]]
    ind_files = [list(shift["path"]).index(f) for f in fname if f in shift["path"]]
    shift = shift[ind_files]

    # align and coadd cleaned data files
    data_coadd = coadd(data, shift["xsh"], shift["ysh"], None)
    obj = np.median(data_coadd, axis=0)

    # write cleaned data to disk
    for i in ind_files:
        fits.writeto(output_dir + f"M57_{band}{i+1:02d}_clean.fit", data_coadd[i], fits.getheader(data_dir + f"M57_{band}{i+1:02d}.fit"), overwrite=True)

    # return combined frame in the given band
    return obj


# In[12]:


# reduce raw data in each band
obj_B = reduce_data("B", flat_B_norm, data_dir + "shift_B.txt")
obj_V = reduce_data("V", flat_V_norm, data_dir + "shift_V.txt")
obj_R = reduce_data("R", flat_R_norm, data_dir + "shift_R.txt")

# find minimum axis length for each dim across all reduced band frames
dims = np.array([x.shape for x in [obj_B, obj_V, obj_R]])
xmin = np.min(dims[:,1])
ymin = np.min(dims[:,0])

# trim frames to be same size
obj_B = obj_B[:ymin, :xmin]
obj_V = obj_V[:ymin, :xmin]
obj_R = obj_R[:ymin, :xmin]

# read shift and scale values for each reduced band image
shift = ascii.read(output_dir + "shift.txt")
xsh = shift["xsh"]
ysh = shift["ysh"]
skyval = shift["skyval"]
starval = shift["starval"]
obj_cube = np.array([obj_B, obj_V, obj_R])

# shift and scale each reduced band image
obj_coadd = coadd(obj_cube, xsh, ysh, skyval, starval)
obj_B, obj_V, obj_R = obj_coadd
obj_intensity = np.median(obj_coadd, axis=0)

# write reduced band frame to disk
fits.writeto(output_dir + f"M57_B.fit", obj_B, fits.getheader(data_dir + f"M57_B01.fit"), overwrite=True)
fits.writeto(output_dir + f"M57_V.fit", obj_V, fits.getheader(data_dir + f"M57_V01.fit"), overwrite=True)
fits.writeto(output_dir + f"M57_R.fit", obj_R, fits.getheader(data_dir + f"M57_R01.fit"), overwrite=True)

# save intensity map
fits.writeto(output_dir + "M57_grey.fit", obj_intensity, header=fits.getheader(data_dir + "M57_B01.fit"), overwrite=True)


# In[13]:


im = plt.imshow(obj_B, cmap="grey_r", vmin=2000, vmax=8000, origin="lower")
plt.title("Reduced B Band Image")
plt.show()
im = plt.imshow(obj_V, cmap="grey_r", vmin=3000, vmax=9000, origin="lower")
plt.title("Reduced V Band Image")
plt.show()
im = plt.imshow(obj_R, cmap="grey_r", vmin=2000, vmax=8000, origin="lower")
plt.title("Reduced R Band Image")
plt.show()
plt.imshow(obj_intensity, cmap="grey_r", vmin=0, vmax=1000, origin="lower")
plt.title("Intensity Map of Object")
plt.show()


# In[14]:


# force FITS header with valid WCS params for object onto shifted and scaled files
#WCS_header = fits.getheader(data_dir + "M57_DSS.fits")
#del(WCS_header["SKEW"])
wcs_B = fits.getheader(output_dir + "M57_B_WCS_header.fits")
wcs_V = fits.getheader(output_dir + "M57_V_WCS_header.fits")
wcs_R = fits.getheader(output_dir + "M57_R_WCS_header.fits")

# write reduced band frame to disk -- overwrite with old header + WCS
fits.writeto(output_dir + f"M57_B.fit", obj_B, wcs_B, overwrite=True)
fits.writeto(output_dir + f"M57_V.fit", obj_V, wcs_V, overwrite=True)
fits.writeto(output_dir + f"M57_R.fit", obj_R, wcs_R, overwrite=True)


# In[15]:


# aplpy make rgb image from shifted and scaled frames with valid WCS coords
obj_fnames = [output_dir + "M57_R.fit", output_dir + "M57_V.fit", output_dir + "M57_B.fit"]
aplpy.make_rgb_cube(obj_fnames, output_dir + "M57.fits") # NOTE! extension must be FITS (S!) because multiple channels (?)
aplpy.make_rgb_image(output_dir + "M57.fits", output_dir + "M57.png")


# In[16]:


# plot using aplpy
f = aplpy.FITSFigure(output_dir + "M57.png")
f.show_rgb()


# In[ ]:




