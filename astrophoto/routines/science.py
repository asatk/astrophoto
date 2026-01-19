import numpy as np
import astroalign as aa
from astropy.stats import SigmaClip, sigma_clipped_stats
import multiprocessing as mp
from photutils.background import SExtractorBackground, Background2D, BackgroundBase

from astrophoto.objects import Frame, Cube
from photutils.detection import DAOStarFinder

from astrophoto.routines.routine import Routine


class CalibrateScienceFrames(Routine):

    def __init__(self,
                 bias: Frame=None,
                 dark: Frame=None,
                 flat: Frame=None):
        self._bias = bias
        self._dark = dark
        self._flat = flat
        ...

    # TODO add arithmetic to Cube
    def run(self, cube: Cube):
        cube -= self._bias
        cube -= self._dark
        cube /= self._flat
        return cube



class BackgroundSubtraction(Routine):

    def __init__(self,
                 box_size: int | tuple[int] = 100,
                 sigma: float=3.0,
                 estimator: BackgroundBase=None):

        self._box_size = box_size
        self._sigma = sigma
        self._estimator = SExtractorBackground() if estimator is None else estimator



    @property
    def box_size(self):
        return self._box_size

    @box_size.setter
    def box_size(self, value):
        self._box_size = value

    @property
    def sigma(self):
        return self._sigma

    @sigma.setter
    def sigma(self, value):
        self._sigma = value

    @property
    def estimator(self):
        return self._estimator

    @estimator.setter
    def estimator(self, value):
        self._estimator = value



    def run(self, cube: Cube, num_workers: int=1, mask_cube: Cube=None):

        def _frame_bsub(frame_: Frame, mask_: Frame=None):
            bkg = Background2D(frame_,
                               box_size=self._box_size,
                               mask=mask_,
                               filter_size=(3, 3),
                               sigma_clip=SigmaClip(self._sigma),
                               bkg_estimator=self._estimator)
            frame_bsub_ = frame_ - bkg.background
            return frame_bsub_

        if mask_cube is None:
            mask_cube = [None] * len(cube)

        if num_workers > 1:
            with mp.Pool(num_workers) as pool:
                result = pool.starmap(_frame_bsub, zip(cube, mask_cube))
                cube_bsub = Cube(result)
        else:
            cube_bsub = Cube([_frame_bsub(frame, mask) for frame, mask in zip(cube, mask_cube)])

        return cube_bsub



class AlignScienceFrames(Routine):
    def __init__(self,
                 sigma: float=3.0,
                 thresh_sigma: float=10.0,
                 fwhm: float=15.0,
                 roundlim: float=0.7,
                 sharplo: float=0.25):
        self._sigma = sigma
        self._thresh_sigma = thresh_sigma
        self._fwhm = fwhm
        self._roundlim = roundlim
        self._sharplo = sharplo

    @property
    def sigma(self):
        return self._sigma

    @sigma.setter
    def sigma(self, value):
        self._sigma = value

    @property
    def thresh_sigma(self):
        return self._thresh_sigma

    @thresh_sigma.setter
    def thresh_sigma(self, value):
        self._thresh_sigma = value

    @property
    def fwhm(self):
        return self._fwhm

    @fwhm.setter
    def fwhm(self, value):
        self._fwhm = value

    @property
    def roundlim(self):
        return self._roundlim

    @roundlim.setter
    def roundlim(self, value):
        self._roundlim = value

    @property
    def sharplo(self):
        return self._sharplo

    @sharplo.setter
    def sharplo(self, value):
        self._sharplo = value



    def run(self, cube: Cube, num_workers: int=1):

        # estimate background RMS to determine threshold level for detection
        bstd_med = np.median([
            sigma_clipped_stats(img_, sigma=self._sigma)[2] for img_ in cube])

        # initialize instance of star finder object
        finder = DAOStarFinder(threshold=self._thresh_sigma * bstd_med,
                               fwhm=self._fwhm,
                               roundlo=-self._roundlim,
                               roundhi=self._roundlim,
                               sharplo=self._sharplo)

        def _find_srcs(frame: Frame):
            data = frame.data
            srcs = finder.find_stars(data)
            xcen = srcs["xcentroid"]
            ycen = srcs["ycentroid"]
            cens = np.array([xcen, ycen]).T
            return cens

        if num_workers > 1:
            with mp.Pool(num_workers) as pool:
                cens_all = pool.map(_find_srcs, cube)
        else:
            cens_all = [_find_srcs(frame) for frame in cube]

        nsrcs = [len(cens_) for cens_ in cens_all]

        ref_idx = np.argmax(nsrcs)
        cens_ref = cens_all[ref_idx]    # TODO option to use centers or just pure image align... or exclude DAOStarFinder at all...
        frame_ref = cube[ref_idx]

        def _find_xform(frame: Frame):
            data = frame.data
            data_ref = frame_ref.data
            xform, _ = aa.find_transform(data, data_ref, detection_sigma=self._thresh_sigma)
            data_xform, _ = aa.apply_transform(xform, data, data_ref)
            shift_xy = xform.params[1::-1, 2]
            shift_xy = np.astype(np.sign(shift_xy) * np.ceil(np.abs(shift_xy)), int)
            return data_xform, shift_xy

        if num_workers > 1:
            with mp.Pool(num_workers) as pool:
                result = pool.map(_find_xform, cube)
        else:
            result = [_find_xform(frame) for frame in cube]

        cube_xform = Cube([res_[0] for res_ in result])
        shifts = np.array([res_[1] for res_ in result])

        # TODO trim cube so it can be stacked

        return cube_xform

