from astrophoto.objects import Frame
from astrophoto.combiners import MedianCombiner, Combiner
from astrophoto.objects import StackableCube
from astrophoto.routines import Routine



class StackingRoutine(Routine):
    """
    Generic routine for stacking cubes. Co-adds frames in a StackableCube
    according to the provided Combiner.
    """

    def __init__(self, combiner: Combiner=None):
        super(StackingRoutine, self).__init__()

        self._combiner = MedianCombiner() if combiner is None else combiner

    @property
    def combiner(self):
        return self._combiner

    def run(self, cube: StackableCube):
        frame = cube.stack(self.combiner)
        return frame



class CalculateBias(StackingRoutine):
    """
    Stacks bias frames.
    """

    def __init__(self,
                 combiner: Combiner=None):
        super(CalculateBias, self).__init__(combiner)

    def run(self, cube: StackableCube):
        frame = cube.stack(self.combiner)
        return frame


# TODO dark current
class CalculateDark(StackingRoutine):

    def __init__(self,
                 bias: Frame|str=None,
                 combiner: Combiner=None):
        super(CalculateDark, self).__init__(combiner)

        if isinstance(bias, Frame):
            self._bias = bias
        elif isinstance(bias, str):
            self._bias = Frame.load_fits(bias)
        elif bias is None:
            self._bias = 0
        else:
            raise TypeError(f"`bias` must be either a Frame, str, or None.")

    def run(self, cube: StackableCube, is_dark_current: bool=False):
        cube -= self._bias
        frame = cube.stack(self.combiner)
        return frame



class CalculateFlat(StackingRoutine):

    def __init__(self,
                 bias: Frame=None,
                 dark: Frame=None,
                 combiner: Combiner=None):
        super(CalculateFlat, self).__init__(combiner)
        self._bias = bias
        self._dark = dark

    def run(self, cube: StackableCube):
        cube -= self._bias
        cube -= self._dark
        frame = cube.stack(self.combiner)
        return frame