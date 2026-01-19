import numpy as np
import multiprocessing as mp

from astrophoto.combiners import Combiner
from astrophoto.filters.filter import Filter
from astrophoto.core.frame import Frame



class Cube:

    def __init__(self, frame_list: list[Frame]=None):
        self._frame_list = frame_list if frame_list is not None else []



    @property
    def frame_list(self) -> list[Frame]:
        return self._frame_list

    @frame_list.setter
    def frame_list(self, frame_list):
        self._frame_list = frame_list


    # def __add__(self, other):
    #     return Cube([frame.data + other for frame in self.frame_list])
    #
    # def __sub__(self, other):
    #     return Cube([frame.data - other for frame in self.frame_list])
    #
    # def __mul__(self, other):
    #     return Cube([frame.data * other for frame in self.frame_list])
    #
    # def __truediv__(self, other):
    #     return Cube([frame.data / other for frame in self.frame_list])

    def __iter__(self):
        return iter(self._frame_list)

    def __len__(self):
        return len(self._frame_list)

    def __getitem__(self, item):
        return self._frame_list[item]

    def __setitem__(self, key, value):
        self._frame_list[key] = value

    def __array__(self):
        return self._frame_list



    def append(self, frame: Frame):
        self._frame_list.append(frame)



    @classmethod
    def load_fits(cls, filename_list: list[str], filter_list: list[Filter]=None, num_workers: int=1):
        filter_list = filter_list if filter_list is not None else []
        if num_workers > 1:
            with mp.Pool(num_workers) as pool:
                frame_list = pool.starmap(Frame.load_fits, zip(filename_list, filter_list))
        else:
            frame_list = [Frame.load_fits(filename, filter_list) for filename in filename_list]

        frame_list = [frame for frame in frame_list if frame is not None]
        return Cube(frame_list)



class StackableCube(Cube):

    def __init__(self, frame_list: list[Frame]|np.ndarray=None):

        super(StackableCube, self).__init__(frame_list)

        if not isinstance(frame_list, np.ndarray):
            if not self._validate_shapes():
                raise ValueError("StackableCube must have homogeneous shapes in all frame")


        self._data_cube = np.array(self._frame_list)
        self._shape = self._data_cube.shape



    def _validate_shapes(self):

        shapes = [frame.shape for frame in self._frame_list]
        shapes_unique = np.unique(shapes, axis=0)
        if len(shapes_unique) > 1:
            return False

        return True



    @property
    def data_cube(self) -> np.ndarray:
        return self._data_cube

    @data_cube.setter
    def data_cube(self, data_cube):
        self._data_cube = data_cube

    @property
    def shape(self):
        return self._shape



    def __add__(self, other):
        return StackableCube(self.data_cube + other)

    def __sub__(self, other):
        return StackableCube(self.data_cube - other)

    def __mul__(self, other):
        return StackableCube(self.data_cube * other)

    def __truediv__(self, other):
        return StackableCube(self.data_cube / other)

    def __radd__(self, other):
        return StackableCube(other + self.data_cube)

    def __rsub__(self, other):
        return StackableCube(other - self.data_cube)

    def __rmul__(self, other):
        return StackableCube(other * self.data_cube)

    def __rtruediv__(self, other):
        return StackableCube(other / self.data_cube)

    def __iadd__(self, other):
        if isinstance(other, StackableCube):
            self.data_cube += other.data_cube
        else:
            self.data_cube += other

    def __isub__(self, other):
        if isinstance(other, StackableCube):
            self.data_cube -= other.data_cube
        else:
            self.data_cube -= other

    def __imul__(self, other):
        if isinstance(other, StackableCube):
            self.data_cube *= other.data_cube
        else:
            self.data_cube *= other

    def __itruediv__(self, other):
        if isinstance(other, StackableCube):
            self.data_cube /= other.data_cube
        else:
            self.data_cube /= other




    # def __add__(self, other):
    #     if isinstance(other, StackableCube):
    #         if other.shape == self._shape:
    #             return StackableCube(self.data_cube + other.data_cube)
    #         else:
    #             raise ValueError("Frame must have same shape as each StackableCube frame")
    #     elif isinstance(other, Frame):
    #         if other.shape == self._shape[1:]:
    #             return StackableCube(self.data_cube + other.data)
    #         else:
    #             raise ValueError("Frame must have same shape as each StackableCube frame")
    #     elif isinstance(other, np.ndarray):
    #         shape = other.shape
    #         if len(shape) == 2 and shape == self._shape[1:]:
    #             return StackableCube(self.data_cube + other)
    #         elif len(shape) == 3 and (shape == self._shape or (shape[0] == 1 and shape[1:] == self._shape[1:])):
    #             return StackableCube(self.data_cube + other)
    #         else:
    #             raise ValueError("StackableCube arithmetic only performed on data with broadcastable types and shapes.")
    #     elif isinstance(other, (int, float)):
    #         return StackableCube(self.data_cube + other)
    #     else:
    #         raise TypeError("Cannot add %s to %s" % (type(other), StackableCube))






    def stack(self, combiner: Combiner, header: dict=None) -> Frame:

        if len(self._frame_list) == 0:
            raise RuntimeError("No frames available to combine.")

        if len(self._frame_list) == 1:
            raise RuntimeError("Only one frame available to combine -- nothing to do.")

        data_stacked = combiner.apply(self._frame_list)
        frame = Frame(data=data_stacked, header=header)
        return frame