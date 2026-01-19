import numpy as np
from astropy.io import fits

from astrophoto.filters.filter import Filter


class Frame:

    def __init__(self, data: np.ndarray=None, header: dict=None):
        self._data = data
        self._header = header



    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        self._data = value

    @property
    def header(self):
        return self._header

    @header.setter
    def header(self, value):
        self._header = value

    def shape(self):
        return self._data.shape

    def dtype(self):
        return self._data.dtype



    def __array__(self):
        return self._data

    def __add__(self, other):
        return Frame(self.data + other)

    def __sub__(self, other):
        return Frame(self.data - other)

    def __mul__(self, other):
        return Frame(self.data * other)

    def __truediv__(self, other):
        return Frame(self.data / other)

    def __radd__(self, other):
        return Frame(other + self.data)

    def __rsub__(self, other):
        return Frame(other - self.data)

    def __rmul__(self, other):
        return Frame(other * self.data)

    def __rtruediv__(self, other):
        return Frame(other / self.data)

    def __iadd__(self, other):
        if isinstance(other, Frame):
            self.data += other.data
        else:
            self.data += other

    def __isub__(self, other):
        if isinstance(other, Frame):
            self.data -= other.data
        else:
            self.data -= other

    def __imul__(self, other):
        if isinstance(other, Frame):
            self.data *= other.data
        else:
            self.data *= other

    def __itruediv__(self, other):
        if isinstance(other, Frame):
            self.data /= other.data
        else:
            self.data /= other

    def __getitem__(self, item):
        if isinstance(item, slice) or isinstance(item, int):
            return self._data[item]
        elif isinstance(item, str):
            return self._header[item]
        else:
            raise TypeError("Key must be slice or int to access the data or str to access the header.")

    def __setitem__(self, key, value):
        if isinstance(key, slice) or isinstance(key, int):
            self._data[key] = value
        elif isinstance(key, str):
            self._header[key] = value
        else:
            raise TypeError("Key must be slice or int to access the data or str to access the header.")

    def __del__(self):
        del self._data
        del self._header



    def save_fits(self, filename: str, overwrite: bool=True):
        fits.writeto(filename, self.data, self.header, overwrite=overwrite)



    @classmethod
    def load_fits(cls, filename: str, filter_list: list[Filter]=None):
        filter_list = [] if filter_list is None else filter_list

        header = fits.getheader(filename)
        is_valid = True
        for fltr in filter_list:
            is_valid = fltr.apply(header)
        if is_valid:
            data = fits.getdata(filename)
            return Frame(data=data, header=header)
        else:
            del header
            return None
