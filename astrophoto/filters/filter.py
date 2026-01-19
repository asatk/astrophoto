import abc



class Filter(metaclass=abc.ABCMeta):

    def __init__(self, key: str, val: str):
        self._key = key
        self._val = val

    @property
    def key(self):
        return self._key

    @property
    def val(self):
        return self._val

    @abc.abstractmethod
    def apply(self, obj):
        ...

class FilterEQ(Filter):

    def apply(self, obj):
        return obj[self._key] == self._val



class FilterNE(Filter):

    def apply(self, obj):
        return obj[self._key] != self._val



class FilterGT(Filter):

    def apply(self, obj):
        return obj[self._key] > self._val



class FilterLT(Filter):

    def apply(self, obj):
        return obj[self._key] < self._val



class FilterGE(Filter):

    def apply(self, obj):
        return obj[self._key] >= self._val



class FilterLE(Filter):

    def apply(self, obj):
        return obj[self._key] <= self._val