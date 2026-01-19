import abc
from typing import Callable

import numpy as np

class Combiner(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def apply(self, frames: np.ndarray|list[np.ndarray]) -> np.ndarray:
        ...



class MeanCombiner(Combiner):

    def apply(self, frames: np.ndarray|list[np.ndarray]) -> np.ndarray:
        return np.mean(frames, axis=0)



class MedianCombiner(Combiner):

    def apply(self, frames: np.ndarray|list[np.ndarray]) -> np.ndarray:
        return np.median(frames, axis=0)



class SumCombiner(Combiner):

    def apply(self, frames: np.ndarray|list[np.ndarray]) -> np.ndarray:
        return np.sum(frames, axis=0)



class CustomCombiner(Combiner):

    def __init__(self, f_combine: Callable[[np.ndarray|list[np.ndarray]], np.ndarray]):
        self._f_combine = f_combine

    @property
    def f_combine(self):
        return self._f_combine

    @f_combine.setter
    def f_combine(self, f_combine: Callable[[np.ndarray|list[np.ndarray]], np.ndarray]):
        self._f_combine = f_combine

    def apply(self, frames: np.ndarray|list[np.ndarray]) -> np.ndarray:
        return self._f_combine(frames)