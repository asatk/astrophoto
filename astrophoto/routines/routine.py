from abc import abstractmethod, ABCMeta

from astrophoto.core import Frame
from astrophoto.core.cube import Cube



class Routine(metaclass=ABCMeta):
    """
    Base class for an arbitrary data-reduction routine on Cubes and Frames.
    """

    @abstractmethod
    def run(self, *args, **kwargs) -> Frame|Cube:
        ...

    def __call__(self, *args, **kwargs) -> Frame|Cube:
        return self.run(*args, **kwargs)




# TODO do trim
class TrimRoutine(Routine):

    def __init__(self):
        ...

    def run(self, *args, **kwargs):
        ...