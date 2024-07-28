from abc import ABC, abstractmethod
from pathlib import Path


class ExecutableBase(ABC):

    def __init__(self, verbose=False):
        self.verbose = verbose

    @abstractmethod
    def run(self): ...

    @property
    @abstractmethod
    def result(self) -> Path: ...
