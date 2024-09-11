from abc import ABC, abstractmethod
from pathlib import Path


class ExecutableBase(ABC):

    def __init__(self, exec: Path, verbose: bool = False):
        self._verbose = verbose
        self._exec = exec

    @abstractmethod
    def run(self) -> "ExecutableBase": ...
