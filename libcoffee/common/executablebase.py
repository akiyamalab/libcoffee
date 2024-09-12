from abc import ABC, abstractmethod
from pathlib import Path


class ExecutableBase(ABC):

    def __init__(self, exec: Path, verbose: bool = False):
        self._verbose = verbose
        self._exec = exec
        self.__done = False

    @abstractmethod
    def _run(self, file: Path) -> None: ...

    def run(self, file: Path) -> "ExecutableBase":
        self._run(file)
        self.__done = True
        return self

    @property
    def done(self) -> bool:
        return self.__done
