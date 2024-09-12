from abc import ABC, abstractmethod
from pathlib import Path


class ExecutableBase(ABC):

    def __init__(self, exec: Path, verbose: bool = False):
        self._verbose = verbose
        self._exec = exec
        self.__done = False

    @abstractmethod
    def _run(self, *args, **kwargs) -> None: ...

    def run(self, *args, **kwargs) -> "ExecutableBase":
        self._run(*args, **kwargs)
        self.__done = True
        return self

    @property
    def done(self) -> bool:
        return self.__done
