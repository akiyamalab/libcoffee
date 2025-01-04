from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Self


class ExecutableBase(ABC):

    def __init__(self: Self, exec: Path, verbose: bool = False):
        self._verbose = verbose
        self._exec = exec
        self.__done = False

    @abstractmethod
    def _run(self: Self, *args: Any, **kwargs: Any) -> None: ...

    def run(self: Self, *args: Any, **kwargs: Any) -> "ExecutableBase":
        self._run(*args, **kwargs)
        self.__done = True
        return self

    @property
    def done(self) -> bool:
        return self.__done
