import shutil
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Self


class ExecutableBase(ABC):

    def __validate_exec(self: Self, exec: Path) -> None:
        if shutil.which(exec) is None:
            raise FileNotFoundError(f"Executable '{exec}' is not found")

    def __init__(self: Self, exec: Path, verbose: bool = False):
        self._verbose = verbose
        self._exec = exec
        self.__done = False
        self.__validate_exec(exec)

    @abstractmethod
    def _run(self: Self, *args: Any, **kwargs: Any) -> None: ...

    def run(self: Self, *args: Any, **kwargs: Any) -> "ExecutableBase":
        self._run(*args, **kwargs)
        self.__done = True
        return self

    @property
    def done(self) -> bool:
        return self.__done
