from abc import abstractmethod
from pathlib import Path

from libcoffee.common.executablebase import ExecutableBase
from libcoffee.common.path import SDFFile


class StateGeneratorBase(ExecutableBase):

    def __init__(
        self,
        exec: Path,
        n_jobs: int = 1,
        enum_ionization: bool = True,
        pH: float = 7.0,
        max_states: int = 32,
        verbose: bool = False,
    ):
        super().__init__(exec, verbose)
        self._n_jobs = n_jobs
        self._enum_ionization = enum_ionization
        self._pH = pH
        self._max_states = max_states

    @abstractmethod
    def save(self, path: SDFFile) -> "StateGeneratorBase": ...
