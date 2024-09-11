from abc import abstractmethod
from pathlib import Path

from libcoffee.common.executablebase import ExecutableBase
from libcoffee.common.path import SDFFile


class ConformerGeneratorBase(ExecutableBase):

    def __init__(
        self,
        exec: Path,
        n_jobs: int = 1,
        max_confs: int = 200,
        min_rmsd: float = 0.5,
        energy_tolerance: float = 10.0,
        verbose: bool = False,
    ):
        super().__init__(exec, verbose)
        self._n_jobs = n_jobs
        self._max_confs = max_confs
        self._min_rmsd = min_rmsd
        self._energy_tolerance = energy_tolerance

    @abstractmethod
    def save(self, path: SDFFile) -> "ConformerGeneratorBase": ...
