from abc import abstractmethod
from pathlib import Path
from typing import Sequence

from libcoffee.common.executablebase import ExecutableBase
from libcoffee.common.path import PDBFile, SDFFile
from libcoffee.docking.docking_region import DockingRegion


class DockingToolBase(ExecutableBase):

    def __init__(
        self,
        exec: Path,
        receptor: PDBFile,
        ligand: SDFFile,
        region: DockingRegion | Sequence[DockingRegion],
        n_jobs: int = 1,
        verbose: bool = False,
    ):
        super().__init__(exec, verbose)
        self._n_jobs = n_jobs
        self._receptor = receptor
        self._ligand = ligand
        self._regions = tuple(region if isinstance(region, Sequence) else [region])

    @abstractmethod
    def save(self, path: SDFFile) -> "DockingToolBase": ...
