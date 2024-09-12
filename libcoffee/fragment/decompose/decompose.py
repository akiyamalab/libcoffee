import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Optional

from libcoffee.common.executablebase import ExecutableBase
from libcoffee.common.path import SDFFile

__PATH_DECOMPOSE = f"{os.path.dirname(__file__)}/decompose"


@dataclass(frozen=True)
class _CmpdDecomposeConfig:
    log: Optional[Path] = None
    capping_atomic_num: int = -1
    enable_carbon_capping: bool = False
    ins_fragment_id: bool = False
    max_ring_size: int = -1
    no_merge_solitary: bool = False

    def __str__(self) -> str:
        options = []
        if self.log is not None:
            options.append(f"--log {self.log}")
        if self.capping_atomic_num != -1:
            options.append(f"--capping_atomic_num {self.capping_atomic_num}")
        if self.enable_carbon_capping:
            options.append("--enable_carbon_capping")
        if self.ins_fragment_id:
            options.append("--ins_fragment_id")
        if self.max_ring_size != -1:
            options.append(f"--max_ring_size {self.max_ring_size}")
        if self.no_merge_solitary:
            options.append("--no_merge_solitary")
        return " ".join(options)


def _decompose(
    config: _CmpdDecomposeConfig, ligand_file: SDFFile, fragment_file: SDFFile, annotated_file: SDFFile, verbose=False
) -> tuple[str, str]:
    """
    execute decompose based on the given config file.
    """
    ret = subprocess.run(
        __PATH_DECOMPOSE + " " + f"-l {ligand_file} " + f"-f {fragment_file} " + f"-o {annotated_file} " + str(config),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
        shell=True,
    )
    return ret.stdout.decode("utf-8").strip(), ret.stderr.decode("utf-8").strip()


class CmpdDecompose(ExecutableBase):
    def __init__(
        self,
        log: Optional[Path] = None,
        capping_atomic_num: int = -1,
        enable_carbon_capping: bool = False,
        ins_fragment_id: bool = False,
        max_ring_size: int = -1,
        no_merge_solitary: bool = False,
        verbose: bool = False,
    ):
        self.__config = _CmpdDecomposeConfig(
            log=log if log is not None else Path("/dev/null"),
            capping_atomic_num=capping_atomic_num,
            enable_carbon_capping=enable_carbon_capping,
            ins_fragment_id=ins_fragment_id,
            max_ring_size=max_ring_size,
            no_merge_solitary=no_merge_solitary,
        )
        self._verbose = verbose
        self.stdout = ""
        self.stderr = ""

    def _run(self, file: SDFFile) -> None:
        self.__fragment_file = NamedTemporaryFile(suffix=".sdf")
        self.__annotated_file = NamedTemporaryFile(suffix=".sdf")
        try:
            self.stdout, self.stderr = _decompose(
                self.__config,
                ligand_file=file,
                fragment_file=SDFFile(self.__fragment_file.name),
                annotated_file=SDFFile(self.__annotated_file.name),
                verbose=self._verbose,
            )
        except subprocess.CalledProcessError as e:
            # TODO treat exceptions more properly
            print(f"Failed to execute {e.cmd}:")
            print(f"  {e.stderr.decode('utf-8')}")
            print(f"Configs are:")
            print(str(self.config))
            raise e

    def save(self, fragmentfile: SDFFile, annotatedfile: SDFFile) -> "CmpdDecompose":
        shutil.copy(self.__fragment_file.name, fragmentfile)
        shutil.copy(self.__annotated_file.name, annotatedfile)
        return self
