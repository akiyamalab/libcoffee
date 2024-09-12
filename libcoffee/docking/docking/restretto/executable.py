import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import numpy as np
import numpy.typing as npt

from libcoffee.common.executablebase import ExecutableBase
from libcoffee.common.path import PDBFile, SDFFile
from libcoffee.docking.docking_region import DockingRegion

from .config import _REstrettoConfig

__PATH_ATOMGRID_GEN = f"{os.path.dirname(__file__)}/atomgrid-gen"
__PATH_CONFORMER_DOCKING = f"{os.path.dirname(__file__)}/conformer-docking"


def _generate_atomgrid(
    config: _REstrettoConfig, proteinfile: PDBFile, grid_folder: Path, verbose: bool = False
) -> tuple[str, str]:
    """
    execute atomgrid-gen based on the given config file.
    """
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as fout:
        fout.write(str(config))
        fout.flush()
        config_path = Path(fout.name)
        ret = subprocess.run(
            [__PATH_ATOMGRID_GEN, str(config_path), "-r", str(proteinfile), "-g", str(grid_folder)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
    return ret.stdout.decode("utf-8").strip(), ret.stderr.decode("utf-8").strip()


def _dock_cmpds(
    config: _REstrettoConfig,
    proteinfile: PDBFile,
    ligandfile: SDFFile,
    grid_folder: Path,
    outputfile: SDFFile,
    verbose=False,
) -> tuple[str, str]:
    """
    conformer-dockingを実行する。
    """
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as fout:
        fout.write(str(config))
        fout.flush()
        config_path = Path(fout.name)
        ret = subprocess.run(
            [
                __PATH_CONFORMER_DOCKING,
                str(config_path),
                "-r",
                str(proteinfile),
                "-l",
                str(ligandfile),
                "-g",
                str(grid_folder),
                "-o",
                str(outputfile),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
    return ret.stdout.decode("utf-8").strip(), ret.stderr.decode("utf-8").strip()


# receptor: Path | None = None,
# ligands: list[Path] = [],
# output: Path | None = None,
# grid_folder: Path | None = None,
#
class REstretto(ExecutableBase):
    def __init__(
        self,
        docking_region: DockingRegion,
        outerbox: Optional[npt.NDArray[np.int32]] = None,
        scoring_pitch: npt.NDArray[np.float64] = np.array([0.25, 0.25, 0.25], dtype=np.float64),
        memory_size: int = 8000,
        no_local_opt: bool = False,
        poses_per_lig: int = 1,
        pose_rmsd: float = 2.0,
        output_score_threshold: float = -3.0,
        poses_per_lig_before_opt: int = 2000,
        log: Path | None = None,
        verbose: bool = False,
    ):
        self.__config = _REstrettoConfig(
            docking_region=docking_region,
            outerbox=outerbox,
            scoring_pitch=scoring_pitch,
            memory_size=memory_size,
            no_local_opt=no_local_opt,
            poses_per_lig=poses_per_lig,
            pose_rmsd=pose_rmsd,
            output_score_threshold=output_score_threshold,
            poses_per_lig_before_opt=poses_per_lig_before_opt,
            log=log if log is not None else Path("/dev/null"),
        )
        self.__grid_folder = tempfile.TemporaryDirectory()
        self.__outputfile = tempfile.NamedTemporaryFile(suffix=".sdf")
        self.verbose = verbose
        self.stdout = {"atom_grid": "", "docking": ""}
        self.stderr = {"atom_grid": "", "docking": ""}

    def _run(self, proteinfile: PDBFile, ligandfile: SDFFile) -> None:
        try:
            self.stdout["atom_grid"], self.stderr["atom_grid"] = _generate_atomgrid(
                self.__config, proteinfile, Path(self.__grid_folder.name), verbose=self.verbose
            )
            self.stdout["docking"], self.stderr["docking"] = _dock_cmpds(
                self.__config,
                proteinfile,
                ligandfile,
                Path(self.__grid_folder.name),
                SDFFile(self.__outputfile.name),
                verbose=self.verbose,
            )
        except subprocess.CalledProcessError as e:
            # TODO treat exceptions more properly
            print(f"Failed to execute {e.cmd}:")
            print(f"  {e.stderr.decode('utf-8')}")
            print(f"Configs are:")
            print(str(self.__config))
            raise e

    def save(self, outputfile: SDFFile) -> "REstretto":
        shutil.copy(self.__outputfile.name, outputfile)
        return self
