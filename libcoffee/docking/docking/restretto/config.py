import tempfile
from pathlib import Path
from typing import Optional

import numpy as np
import numpy.typing as npt

from libcoffee.docking.docking_region import DockingRegion


class _REstrettoConfig:
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
    ):
        self.__docking_region = docking_region
        self.outerbox = outerbox if outerbox is not None else np.array(docking_region.size) + 10
        self.scoring_pitch = scoring_pitch
        self.memory_size = memory_size
        self.no_local_opt = no_local_opt
        self.poses_per_lig = poses_per_lig
        self.pose_rmsd = pose_rmsd
        self.output_score_threshold = output_score_threshold
        self.poses_per_lig_before_opt = poses_per_lig_before_opt
        self.log = log if log is not None else Path("/dev/null")

    @property
    def outerbox(self) -> npt.NDArray[np.int32]:
        return self._outerbox

    @outerbox.setter
    def outerbox(self, value: npt.NDArray[np.int32]) -> None:
        if value.shape != (3,):
            raise ValueError("outerbox must be a 3D information.")
        if not np.all(value > 0):
            raise ValueError("outerbox must be positive.")
        self._outerbox = value

    @property
    def scoring_pitch(self) -> npt.NDArray[np.float64]:
        return self._scoring_pitch

    @scoring_pitch.setter
    def scoring_pitch(self, value: npt.NDArray[np.float64]) -> None:
        if value.shape != (3,):
            raise ValueError("scoring_pitch must be a 3D information.")
        if not np.all(value > 0):
            raise ValueError("scoring_pitch must be positive.")
        self._scoring_pitch = value

    @property
    def memory_size(self) -> int:
        return self._memory_size

    @memory_size.setter
    def memory_size(self, value: int) -> None:
        if value <= 0:
            raise ValueError("memory_size must be positive.")
        self._memory_size = value

    @property
    def poses_per_lig(self) -> int:
        return self._poses_per_lig

    @poses_per_lig.setter
    def poses_per_lig(self, value: int) -> None:
        if value <= 0:
            raise ValueError("poses_per_lig must be positive.")
        self._poses_per_lig = value

    @property
    def pose_rmsd(self) -> float:
        return self._pose_rmsd

    @pose_rmsd.setter
    def pose_rmsd(self, value: float) -> None:
        if value <= 0:
            raise ValueError("pose_rmsd must be positive.")
        self._pose_rmsd = value

    @property
    def poses_per_lig_before_opt(self) -> int:
        return self._poses_per_lig_before_opt

    @poses_per_lig_before_opt.setter
    def poses_per_lig_before_opt(self, value: int) -> None:
        if value <= 0:
            raise ValueError("poses_per_lig_before_opt must be positive.")
        self._poses_per_lig_before_opt = value

    @property
    def log(self) -> Path | None:
        return self._log

    @log.setter
    def log(self, value: Path | None) -> None:
        if value is not None and value.is_dir():
            raise IsADirectoryError(f"Log is a directory: {value}")
        self._log = value

    def __str__(self) -> str:

        ret: list[str] = []
        ret.append(f"INNERBOX {', '.join(str(elem) for elem in self.__docking_region.size)}")
        ret.append(f"OUTERBOX {', '.join(str(elem) for elem in self.outerbox)}")
        ret.append(f"BOX_CENTER {', '.join(str(elem) for elem in self.__docking_region.center)}")
        ret.append(f"SEARCH_PITCH {', '.join(str(elem) for elem in self.__docking_region.pitch)}")
        ret.append(f"SCORING_PITCH {', '.join(str(elem) for elem in self.scoring_pitch)}")
        ret.append(f"MEMORY_SIZE {str(self.memory_size)}")
        ret.append(f"NO_LOCAL_OPT {str(self.no_local_opt)}")
        ret.append(f"POSES_PER_LIG {str(self.poses_per_lig)}")
        ret.append(f"POSE_RMSD {str(self.pose_rmsd)}")
        ret.append(f"POSES_PER_LIG_BEFORE_OPT {str(self.poses_per_lig_before_opt)}")
        ret.append(f"OUTPUT_SCORE_THRESHOLD {str(self.output_score_threshold)}")
        if self.log is not None:
            ret.append(f"LOG {str(self.log)}")
        return "\n".join(ret)

    @staticmethod
    def parse_file(filepath: Path) -> "_REstrettoConfig":
        ret = _REstrettoConfig()
        with open(filepath, "r") as fin:
            for line in fin:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                k, v = line.split(maxsplit=1)
                if k == "INNERBOX":
                    ret.innerbox = np.array([int(elem) for elem in v.split(",")], dtype=np.int32)
                elif k == "OUTERBOX":
                    ret.outerbox = np.array([int(elem) for elem in v.split(",")], dtype=np.int32)
                elif k == "BOX_CENTER":
                    ret.box_center = np.array([float(elem) for elem in v.split(",")], dtype=np.float64)
                elif k == "SEARCH_PITCH":
                    ret.search_pitch = np.array([float(elem) for elem in v.split(",")], dtype=np.float64)
                elif k == "SCORING_PITCH":
                    ret.scoring_pitch = np.array([float(elem) for elem in v.split(",")], dtype=np.float64)
                elif k == "MEMORY_SIZE":
                    ret.memory_size = int(v)
                elif k == "RECEPTOR":
                    ret.receptor = Path(v)
                elif k == "LIGAND":
                    ret.ligands.append(Path(v))
                elif k == "OUTPUT":
                    ret.output = Path(v)
                elif k == "NO_LOCAL_OPT":
                    ret.no_local_opt = bool(v)
                elif k == "POSES_PER_LIG":
                    ret.poses_per_lig = int(v)
                elif k == "POSE_RMSD":
                    ret.pose_rmsd = float(v)
                elif k == "POSES_PER_LIG_BEFORE_OPT":
                    ret.poses_per_lig_before_opt = int(v)
                elif k == "OUTPUT_SCORE_THRESHOLD":
                    ret.output_score_threshold = float(v)
                else:
                    raise KeyError(f"Unknown parameter: {k}")
        return ret
