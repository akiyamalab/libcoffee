import tempfile
from pathlib import Path

import numpy as np
import numpy.typing as npt


class REstrettoConfig:
    def __init__(
        self,
        innerbox: npt.NDArray[np.int_] = np.array([10, 10, 10], dtype=np.int_),
        outerbox: npt.NDArray[np.int_] = np.array([20, 20, 20], dtype=np.int_),
        box_center: npt.NDArray[np.float_] = np.array([0, 0, 0], dtype=np.float_),
        receptor: Path | None = None,
        ligands: list[Path] = [],
        output: Path | None = None,
        search_pitch: npt.NDArray[np.float_] = np.array([1.0, 1.0, 1.0], dtype=np.float_),
        scoring_pitch: npt.NDArray[np.float_] = np.array([0.25, 0.25, 0.25], dtype=np.float_),
        memory_size: int = 8000,
        grid_folder: Path | None = None,
        no_local_opt: bool = False,
        poses_per_lig: int = 1,
        pose_rmsd: float = 2.0,
        make_parent_dirs: bool = True,
        overwrite: bool = False,
        output_score_threshold: float = -3.0,
        poses_per_lig_before_opt: int = 2000,
        log: Path | None = None,
    ):
        self.innerbox = innerbox
        self.outerbox = outerbox
        self.box_center = box_center
        self.receptor = receptor
        self.ligands = ligands
        self.output = output
        self.search_pitch = search_pitch
        self.scoring_pitch = scoring_pitch
        self.memory_size = memory_size
        self.grid_folder = grid_folder
        self.no_local_opt = no_local_opt
        self.poses_per_lig = poses_per_lig
        self.pose_rmsd = pose_rmsd
        self.make_parent_dirs = make_parent_dirs
        self.overwrite = overwrite
        self.output_score_threshold = output_score_threshold
        self.poses_per_lig_before_opt = poses_per_lig_before_opt
        self.log = log

    @property
    def innerbox(self) -> npt.NDArray[np.int_]:
        return self._innerbox

    @innerbox.setter
    def innerbox(self, value: npt.NDArray[np.int_]) -> None:
        if value.shape != (3,):
            raise ValueError("innerbox must be a 3D information.")
        if not np.all(value > 0):
            raise ValueError("innerbox must be positive.")
        self._innerbox = value

    @property
    def outerbox(self) -> npt.NDArray[np.int_]:
        return self._outerbox

    @outerbox.setter
    def outerbox(self, value: npt.NDArray[np.int_]) -> None:
        if value.shape != (3,):
            raise ValueError("outerbox must be a 3D information.")
        if not np.all(value > 0):
            raise ValueError("outerbox must be positive.")
        self._outerbox = value

    @property
    def receptor(self) -> Path | None:
        return self._receptor

    @receptor.setter
    def receptor(self, value: Path | None) -> None:
        if value is not None and not value.exists():
            raise FileNotFoundError(f"Receptor file not found: {value}")
        self._receptor = value

    @property
    def ligands(self) -> list[Path]:
        return self._ligands

    @ligands.setter
    def ligands(self, value: list[Path]) -> None:
        for elem in value:
            if not elem.exists():
                raise FileNotFoundError(f"Ligand file not found: {elem}")
        self._ligands = value

    @property
    def output(self) -> Path | None:
        return self._output

    @output.setter
    def output(self, value: Path | None) -> None:
        if value is not None and value.exists() and not self.overwrite:
            raise FileExistsError(f"Output file already exists: {value}")
        if value is not None and value.is_dir():
            raise IsADirectoryError(f"Output is a directory: {value}")
        self._output: Path | None = value

    @property
    def box_center(self) -> npt.NDArray[np.float_]:
        return self._box_center

    @box_center.setter
    def box_center(self, value: npt.NDArray[np.float_]) -> None:
        if value.shape != (3,):
            raise ValueError("box_center must be a 3D information.")
        self._box_center = value

    @property
    def search_pitch(self) -> npt.NDArray[np.float_]:
        return self._search_pitch

    @search_pitch.setter
    def search_pitch(self, value: npt.NDArray[np.float_]) -> None:
        if value.shape != (3,):
            raise ValueError("search_pitch must be a 3D information.")
        if not np.all(value > 0):
            raise ValueError("search_pitch must be positive.")
        self._search_pitch = value

    @property
    def scoring_pitch(self) -> npt.NDArray[np.float_]:
        return self._scoring_pitch

    @scoring_pitch.setter
    def scoring_pitch(self, value: npt.NDArray[np.float_]) -> None:
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
    def grid_folder(self) -> Path | None:
        return self._grid_folder

    @grid_folder.setter
    def grid_folder(self, value: Path | None) -> None:
        if value is None:
            self._grid_folder = Path(tempfile.mkdtemp())
        else:
            self._grid_folder = value

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
        if self.receptor is None:
            raise ValueError("Receptor file must be specified.")
        if len(self.ligands) == 0:
            raise ValueError("Ligand files must be specified.")
        if self.output is None:
            raise ValueError("Output file must be specified.")

        if self.grid_folder is None:
            self.grid_folder = Path(tempfile.mkdtemp())
        ret: list[str] = []
        ret.append(f"INNERBOX {', '.join(str(elem) for elem in self.innerbox)}")
        ret.append(f"OUTERBOX {', '.join(str(elem) for elem in self.outerbox)}")
        ret.append(f"BOX_CENTER {', '.join(str(elem) for elem in self.box_center)}")
        ret.append(f"RECEPTOR {str(self.receptor)}")
        for ligand in self.ligands:
            ret.append(f"LIGAND {str(ligand)}")
        ret.append(f"OUTPUT {str(self.output)}")

        ret.append(f"SEARCH_PITCH {', '.join(str(elem) for elem in self.search_pitch)}")
        ret.append(f"SCORING_PITCH {', '.join(str(elem) for elem in self.scoring_pitch)}")
        ret.append(f"MEMORY_SIZE {str(self.memory_size)}")
        ret.append(f"GRID_FOLDER {str(self.grid_folder)}")
        ret.append(f"NO_LOCAL_OPT {str(self.no_local_opt)}")
        ret.append(f"POSES_PER_LIG {str(self.poses_per_lig)}")
        ret.append(f"POSE_RMSD {str(self.pose_rmsd)}")
        ret.append(f"POSES_PER_LIG_BEFORE_OPT {str(self.poses_per_lig_before_opt)}")
        ret.append(f"OUTPUT_SCORE_THRESHOLD {str(self.output_score_threshold)}")
        if self.log is not None:
            ret.append(f"LOG {str(self.log)}")
        return "\n".join(ret)

    @staticmethod
    def parse_file(filepath: Path) -> "REstrettoConfig":
        ret = REstrettoConfig()
        with open(filepath, "r") as fin:
            for line in fin:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                k, v = line.split(maxsplit=1)
                if k == "INNERBOX":
                    ret.innerbox = np.array([int(elem) for elem in v.split(",")], dtype=np.int_)
                elif k == "OUTERBOX":
                    ret.outerbox = np.array([int(elem) for elem in v.split(",")], dtype=np.int_)
                elif k == "BOX_CENTER":
                    ret.box_center = np.array([float(elem) for elem in v.split(",")], dtype=np.float_)
                elif k == "SEARCH_PITCH":
                    ret.search_pitch = np.array([float(elem) for elem in v.split(",")], dtype=np.float_)
                elif k == "SCORING_PITCH":
                    ret.scoring_pitch = np.array([float(elem) for elem in v.split(",")], dtype=np.float_)
                elif k == "MEMORY_SIZE":
                    ret.memory_size = int(v)
                elif k == "RECEPTOR":
                    ret.receptor = Path(v)
                elif k == "LIGAND":
                    ret.ligands.append(Path(v))
                elif k == "OUTPUT":
                    ret.output = Path(v)
                elif k == "GRID_FOLDER":
                    ret.grid_folder = Path(v)
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
