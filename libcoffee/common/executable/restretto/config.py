import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import numpy.typing as npt


@dataclass
class REstrettoConfig:
    innerbox: npt.NDArray[np.int_] = field(default_factory=lambda: np.array([10, 10, 10], dtype=np.int_))
    outerbox: npt.NDArray[np.int_] = field(default_factory=lambda: np.array([20, 20, 20], dtype=np.int_))
    box_center: npt.NDArray[np.float_] = field(default_factory=lambda: np.array([0, 0, 0], dtype=np.float_))
    receptor: Path = Path()
    ligands: list[Path] = field(default_factory=list)
    output: Path = Path()

    search_pitch: npt.NDArray[np.float_] = field(default_factory=lambda: np.array([1.0, 1.0, 1.0], dtype=np.float_))
    scoring_pitch: npt.NDArray[np.float_] = field(default_factory=lambda: np.array([0.25, 0.25, 0.25], dtype=np.float_))
    memory_size: int = 8000
    grid_folder: Path | None = None
    no_local_opt: bool = False
    poses_per_lig: int = 1
    pose_rmsd: float = 2.0
    make_parent_dirs: bool = True  # TODO outputする先のディレクトリが存在しない時に、ディレクトリを作成するかどうか
    overwrite: bool = False  # TODO outputする先にファイルが存在する時に、上書きするかどうか
    output_score_threshold: float = -3.0
    poses_per_lig_before_opt: int = 2000
    log: Path | None = None

    def __str__(self) -> str:
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
        ret = REstrettoConfig(ligands=[])
        with open(filepath, "r") as fin:
            for line in fin:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                k, v = line.split(maxsplit=1)
                if k == "INNERBOX":
                    ret.innerbox = np.array([int(elem) for elem in v.split(",")])
                elif k == "OUTERBOX":
                    ret.outerbox = np.array([int(elem) for elem in v.split(",")])
                elif k == "BOX_CENTER":
                    ret.box_center = np.array([float(elem) for elem in v.split(",")])
                elif k == "SEARCH_PITCH":
                    ret.search_pitch = np.array([float(elem) for elem in v.split(",")])
                elif k == "SCORING_PITCH":
                    ret.scoring_pitch = np.array([float(elem) for elem in v.split(",")])
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
                elif k == "OUTPUT_SCORE_THRESHOLD":
                    ret.output_score_threshold = float(v)
                else:
                    raise KeyError(f"Unknown parameter: {k}")
        return ret
