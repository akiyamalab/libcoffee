from dataclasses import dataclass
from pathlib import Path


@dataclass
class CmpdDecomposeConfig:
    fragment_path: Path = Path()
    ligand_path: Path = Path()
    output_path: Path = Path()
    log: Path | None = None
    capping_atomic_num: int = -1
    enable_carbon_capping: bool = False
    ins_fragment_id: bool = False
    max_ring_size: int = -1
    no_merge_solitary: bool = False

    def __str__(self) -> str:
        options = []
        options.append(f"-f {self.fragment_path}")
        options.append(f"-l {self.ligand_path}")
        options.append(f"-o {self.output_path}")
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
