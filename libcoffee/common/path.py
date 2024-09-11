from pathlib import Path


class SDFFile(Path):
    def __init__(self, file: Path):
        super().__init__(file)
        if self.suffix != ".sdf":
            raise ValueError(f"File {file} is not an SDF file")


class PDBFile(Path):
    def __init__(self, file: Path):
        super().__init__(file)
        if self.suffix != ".pdb":
            raise ValueError(f"File {file} is not a PDB file")


class FBDBFile(Path):
    def __init__(self, file: Path):
        super().__init__(file)
        if self.suffix != ".fbdb":
            raise ValueError(f"File {file} is not a FBDB file")
