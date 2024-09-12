from pathlib import Path


class SDFFile(Path):
    _flavour = type(Path())._flavour

    def __init__(self, file: Path):
        super().__init__()
        if self.suffix != ".sdf":
            raise ValueError(f"File {file} is not an SDF file")


class PDBFile(Path):
    _flavour = type(Path())._flavour

    def __init__(self, file: Path):
        super().__init__()
        if self.suffix != ".pdb":
            raise ValueError(f"File {file} is not a PDB file")


class FBDBFile(Path):
    _flavour = type(Path())._flavour

    def __init__(self, file: Path):
        super().__init__()
        if self.suffix != ".fbdb":
            raise ValueError(f"File {file} is not a FBDB file")
