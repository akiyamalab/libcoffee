from pathlib import Path
from typing import Self


class SDFFile(Path):
    _flavour = type(Path())._flavour  # type: ignore[attr-defined]

    def __new__(cls, file: Path | str) -> Self:
        if Path(file).suffix != ".sdf":
            raise ValueError(f"File {file} is not a SDF file")
        return super().__new__(cls, file)


class PDBFile(Path):
    _flavour = type(Path())._flavour  # type: ignore[attr-defined]

    def __new__(cls, file: Path | str) -> Self:
        if Path(file).suffix != ".pdb":
            raise ValueError(f"File {file} is not a PDB file")
        return super().__new__(cls, file)


class FBDBFile(Path):
    _flavour = type(Path())._flavour  # type: ignore[attr-defined]

    def __new__(cls, file: Path | str) -> Self:
        if Path(file).suffix != ".fbdb":
            raise ValueError(f"File {file} is not a FBDB file")
        return super().__new__(cls, file)
