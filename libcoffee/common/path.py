import os
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar


@dataclass(frozen=True, init=False)
class _TypedFile(os.PathLike[str]):
    """
    A file path whose suffix is restricted to `_suffix`.
    This class does not inherit Path; use `.path` for path operations (their results are plain Path objects).
    It can be passed to open(), shutil, subprocess, etc. via os.PathLike.
    """

    _suffix: ClassVar[str]
    _format: ClassVar[str]
    path: Path

    def __init__(self, file: str | os.PathLike[str]) -> None:
        path = Path(file)
        if path.suffix != self._suffix:
            raise ValueError(f"File {file} is not a {self._format} file")
        object.__setattr__(self, "path", path)

    def __fspath__(self) -> str:
        return os.fspath(self.path)

    def __str__(self) -> str:
        # used when embedded in command lines via str() or f-strings
        return str(self.path)


class SDFFile(_TypedFile):
    _suffix = ".sdf"
    _format = "SDF"


class PDBFile(_TypedFile):
    _suffix = ".pdb"
    _format = "PDB"


class FBDBFile(_TypedFile):
    _suffix = ".fbdb"
    _format = "FBDB"
