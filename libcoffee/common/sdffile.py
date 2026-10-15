from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import IO, Generator

from libcoffee.common.path import SDFFile


def SubSDFFileGenerator(
    file: str | Path | SDFFile, n_cmpds_per_file: int = -1
) -> Generator[IO[str], None, None]:
    """
    SDFファイルを分割するためのジェネレータ。tempfile.NamedTemporaryFileを返す。
    """
    file = SDFFile(file)
    with open(file) as fin:
        n_cmpds = 0
        fout = NamedTemporaryFile(mode="w", suffix=".sdf")
        has_content = False  # whether fout has non-blank lines that are not yielded yet
        for line in fin:
            fout.write(line)
            has_content = has_content or bool(line.strip())
            if line.strip() == "$$$$":
                n_cmpds += 1
                if n_cmpds_per_file != -1 and n_cmpds % n_cmpds_per_file == 0:
                    fout.flush()
                    yield fout
                    fout = NamedTemporaryFile(mode="w", suffix=".sdf")
                    has_content = False
        if has_content:
            fout.flush()
            yield fout
        else:
            fout.close()
