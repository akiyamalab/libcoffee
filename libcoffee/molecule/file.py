from tempfile import NamedTemporaryFile, _TemporaryFileWrapper
from typing import Generator

from libcoffee.common.path import SDFFile


def SubSDFFileGenerator(file: SDFFile, n_cmpds_per_file: int = -1) -> Generator[_TemporaryFileWrapper[str], None, None]:
    """
    SDFファイルを分割するためのジェネレータ。tempfile.NamedTemporaryFileを返す。
    """
    with open(file) as fin:
        n_cmpds = 0
        fout = NamedTemporaryFile(mode="w", suffix=".sdf")
        for line in fin:
            fout.write(line)
            if line.strip() == "$$$$":
                n_cmpds += 1
                if n_cmpds_per_file != -1 and n_cmpds % n_cmpds_per_file == 0:
                    fout.flush()
                    yield fout
                    fout = NamedTemporaryFile(mode="w", suffix=".sdf")
        fout.flush()
        yield fout
