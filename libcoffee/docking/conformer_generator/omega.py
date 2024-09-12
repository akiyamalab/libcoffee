import shutil
import subprocess
import tempfile

from libcoffee.common.path import SDFFile
from libcoffee.docking.conformer_generator import ConformerGeneratorBase


class Omega(ConformerGeneratorBase):

    def run(self, file: SDFFile) -> "Omega":
        self.__outputfile = tempfile.NamedTemporaryFile(suffix=".sdf")
        subprocess.run(
            [
                str(self.exec),
                "-mpi_np",
                self._n_jobs,
                "-in",
                file,
                "-out",
                self.__outputfile,
                "-maxConfs",
                self._max_confs,
                "-rms",
                self._min_rmsd,
                "-eWindoew",
                self._energy_tolerance,
            ],
            check=True,
        )
        return self

    def save(self, path: SDFFile) -> "Omega":
        shutil.copy(self.__outputfile.name, path)
        return self
