import shutil
import subprocess
import tempfile

from libcoffee.common.path import SDFFile
from libcoffee.docking.conformer_generator.conformergeneratorbase import ConformerGeneratorBase


class Omega(ConformerGeneratorBase):

    def _run(self, file: SDFFile) -> None:
        self.__outputfile = tempfile.NamedTemporaryFile(suffix=".sdf")
        subprocess.run(
            [
                str(self._exec),
                "-mpi_np",
                str(self._n_jobs),
                "-in",
                file,
                "-out",
                self.__outputfile.name,
                "-maxConfs",
                str(self._max_confs),
                "-rms",
                str(self._min_rmsd),
                "-eWindow",
                str(self._energy_tolerance),
            ],
            check=True,
        )

    def save(self, path: SDFFile) -> "Omega":
        shutil.copy(self.__outputfile.name, path)
        return self
        # TODO: below three files
        # -rw-r--r-- 1 ud02114 tga-pharma 2.1K Sep 12 14:12 oeomega_classic_parm.txt
        # -rw-r--r-- 1 ud02114 tga-pharma  193 Sep 12 14:12 oeomega_classic_rpt.csv
        # -rw-r--r-- 1 ud02114 tga-pharma 4.0K Sep 12 14:12 oeomega_classic_log.txt
