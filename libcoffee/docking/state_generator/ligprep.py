import shutil
import subprocess
from tempfile import NamedTemporaryFile

from libcoffee.common.path import SDFFile
from libcoffee.docking.state_generator.stategeneratorbase import StateGeneratorBase

_inputfile_format = """
INPUT_FILE_NAME   {inputsdf}
OUT_SD   {outputsdf}
FORCE_FIELD   16
EPIK   {enum_ionization_yes_or_no}
DETERMINE_CHIRALITIES   no
IGNORE_CHIRALITIES   no
NUM_STEREOISOMERS   {max_states}
"""


class Ligprep(StateGeneratorBase):

    def _run(self, file: SDFFile) -> None:
        self.__inputfile = NamedTemporaryFile(suffix=".inp")
        self.__outputfile = NamedTemporaryFile(dir=".", prefix=".", suffix=".sdf")
        # Ligprep requires output file must under the currect working directory
        with open(self.__inputfile.name, "w") as f:
            f.write(
                _inputfile_format.format(
                    inputsdf=file,
                    outputsdf=self.__outputfile.name,
                    enum_ionization_yes_or_no="yes" if self._enum_ionization else "no",
                    max_states=self._max_states,
                )
            )
            f.flush()
        n_subtask = self._n_jobs * 10 if self._n_jobs > 1 else 1  # 10 is a magic number
        subprocess.run(
            " ".join(
                [
                    str(self._exec),
                    "-inp",
                    self.__inputfile.name,
                    "-NJOBS",
                    str(n_subtask),
                    "-HOST",
                    f"localhost:{self._n_jobs}",
                    "-WAIT",
                ]
            ),
            check=True,
            shell=True,
        )
        # TODO: raise error with log file content if ligprep failed

    def save(self, path: SDFFile) -> "Ligprep":
        # dropped-file "h8iu2tb1-dropped.sdf"
        # log file "h8iu2tb1.log"
        shutil.copy(self.__outputfile.name, path)
        return self
