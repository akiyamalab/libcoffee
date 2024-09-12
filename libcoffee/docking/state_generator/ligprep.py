import shutil
import subprocess
from tempfile import NamedTemporaryFile

from libcoffee.common.path import SDFFile
from libcoffee.docking.state_generator import StateGeneratorBase

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

    def run(self, file: SDFFile) -> "Ligprep":  # type: ignore[override]
        self.__inputfile = NamedTemporaryFile(suffix=".inp")
        self.__outputfile = NamedTemporaryFile(suffix=".sdf")
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
            [
                str(self._exec),
                self.__inputfile.name,
                "-NJOBS",
                str(n_subtask),
                "-HOST",
                f"localhost:{self._n_jobs}",
                "WAIT",
            ],
            check=True,
        )
        return self

    def save(self, path: SDFFile) -> "Ligprep":
        shutil.copy(self.__outputfile.name, path)
        return self
