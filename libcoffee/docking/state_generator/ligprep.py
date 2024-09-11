from libcoffee.common.path import SDFFile
from libcoffee.docking.state_generator import StateGeneratorBase


class Ligprep(StateGeneratorBase):

    def run(self) -> "Ligprep":
        raise NotImplementedError

    def save(self, path: SDFFile) -> "Ligprep":
        raise NotImplementedError
