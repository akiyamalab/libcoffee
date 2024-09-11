from libcoffee.common.path import SDFFile
from libcoffee.docking.conformer_generator import ConformerGeneratorBase


class Omega(ConformerGeneratorBase):

    def run(self) -> "Omega":
        raise NotImplementedError

    def save(self, path: SDFFile) -> "Omega":
        raise NotImplementedError
