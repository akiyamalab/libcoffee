import os
import subprocess
import tempfile
from pathlib import Path

from .config import REstrettoConfig

__PATH_ATOMGRID_GEN = f"{os.path.dirname(__file__)}/atomgrid-gen"
__PATH_CONFORMER_DOCKING = f"{os.path.dirname(__file__)}/conformer-docking"


def _generate_atomgrid(config: REstrettoConfig, verbose=False) -> tuple[str, str]:
    """
    execute atomgrid-gen based on the given config file.
    """
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as fout:
        fout.write(str(config))
        config_path = Path(fout.name)
        ret = subprocess.run(
            [__PATH_ATOMGRID_GEN, str(config_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
    return ret.stdout.decode("utf-8").strip(), ret.stderr.decode("utf-8").strip()


def _dock_cmpds(config: REstrettoConfig, verbose=False) -> tuple[str, str]:
    """
    conformer-dockingを実行する。
    """
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as fout:
        fout.write(str(config))
        config_path = Path(fout.name)
        ret = subprocess.run(
            [__PATH_CONFORMER_DOCKING, str(config_path)],
            stdout=(subprocess.PIPE if verbose else subprocess.DEVNULL),
            stderr=subprocess.PIPE,
            check=True,
        )
    return ret.stdout.decode("utf-8").strip(), ret.stderr.decode("utf-8").strip()


class REstretto:
    def __init__(self, config: REstrettoConfig, verbose: bool = False):
        self.config = config
        self.verbose = verbose
        self.done = False
        self.stdout = {"atom_grid": "", "docking": ""}
        self.stderr = {"atom_grid": "", "docking": ""}

    def run(self) -> "REstretto":
        try:
            self.stdout["atom_grid"], self.stderr["atom_grid"] = _generate_atomgrid(self.config, verbose=self.verbose)
            self.stdout["docking"], self.stderr["docking"] = _dock_cmpds(self.config, verbose=self.verbose)
            self.done = True
        except subprocess.CalledProcessError as e:
            # TODO treat exceptions more properly
            print(f"Failed to execute {e.cmd}:")
            print(f"  {e.stderr.decode('utf-8')}")
            print(f"Configs are:")
            print(str(self.config))
            raise e
        return self

    @property
    def result(self) -> Path:
        if not self.done:
            raise ValueError("REstretto has not been run yet.")
        return self.config.output
