import os
import subprocess
from pathlib import Path

from .config import CmpdDecomposeConfig

__PATH_DECOMPOSE = f"{os.path.dirname(__file__)}/decompose"


def _decompose(config: CmpdDecomposeConfig, verbose=False) -> tuple[str, str]:
    """
    execute decompose based on the given config file.
    """
    ret = subprocess.run(
        [__PATH_DECOMPOSE, str(config)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    return ret.stdout.decode("utf-8").strip(), ret.stderr.decode("utf-8").strip()


class CmpdDecompose:
    def __init__(self, config: CmpdDecomposeConfig, verbose: bool = False):
        self.config = config
        self.verbose = verbose
        self.done = False
        self.stdout = ""
        self.stderr = ""

    def run(self) -> "CmpdDecompose":
        try:
            self.stdout, self.stderr = _decompose(self.config, verbose=self.verbose)
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
            raise ValueError("Compound decomposition has not been run yet.")
        return self.config.output_path
