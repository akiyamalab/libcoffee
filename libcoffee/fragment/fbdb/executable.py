from libcoffee.common.path import SDFFile


class FBDB:
    def __init__(self, path: SDFFile, verbose: bool = False):
        raise NotImplementedError

    def construct(self):  # type: ignore
        raise NotImplementedError

    def search(self, query):  # type: ignore
        raise NotImplementedError
