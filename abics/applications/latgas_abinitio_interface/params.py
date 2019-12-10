import os

from ...util import expand_path

class DFTParams:
    def __init__(self):
        pass

    @classmethod
    def from_dict(cls, d):
        if 'solver' in d:
            d = d['solver']
        params = cls()
        params.base_input_dir = expand_path(d.get('base_input_dir', './baseinput'), os.getcwd())
        params.solver = d['type']
        params.path = expand_path(d['path'], os.getcwd())
        params.perturb = d.get('perturb', 0.1)
        params.solver_run_scheme = d.get('run_scheme',
                                         'mpi_spawn_ready')
        return params

    @classmethod
    def from_toml(cls, f):
        import toml
        return cls.from_dict(toml.load(f))
