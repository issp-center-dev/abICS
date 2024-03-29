from __future__ import annotations

from typing import Any
from collections.abc import MutableMapping

class SamplerParams:
    """Parameter set for specifying sampling algorithm

    Attributes
    ----------
    sampler : str
        Sampler name
    """

    def __init__(self) -> None:
        self.sampler = "RXMC"

    @classmethod
    def from_dict(cls, d: MutableMapping[str, Any]) -> "SamplerParams":
        """
        Read information from dictionary

        Parameters
        ----------
        d: dict
            Dictionary including parameters for specifying sampling algorithm

        Returns
        -------
        params: SamplerParams object
            self
        """
        if "sampler" in d:
            d = d["sampler"]
        params = cls()
        params.sampler = d.get("sampler", "RXMC")
        return params

    @classmethod
    def from_toml(cls, fname: str) -> "SamplerParams":
        """
        Read information from toml file

        Parameters
        ----------
        f: str
            The name of input toml File

        Returns
        -------
        SamplerParams: SamplerParams object
            self
        """
        import toml

        return cls.from_dict(toml.load(fname))


class ParallelRandomParams:
    """Parameter set for parallel random sampling

    Attributes
    ----------
    nreplicas : int
        The number of replicas
    nprocs_per_replica : int
        The number of processes which a replica uses
    nsteps : int
        The number of MC steps
    sample_frequency :
        The number of MC steps between measurements observables
    print_frequency :
        The number of MC steps between show information
    reload : bool
        Whether to restart simulation or not
    seed : int
        The seed of the random number generator
        If 0, some random number is used (e.g., system time or some random noise).
    """

    def __init__(self) -> None:
        self.nreplicas = 1
        self.nprocs_per_replica = 1
        self.nsteps = None
        self.sample_frequency = 1
        self.print_frequency = 1
        self.reload = False
        self.seed = 0

    @classmethod
    def from_dict(cls, d: MutableMapping[str, Any]) -> "ParallelRandomParams":
        """
        Read information from dictionary

        Parameters
        ----------
        d: dict
            Dictionary including parameters for parallel random sampling

        Returns
        -------
        params: DFTParams object
            self
        """
        if "replica" in d:
            d = d["replica"]
        params = cls()
        params.nreplicas = d["nreplicas"]
        params.nprocs_per_replica = d["nprocs_per_replica"]
        params.nsteps = d["nsteps"]
        params.sample_frequency = d.get("sample_frequency", 1)
        params.print_frequency = d.get("print_frequency", 1)
        params.reload = d.get("reload", False)
        params.seed = d.get("seed", 0)
        return params

    @classmethod
    def from_toml(cls, fname: str) -> "ParallelRandomParams":
        """
        Read information from toml file

        Parameters
        ----------
        f: str
            The name of input toml File

        Returns
        -------
        DFTParams: DFTParams object
            self
        """
        import toml

        return cls.from_dict(toml.load(fname))


class ParallelMCParams:
    """Parameter set for embarrasingly parallel Monte Carlo

    Attributes
    ----------
    nreplicas : int
        The number of replicas
    nprocs_per_replica : int
        The number of processes which a replica uses
    kTstart : float
        The lower bound of temperature range
    kTend : float
        The upper bound of temperature range
    nsteps : int
        The number of MC steps
    sample_frequency :
        The number of MC steps between measurements observables
    print_frequency :
        The number of MC steps between show information
    reload : bool
        Whether to restart simulation or not
    seed : int
        The seed of the random number generator
        If 0, some random number is used (e.g., system time or some random noise).
    """

    def __init__(self) -> None:
        self.nreplicas = -1
        self.nprocs_per_replica = 1
        self.kTstart = -1.0
        self.kTend = -1.0
        self.nsteps = -1
        self.sample_frequency = 1
        self.print_frequency = 1
        self.reload = False
        self.seed = 0

    @classmethod
    def from_dict(cls, d: MutableMapping[str, Any]) -> "ParallelMCParams":
        """
        Read information from dictionary

        Parameters
        ----------
        d: dict
            Dictionary including parameters for embarrassingly parallel Monte Carlo method

        Returns
        -------
        params: DFTParams object
            self
        """
        if "replica" in d:
            d = d["replica"]
        params = cls()
        params.nreplicas = d["nreplicas"]
        params.nprocs_per_replica = d["nprocs_per_replica"]
        params.kTstart = d["kTstart"]
        params.kTend = d["kTend"]
        params.nsteps = d["nsteps"]
        params.sample_frequency = d.get("sample_frequency", 1)
        params.print_frequency = d.get("print_frequency", 1)
        params.reload = d.get("reload", False)
        params.seed = d.get("seed", 0)
        return params

    @classmethod
    def from_toml(cls, fname: str) -> "ParallelMCParams":
        """
        Read information from toml file

        Parameters
        ----------
        f: str
            The name of input toml File

        Returns
        -------
        DFTParams: DFTParams object
            self
        """
        import toml

        return cls.from_dict(toml.load(fname))


class RXParams:
    """Parameter set for replica exchange Monte Carlo

    Attributes
    ----------
    nreplicas : int
        The number of replicas
    nprocs_per_replica : int
        The number of processes which a replica uses
    kTstart : float
        The lower bound of temperature range
    kTend : float
        The upper bound of temperature range
    nsteps : int
        The number of MC steps
    RXtrial_frequency :
        The number of MC steps between replica exchange operations
    sample_frequency :
        The number of MC steps between measurements observables
    print_frequency :
        The number of MC steps between show information
    reload : bool
        Whether to restart simulation or not
    seed : int
        The seed of the random number generator
        If 0, some random number is used (e.g., system time or some random noise).
    """

    def __init__(self) -> None:
        self.nreplicas = None
        self.nprocs_per_replica = 1
        self.kTstart = None
        self.kTend = None
        self.nsteps = None
        self.RXtrial_frequency = 1
        self.sample_frequency = 1
        self.print_frequency = 1
        self.reload = False
        self.seed = 0

    @classmethod
    def from_dict(cls, d: MutableMapping[str, Any]) -> "RXParams":
        """
        Read information from dictionary

        Parameters
        ----------
        d: dict
            Dictionary including parameters for replica exchange Monte Carlo method

        Returns
        -------
        params: DFTParams object
            self
        """
        if "replica" in d:
            d = d["replica"]
        params = cls()
        params.nreplicas = d["nreplicas"]
        params.nprocs_per_replica = d["nprocs_per_replica"]
        params.kTstart = d["kTstart"]
        params.kTend = d["kTend"]
        params.nsteps = d["nsteps"]
        params.RXtrial_frequency = d.get("RXtrial_frequency", 1)
        params.sample_frequency = d.get("sample_frequency", 1)
        params.print_frequency = d.get("print_frequency", 1)
        params.reload = d.get("reload", False)
        params.seed = d.get("seed", 0)
        return params

    @classmethod
    def from_toml(cls, fname: str) -> "RXParams":
        """
        Read information from toml file

        Parameters
        ----------
        f: str
            The name of input toml File

        Returns
        -------
        DFTParams: DFTParams object
            self
        """
        import toml

        return cls.from_dict(toml.load(fname))
