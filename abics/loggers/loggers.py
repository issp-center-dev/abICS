"""Logger initialization for package."""

import logging
import os, sys
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from pathlib import Path

    from mpi4py import MPI

    _MPI_APPEND_MODE = MPI.MODE_CREATE | MPI.MODE_APPEND | MPI.MODE_WRONLY

logging.getLogger(__name__)

__all__ = ["set_log_handles"]

# logger formater
FFORMATTER = logging.Formatter(
    "[%(asctime)s] %(app_name)s %(levelname)s %(name)s %(message)s"
)
CFORMATTER = logging.Formatter(
    "%(app_name)s %(levelname)s %(message)s"
)
FFORMATTER_MPI = logging.Formatter(
    "[%(asctime)s] %(app_name)s rank:%(rank)s %(levelname)s %(name)s %(message)s"
)
CFORMATTER_MPI = logging.Formatter(
    "%(app_name)s rank:%(rank)s %(levelname)s %(message)s"
)


class _AppFilter(logging.Filter):
    """Add field `app_name` to log messages."""

    def __init__(self, app_name: str) -> None:
        super().__init__(name="App_name")
        self.app_name = app_name

    def filter(self, record):
        record.app_name = self.app_name
        return True


class _MPIRankFilter(logging.Filter):
    """Add MPI rank number to log messages, adds field `rank`."""

    def __init__(self, rank: int) -> None:
        super().__init__(name="MPI_rank_id")
        self.mpi_rank = str(rank)

    def filter(self, record):
        record.rank = self.mpi_rank
        return True


class _MPIMasterFilter(logging.Filter):
    """Filter that lets through only messages emited from rank==0."""

    def __init__(self, rank: int) -> None:
        super().__init__(name="MPI_master_log")
        self.mpi_rank = rank

    def filter(self, record):
        record.rank = self.mpi_rank
        return self.mpi_rank == 0

class _MPIRankSelectFilter(logging.Filter):
    """Filter that lets through only messages emited from rank if flag is set."""

    def __init__(self, rank: int, output_flag: bool) -> None:
        super().__init__(name="MPI_rank_log")
        self.mpi_rank = rank
        self.output_flag = output_flag

    def filter(self, record):
        record.rank = self.mpi_rank
        return self.output_flag

class _MPIFileStream:
    """Wrap MPI.File` so it has the same API as python file streams.

    Parameters
    ----------
    filename : Path
        disk location of the file stream
    MPI : MPI
        MPI communicator object
    mode : str, optional
        file write mode, by default _MPI_APPEND_MODE
    """

    def __init__(
        self, filename: "Path", MPI: "MPI", mode: int = None
    ) -> None:
        if mode is None:
            mode = MPI.MODE_CREATE | MPI.MODE_APPEND | MPI.MODE_WRONLY
        self.stream = MPI.File.Open(MPI.COMM_WORLD, filename, mode)
        self.stream.Set_atomicity(True)
        self.name = "MPIfilestream"

    def write(self, msg: str):
        """Write to MPI shared file stream.

        Parameters
        ----------
        msg : str
            message to write
        """
        b = bytearray()
        b.extend(map(ord, msg))
        self.stream.Write_shared(b)

    def close(self):
        """Synchronize and close MPI file stream."""
        self.stream.Sync()
        self.stream.Close()


class _MPIHandler(logging.FileHandler):
    """Emulate `logging.FileHandler` with MPI shared File that all ranks can write to.

    Parameters
    ----------
    filename : Path
        file path
    MPI : MPI
        MPI communicator object
    mode : str, optional
        file access mode, by default "_MPI_APPEND_MODE"
    """

    def __init__(
        self,
        filename: "Path",
        MPI: "MPI",
        mode: int = None,
    ) -> None:
        self.MPI = MPI
        self.mpi_mode = mode
        super().__init__(filename, mode="b", encoding=None, delay=False)

    def _open(self):
        return _MPIFileStream(self.baseFilename, self.MPI, self.mpi_mode)

    def setStream(self, stream):
        """Stream canot be reasigned in MPI mode."""
        raise NotImplementedError("Unable to do for MPI file handler!")


def set_log_handles(
    app_name: str,
    level: int = logging.INFO,
    console: str = "default",
    console_level: int = None,
    logfile_path: Optional["Path"] = None,
    logfile_mode: Optional[str] = None,
    logfile_level: int = None,
    logfile_rank = None,
    params: Optional['dict'] = None,
):
    """Set desired level for package loggers and add file handlers.

    Parameters
    ----------
    app_name: str
        application name
    level: int
        logging level
    log_path: Optional[str]
        path to log file, if None logs will be send only to console. If the parent
        directory does not exist it will be automatically created, by default None
    mpi_log: Optional[str], optional
        mpi log type. Has four options.
         - `master` will output logs to file and console only from rank==0.
         - `collect` will write messages from all ranks to one file opened under
           rank==0 and to console.
         - `workers` will open one log file for each worker designated by its rank,
           console behaviour is the same as for `collect`.
         - `rank` will write messages from ranks specified by rank parameter to
           one file similar to `collect` mode.
        If this argument is specified, package 'mpi4py' must be already installed.
        by default None
    rank: int or list of int
        rank number or list of rank numbers from which logs are written to console
        and/or file when mpi_log is set to `rank` mode.

    Raises
    ------
    RuntimeError
        If the argument `mpi_log` is specified, package `mpi4py` is not installed.

    References
    ----------
    https://groups.google.com/g/mpi4py/c/SaNzc8bdj6U
    https://stackoverflow.com/questions/35869137/avoid-tensorflow-print-on-standard-error
    https://stackoverflow.com/questions/56085015/suppress-openmp-debug-messages-when-running-tensorflow-on-cpu

    Notes
    -----
    Logging levels:

    +---------+--------------+----------------+----------------+----------------+
    |         | our notation | python logging | tensorflow cpp | OpenMP         |
    +=========+==============+================+================+================+
    | debug   | 10           | 10             | 0              | 1/on/true/yes  |
    +---------+--------------+----------------+----------------+----------------+
    | info    | 20           | 20             | 1              | 0/off/false/no |
    +---------+--------------+----------------+----------------+----------------+
    | warning | 30           | 30             | 2              | 0/off/false/no |
    +---------+--------------+----------------+----------------+----------------+
    | error   | 40           | 40             | 3              | 0/off/false/no |
    +---------+--------------+----------------+----------------+----------------+

    """
    # # silence logging for OpenMP when running on CPU if level is any other than debug
    # if level <= 10:
    #     os.environ["KMP_WARNINGS"] = "FALSE"

    # # set TF cpp internal logging level
    # os.environ['TF_CPP_MIN_LOG_LEVEL'] = str(int((level / 10) - 1))

    level_errlog = logging.ERROR

    # expand parameter pack
    if params:
        level = params.get("level", level)
        console = params.get("console", console)
        console_level = params.get("console_level", console_level)
        _logfile_path = params.get("logfile_path", None)
        if _logfile_path:
            from pathlib import Path
            logfile_path = Path(_logfile_path)
        logfile_mode = params.get("logfile_mode", logfile_mode)
        logfile_level = params.get("logfile_level", logfile_level)
        logfile_rank = params.get("logfile_rank", logfile_rank)

    # get root logger
    root_log = logging.getLogger()

    # remove all old handlers
    root_log.setLevel(level)
    for hdlr in root_log.handlers[:]:
        root_log.removeHandler(hdlr)

    # check if MPI environment is available
    MPI = None
    if console in [ 'serial', 'none' ] or logfile_path == None or logfile_mode == 'serial':
        # skip
        pass
    else:
        try:
            from mpi4py import MPI
        except ImportError as e:
            pass

    # check mode
    if console == "default":
        if MPI:
            console = "mpi"
        else:
            console = "serial"

    if console == "mpi" and not MPI:
        raise RuntimeError("You cannot specify 'mpi' for console when mpi4py not available")

    if logfile_path:
        if logfile_mode in [ 'master', 'collect', 'workers' ] and not MPI:
            raise RuntimeError("You cannot specify '{}' for logfile_mode when mpi4py not available".format(logfile_mode))

    # check level
    console_level = console_level if console_level else level
    if console != "none" and not console_level:
        raise RuntimeError("You must specify either level or console_level")

    logfile_level = logfile_level if logfile_level else level
    if logfile_path and not logfile_level:
        raise RuntimeError("You must specify either level or logfile_level")

    # check rank
    if logfile_rank:
        logfile_rank = [ logfile_rank ] if type(logfile_rank) is not list else logfile_rank

    # * add console handler ************************************************************
    if console == "mpi":
        _rank = MPI.COMM_WORLD.Get_rank()

        # - set console log handler
        ch_out = logging.StreamHandler(sys.stdout)

        ch_out.setLevel(console_level)
        ch_out.addFilter(_MPIMasterFilter(_rank))
        ch_out.setFormatter(CFORMATTER)

        ch_out.addFilter(_AppFilter(app_name))
        root_log.addHandler(ch_out)

        # - set error log handler
        ch_err = logging.StreamHandler(sys.stderr)

        ch_err.setLevel(level_errlog)
        ch_err.addFilter(_MPIRankFilter(_rank))
        ch_err.setFormatter(CFORMATTER_MPI)

        ch_err.addFilter(_AppFilter(app_name))
        root_log.addHandler(ch_err)

    elif console == "serial":
        # - set console log handler
        ch_out = logging.StreamHandler(sys.stdout)

        ch_out.setLevel(console_level)
        ch_out.setFormatter(CFORMATTER)

        ch_out.addFilter(_AppFilter(app_name))
        root_log.addHandler(ch_out)

        # - et err log handler
        ch_err = logging.StreamHandler(sys.stderr)

        ch_err.setLevel(level_errlog)
        ch_err.setFormatter(CFORMATTER)

        ch_err.addFilter(_AppFilter(app_name))
        root_log.addHandler(ch_err)

    elif console == "none":
        # - suppress console output
        pass

    else:
        raise RuntimeError("Unsupported console mode {}".format(console))

    # * add file handler ***************************************************************
    if logfile_path:
        # create directory
        logfile_path.parent.mkdir(exist_ok=True, parents=True)

        fh = None
        if logfile_mode is None or logfile_mode == "serial":
            fh = logging.FileHandler(logfile_path, mode="w")
            fh.setFormatter(FFORMATTER)
        elif logfile_mode == "master":
            _rank = MPI.COMM_WORLD.Get_rank()
            if _rank == 0:
                fh = logging.FileHandler(logfile_path, mode="w")
                fh.addFilter(_MPIMasterFilter(_rank))
                fh.setFormatter(FFORMATTER)
        elif logfile_mode == "collect":
            _rank = MPI.COMM_WORLD.Get_rank()
            fh = _MPIHandler(logfile_path, MPI, mode=MPI.MODE_WRONLY | MPI.MODE_CREATE)
            fh.addFilter(_MPIRankSelectFilter(_rank, logfile_rank is None or _rank in logfile_rank))
            fh.setFormatter(FFORMATTER_MPI)
        elif logfile_mode == "workers":
            _rank = MPI.COMM_WORLD.Get_rank()
            # if file has suffix than inser rank number before suffix
            # e.g deepmd.log -> deepmd_<rank>.log
            # if no suffix is present, insert rank as suffix
            # e.g. deepmdlog -> deepmdlog.<rank>
            if logfile_path.suffix:
                worker_log = (logfile_path.parent / f"{logfile_path.stem}_{_rank}").with_suffix(logfile_path.suffix)
            else:
                worker_log = logfile_path.with_suffix(f".{_rank}")
            if logfile_rank is None or _rank in logfile_rank:
                fh = logging.FileHandler(worker_log, mode="w")
                fh.setFormatter(FFORMATTER)
        else:
            raise RuntimeError("Unsupported logfile mode {}".format(logfile_mode))

        if fh:
            fh.setLevel(logfile_level)
            fh.addFilter(_AppFilter(app_name))
            root_log.addHandler(fh)

    # **********************************************************************************
