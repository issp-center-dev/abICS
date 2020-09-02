# abICS
ab-Initio Configuration Sampling tool kit (abICS) is a Python code (library and application) for Metropolice Monte Carlo.

## Requirement

- python3 (>=3.6)
- numpy
- scipy
- toml (for parsing input files)
- mpi4py (for parallel tempering)
  - This requires one of the MPI implementation
- pymatgen (>=2019.12.3) (for using Structure as a configuration)
  - This requires Cython
- qe-tools (for parsing QE I/O)

## Install abICS

Pymatgen requires Cython but Cython will not be installed automatically,
please make sure that this is installed,

``` bash
$ pip3 install Cython
```

mpi4py requires one of the MPI implementations such as OpenMPI,
please make sure that this is also installed.
In the case of using homebrew on macOS, for example,

``` bash
$ brew install open-mpi
```

After installing Cython and MPI,

``` bash
$ pip3 install abics
```

will install abICS and dependencies.

If you want to change the directory where installed,
add `--user` option or `--prefix=DIRECTORY` option into the above command as

``` bash
$ pip3 install --user abics
```

For details of `pip` , see the manual of `pip` by `pip3 help install`

If you want to install abICS from source, see [wiki page](https://github.com/issp-center-dev/abICS/wiki/Install)

## License

The distribution of the program package and the source codes follow GNU General Public License version 3 ([GPL v3](http://www.gnu.org/licenses/gpl-3.0.en.html)). 

## Official page

https://www.pasums.issp.u-tokyo.ac.jp/abics

## Author

Shusuke Kasamatsu, Yuichi Motoyama, Kazuyoshi Yoshimi

## Manual

[English online manual](https://issp-center-dev.github.io/abICS/docs/master/en/html/index.html)

[Japnese online manual](https://issp-center-dev.github.io/abICS/docs/master/ja/html/index.html)

[API reference](https://issp-center-dev.github.io/abICS/docs/api/master/html/index.html)
