# abICS
ab-Initio Configuration Sampling tool kit (abICS) is a Python code (library and application) for Metropolice Monte Carlo.

## Requirement

- python3 (>=3.6)
- numpy
- scipy
- toml (for parsing input files)
- mpi4py (for parallel tempering)
- pymatgen (for parsing vasp I/O)
- qe-tools (for parsing QE I/O)

## Install abICS

``` bash
$ pip3 install abics
```

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

Under construction

## Author

Shusuke Kasamatsu, Yuichi Motoyama, Kazuyoshi Yoshimi

## Manual

[English online manual](https://issp-center-dev.github.io/abICS/docs/sphinx/en/build/html/index.html)

[Japnese online manual](https://issp-center-dev.github.io/abICS/docs/sphinx/ja/build/html/index.html)
