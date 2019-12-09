# abICS
ab-Initio Configuration Sampling tool kit (abICS) is a Python code (library and application) for Metropolice Monte Carlo.

## Requirement

- python3
- numpy
- scipy
- toml (for parsing input files)
- mpi4py (for parallel tempering)
- pymatgen (for parsing vasp I/O)
- qe-tools (for parsing QE I/O)

## Get ABICS

### With Git (recommended)

``` bash
$ git clone https://github.com/issp-center-dev/abics
```

### Without Git

``` bash
$ wget https://github.com/issp-center-dev/abics/archive/master.zip
$ unzip master.zip
```

## Install abICS

### With `pip`

Make a wheel file by

``` bash
$ ./make_wheel.sh
```

and install this as

``` bash
$ pip install dist/abics-*.whl
```

If you want to change the directory where installed,
add `--user` option or `--prefix=DIRECTORY` option into the above command as

``` bash
$ pip install --user dist/abics-*.whl
```

For details, see the manual of `pip` by `pip help install`

### Without `pip`

``` bash
$ export PYTHONPATH=$Path_To_abics:$PYTHONPATH
```

## Tutorial

See [py_mc's wiki](https://github.com/skasamatsu/py_mc/wiki/Getting-Started) (temporary).

## License

The distribution of the program package and the source codes follow GNU General Public License version 3 ([GPL v3](http://www.gnu.org/licenses/gpl-3.0.en.html)). 

## Official page

Under construction

## Author

Shusuke Kasamatsu, Yuichi Motoyama, Kazuyoshi Yoshimi

## Manual

[English online manual](https://issp-center-dev.github.io/abICS/docs/sphinx/en/build/html/index.html)

[Japnese online manual](https://issp-center-dev.github.io/abICS/docs/sphinx/ja/build/html/index.html)
