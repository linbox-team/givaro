Givaro
======

[![Build Status](https://ci.inria.fr/linbox/buildStatus/icon?job=Givaro)](https://ci.inria.fr/linbox/job/Givaro/)

Download and install
--------------------

For lastest releases, please check out [this website](https://forge.imag.fr/frs/?group_id=187).
Then, you can install doing:

```
> tar -zxvf givaro-*.tar.gz
> cd givaro-*
> ./configure --prefix=##GIVAROROOT#
> make install
```

*Configuration can be adapted. Check `configure --help` to print the parameter choices.*

*In particular if GMP is not installed to the default location you might need to add for instance `--with-gmp=##GMPROOT#/gmp-x-y-z` to the configure line.*

*Also, on non-Linux systems you might need to use `gmake` instead of `make`.*

Compile your own files
----------------------

An optional compilation help file is provided: just add the following line to your Makefile. Then a simple call will compile your C and C++ files.

```
include ##GIVAROROOT##/bin/givaro-makefile
```

However, if you want to do it without this tool, you should add `-I##GIVAROROOT##/include` to the CXX compilation flags, and `-L##GIVAROROOT##/lib -lgivaro` to the LD link flags, along those for GMP.
