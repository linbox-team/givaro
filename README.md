Givaro
======

CI Inria: [![Build Status](https://ci.inria.fr/linbox/buildStatus/icon?job=Givaro)](https://ci.inria.fr/linbox/job/Givaro/)

Travis: [![Build Status](https://travis-ci.org/linbox-team/givaro.svg?branch=master)](https://travis-ci.org/linbox-team/givaro/)

Download and install
--------------------

For lastest releases, please check out [this website](http://github.com/linbox-team/givaro); older releases can be found on [that website](https://forge.imag.fr/frs/?group_id=187).
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

Givaro uses pkgconfig to expose the compilation flags it requires.

You will get the compilation flags by calling 
```
pkg-config --cflags givaro
``` 
and the linking flags by calling 
```
pkg-config --libs givaro
```

If you have installed givaro in a non-standard directory (such as `/usr/local`), make sure to have added the path where to find givaro's .pc file to the PKG_CONFIG_PATH environment variable.
```
PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:<path to your givaro install>/lib/pkgconfig
```

An alternative option is to just add the following line to your Makefile. Then a simple call will compile your C and C++ files.
```
include ##GIVAROROOT##/bin/givaro-makefile
```
