# LIBCOULCC

Calculate complex-valued Coulomb and Bessel functions for real-valued arguments.
The original `COULCC()` has been repackaged here as an
[fpm](https://github.com/fortran-lang/fpm) package.
See the original author's article [[1]](#1) for more detail.
This is largely not my work.
No intentional changes to the algorithm have been made.
Test suite not yet implemented.
Untested.

## Building with [fpm](https://github.com/fortran-lang/fpm)
In the package directory, just run

    $ fpm build --profile release

The archive file `liblibcoulcc.a` and `.mod` files will be placed in the generated `build` subdirectory.
These files are necessary to compile another program that uses this library.

## Using this version of `COULCC()`
The module `libcoulcc` contains the following public procedures :

- `coulf(l,eta,x)` : return the regular Coulomb function $F_l(\eta,x)$
- `coulg(l,eta,x)` : return the irregular Coulomb function $G_l(\eta,x)$
- `coulcc_wrapper(zlmin, nl, eta, xx, f, fp, g, gp, sig, kfn, mode)` : a wrapper used to call `COULCC()`, used by the above-mentioned procedures
- `COULCC(X,ETA1,ZLMIN,NL,FC,GC,FCP,GCP,SIG,MODE1,KFN,IFAIL)` : the original code described in [[1]](#1), with some minor modernizations.

Above, `nl` is the number of $l$ values.
The variables `l`, `x`, and `eta` are `complex(real64)`.
The kind `real64` (64 bits / 8 bytes) is defined in the intrinsic module `iso_fortran_env`.
All public procedures other than `COULCC()` are superfluous — they're provided for convenience or as a simple example of calling `COULCC()`.

The following example program

    program test

      use libcoulcc, only: coulf
      use iso_fortran_env, only: rp => real64, stdout => output_unit

      complex(rp) :: l = 0
      complex(rp) :: eta = (-0.5_rp, 0.0_rp)
      complex(rp) :: x = (20.0_rp, 0.0_rp)

      write(stdout, '("F_[", 2(F0.2, X), "i](", 2(F0.2, X), "i,", 2(F0.2, X), "i) = ", &
        & 2(e0.15, X))') l, eta, x, coulf(l, eta, x)

    end program test

should print the following:

    F_[.00 .00 i](-.50 .00 i,20.00 .00 i) = -0.102372301807428 0.000000000000000

## Reference(s)

<a id="1">[1]</a>
Thompson, I. J., and A. R. Barnett.
*COULCC: A continued-fraction algorithm for Coulomb functions of complex order with complex arguments.*
Computer physics communications 36.4 (1985): 363-372.
URL: [http://www.fresco.org.uk/papers/Coulcc-CPC.pdf](http://www.fresco.org.uk/papers/Coulcc-CPC.pdf)
