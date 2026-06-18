# nec2dXS

A modified NEC-2 double-precision engine for the 4nec2 antenna modelling
package. It is based on the original NEC-2 sources by G. Burke and A. Poggio
of Lawrence Livermore Laboratory, modified and maintained by Arie Voors for use
with 4nec2.

This repository packages the engine for a unified NEC-2 validation suite. The
suite draws roughly 1800 sample decks from the
[Cebik, W4RNL](https://github.com/antenna2/cebik) corpus (silent key), together
with the example decks from [nec2c](https://github.com/KJ7LNW/nec2c) and xnec2c,
then validates results across the original Fortran engine and the independent
C++ implementation [necpp](https://github.com/tmolteno/necpp). The comparison
forms a four-way graph across implementations — `nec2dx` ↔ `nec2c` ↔ `xnec2c` ↔
`necpp` — so a discrepancy on any edge surfaces the responsible bug. The deck
corpus and harness live in
[`xnec2c-validation-data`](https://github.com/KJ7LNW/xnec2c-validation-data).

## Name lineage

The `nec2d`/`nec2dx`/`nec2dxs` names mark a line of descent:

- `nec2d` — LLNL's double-precision NEC-2; the trailing `d` denotes double
  precision.
- `nec2dx` — the base file as distributed, `x` indicating the extended
  double-precision variant.
- `nec2dxs` — Arie Voors' work beginning at `av03` (15-mar-02): merging the
  SomNec routines into `nec2dx`, yielding a single unified executable.

## Maintenance and contributing

The intent is to keep this Fortran reference buildable on modern systems, so it
can serve as a clean, maintained canonical reference engine for posterity. The
4nec2 distribution publishes its sources outside version control; tracking the
engine in git provides the full history for auditing back to the origin. The
upstream Nec2dXS sources ship from the 4nec2 support-files page
([supfiles.htm](https://www.qsl.net/4nec2/supfiles.htm),
[Nec2dXS_src.zip](https://qsl.net/4nec2/Nec2dXS_src.zip)).

Other Fortran NEC-2 ports exist on GitHub. They are linked here for third-party
reference, not to be replaced — each is largely unmaintained, partial, or
correctness-unverified at publication:

- [`Jheengut/nec2dxs`](https://github.com/Jheengut/nec2dxs) — a Fortran 90 port
  its author declares barely compiles and does not process models correctly.
- [`yeti01/nec2`](https://github.com/yeti01/nec2) — a minimal "compiles under
  gfortran" wrapper retaining the static per-size `.INC` build.
- [`Levitator1/nec2dXS`](https://github.com/Levitator1/nec2dXS) — minimal changes
  to compile under an older Debian 10 / gfortran 8.3 toolchain, correctness
  unverified.

Broad build compatibility is the target; the engine currently builds with
`gfortran`. Merge requests are welcome when they stay true to the original
Fortran intention, fixing bugs only as necessary along the way.

## Overview

The principal change from stock NEC-2 is that the SomNec ground (Sommerfeld/
Norton) calculations are integrated directly into the NEC-2 engine. No external
`som2d.nec` file is generated; the required data is held internally.

Consequences of the integrated ground model:

- On modern hardware the SomNec ground calculation is fast (typically under one
  second), so pre-generated ground data is no longer needed.
- The SomNec data can be recalculated at every step of a frequency loop or
  sweep, which earlier nec2d(x) engines could not do.
- A frequency sweep over many steps adds roughly one second of ground
  calculation per step. A 3-30 MHz sweep at 0.1 MHz increment adds about
  `27 * 10 * 1 = 270` seconds.
- The legacy negative-conductivity trick still works: when specified, the
  SomNec calculation runs once before the frequency loop, restoring the original
  nec2d precision and speeding up the run.

## Source files

| File                   | Purpose                                       |
| ---------------------- | --------------------------------------------- |
| `nec2dxs.f`            | Main NEC-2 engine source                      |
| `G77PORT.INC`          | Compiler/port identification string           |
| `v<size>/NEC2DPAR.INC` | Per-variant array-size parameters             |

The source hardcodes `INCLUDE 'NEC2DPAR.INC'` and `'G77PORT.INC'`. On a
case-sensitive filesystem each variant therefore keeps a file named exactly
`NEC2DPAR.INC` in its own directory, and the build points `-I` at the chosen
directory. Nothing is copied or linked.

### Array-size variants

Each variant lives in `src/v<size>/NEC2DPAR.INC`. `MAXSEG` and `MAXMAT` set the
segment count and in-core matrix allocation; `NSMAX`, `NETMX`, and `LOADMX`
bound the EX, NT/TL, and LD card counts. The `Matrix CM` column is the reserved
size of the interaction matrix, `MAXMAT**2 * 16` bytes (double-precision
complex); it dominates the memory footprint.

| Variant directory | MAXSEG | NSMAX | NETMX | Matrix `CM` |
| ----------------- | ------ | ----- | ----- | ----------- |
| `src/v500`        | 500    | 64    | 64    | 4 MB        |
| `src/v1k5`        | 1500   | 99    | 128   | 36 MB       |
| `src/v3k0`        | 3000   | 128   | 256   | 144 MB      |
| `src/v5k0`        | 5000   | 128   | 256   | 400 MB      |
| `src/v8k0`        | 8000   | 128   | 256   | 1.0 GB      |
| `src/v11k`        | 11000  | 256   | 256   | 1.94 GB     |
| `src/v45k3`       | 45300  | 256   | 256   | 32.8 GB     |

### Memory model and size limits

`MAXMAT` drives the dominant allocation: the complex interaction matrix
`COMMON /CMB/ CM(MAXMAT**2)`, held at 16 bytes per element. Every other array
scales with `MAXSEG` and is comparatively small.

That matrix lands in the zero-initialised BSS segment, which occupies no file
space and is mapped demand-zero by the kernel. Two consequences follow:

- The on-disk binary is the same length for every variant; only the reserved
  virtual size differs. `v500` and `v45k3` are byte-for-byte identical in size
  on disk.
- Resident memory tracks the portion of the matrix actually written — roughly
  `16 * N**2` for an `N`-segment deck — not `MAXMAT**2`. A small deck on a large
  variant commits only the pages it touches, so there is no resident penalty for
  building a larger variant than a deck needs.

Two ceilings bound `MAXMAT`:

- Static-data addressing. Under the default x86-64 small code model all static
  data must sit below 2 GB so 32-bit relocations reach it. `CM` plus the trailing
  common blocks crosses that boundary near `MAXMAT = 11585`, which is why `11000`
  is the largest stock variant. `-mcmodel=medium` routes the matrix to a
  large-data section addressed by 64-bit relocations, lifting the wall; the
  Makefile sets it for all variants, which is what lets `v45k3` link.
- Integer index width. `IRESRV = MAXMAT**2` is a default 4-byte integer, so
  `MAXMAT**2` must not exceed `2147483647`. The hard ceiling is therefore
  `MAXMAT = 46340` (`46340**2 = 2147395600`); `46341` fails to compile. Going
  beyond requires `-fdefault-integer-8`, which widens every default integer and
  changes unformatted-record layout — the saved matrix files read by `GFIL`
  become incompatible with 4-byte-integer builds.

Past `46340` the limit is physical memory and the `O(N**3)` solve cost rather
than any fixed parameter.

## Building

Run `make` to build every variant into `bin/`:

```sh
make                  # all variants
make bin/nec2dxs1k5   # a single variant
```

The Makefile invokes `gfortran` with `-std=legacy -w -O0 -ffp-contract=off
-fno-automatic -mcmodel=medium`, adding `-Isrc/v<size>` so `INCLUDE
'NEC2DPAR.INC'` resolves to that variant. The `-w` flag suppresses the legacy
"might be used uninitialized" warnings noted in `_Compile.txt`. `-mcmodel=medium`
lets the larger variants place their multi-gigabyte matrix in static storage
without overflowing 32-bit relocations; see Memory model and size limits.

To add a custom size, create `src/v<size>/NEC2DPAR.INC` and append `<size>` to
the `VARIANTS` list in the Makefile.

## Running

Each binary accepts non-interactive command-line arguments, falling back to the
interactive filename prompts when none are given:

```sh
bin/nec2dxs500 -i input.nec -o output.txt                 # input and output
bin/nec2dxs500 -i input.nec -o output.txt -p plot.plt     # also write plot data
```

| Flag | Argument   | Purpose                                              |
| ---- | ---------- | ---------------------------------------------------- |
| `-i` | input deck | NEC-2 input file (required)                          |
| `-o` | output     | computed output listing (required)                   |
| `-p` | plot file  | plot-data destination; without it a PL card in the deck writes to `PLTDAT.NEC` |

Run with no arguments to be prompted for the input and output filenames.

## Historical Windows packages

The original distribution shipped two Windows builds compiled with the G77
version 3.2 ports:

### nec2dXS_VM (MinGW port)

Executables that use virtual memory automatically when on-board RAM is
insufficient for the segment count, via standard Windows DLLs. Smaller
executables (about 333 KB) with slightly lower performance. Most builds (500 to
11k) run on common systems; with 64 MB or less the 8k or 11k builds may fail to
start.

### nec2dXS_FB (DJGPP port)

Executables optimised for speed, especially at lower segment counts, using
static DJGPP libraries. Larger executables (about 418 KB). These cannot use
virtual memory; when physical RAM is insufficient they report
`Load error, no DPMI memory`. As of February 2006, the DOS nature of these
builds means they do not always run on later Windows-XP systems.

## Performance

Sample timings (seconds) comparing the stock `nec2d960` engine with `nec2dXS1k5`:

- Test A — single run: `4x20ra.nec`, far-field pattern at 5-degree resolution.
- Test B — optimization: `36dip.nec`, hill-climbing resonance optimization over
  27 steps.

| Test | W-95 nec2d | W-95 nec2dXS | W-Me nec2d | W-Me nec2dXS | W-XP nec2dXS |
| ---- | ---------- | ------------ | ---------- | ------------ | ------------ |
| A    | 30.98      | 11.37        | -          | 13.1         | 5.82         |
| B    | 23.9       | 28.02        | 8.7        | 8.6          | 19.5         |

## Change history

Engine modifications are tagged `avNN` in the source headers and `_Builds.txt`:

- `av00` (2002-03-01) First G77 compile for Windows.
- `av03` (2002-03-15) SomNec routines integrated into the engine.
- `av05`-`av07` (2002-10-21) Increased LD, NT, and EX card limits.
- `av010` (2003-01-30) DJGPP port, 30-60% speed gain at low segment counts.
- `av012` (2003-09-29) User-specified NGF file name.
- `av013` (2003-09-29) MinGW port for 11k segments and virtual memory.
- `av016` (2006-11-09) Official NEC-2 bugfix by J. Burke.
- `av018` (2008-10-10) Corrected default NGF name handling from `av015`.

See `_Builds.txt` for the full list and `_Compile.txt` for compilation notes.

## Credits

- Original NEC-2 engine: G. Burke and A. Poggio, Lawrence Livermore Laboratory.
- Original packager of the 4nec2 NEC-2 package, and author of the nec2dXS
  modifications: Arie Voors, PA2B (`4nec2@gmx.net`, per `src/nec2dxs.f:155`).
- Initial G77 compiler guidance: Raymond Anderson.
