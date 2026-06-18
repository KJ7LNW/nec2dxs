# Build the nec2dXS NEC-2 engine on Linux for every segment-size variant.
#
# The Fortran source includes 'NEC2DPAR.INC' and 'G77PORT.INC' verbatim, so on
# a case-sensitive filesystem files of exactly those names must sit on the
# include search path. Each variant keeps its parameters in its own directory
# as NEC2DPAR.INC; the shared G77PORT.INC sits beside the source. No include is
# copied or linked: -I aims the compiler at the directory where each file
# already lives. Selecting a variant means choosing which directory -I names.

FC     := gfortran
FFLAGS := -std=legacy -w -O0 -ffp-contract=off -fno-automatic -mcmodel=medium

SRCDIR := src
BINDIR := bin
SRC    := $(SRCDIR)/nec2dxs.f

# Binary suffix == variant directory suffix (src/v<suffix>/NEC2DPAR.INC).
VARIANTS := 500 1k5 3k0 5k0 8k0 11k 45k3
BINS     := $(addprefix $(BINDIR)/nec2dxs,$(VARIANTS))

.PHONY: all clean
all: $(BINS)

$(BINS): $(BINDIR)/nec2dxs%: $(SRC) $(SRCDIR)/v%/NEC2DPAR.INC $(SRCDIR)/G77PORT.INC | $(BINDIR)
	$(FC) $(FFLAGS) -I$(SRCDIR)/v$* -I$(SRCDIR) $(SRC) -o $@

$(BINDIR):
	mkdir -p $@

clean:
	$(RM) -r $(BINDIR)
