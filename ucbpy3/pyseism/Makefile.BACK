CC	= gcc
LD	= gcc

CFLAGS	= -O2 -fPIC -Isrc/include

PREFIX	= test

O	= obj

OBJS	= $O/pyseism.o \
	  $O/seism2y_ellrot_auto.o \
	  $O/loadmode4y_ellrot_auto.o \
	  $O/addmodes.o \
	  $O/geteigys_ellrot.o \
	  $O/resample.o \
	  $O/cfft.o \
	  $O/PNm.o \
	  $O/rotmatrix.o \
	  $O/rvector.o \
	  $O/rvec.o \
	  $O/svec.o \
	  $O/array.o \
	  $O/cmplx.o \
	  $O/min2n.o

# example PYTHONROOT on OS X:   /Library/Frameworks/EPD64.framework/Versions/Current (2.7)
# example PYTHONROOT on shasta: /work/shasta/sfrench/epd/current (2.7)
PYTHONROOT	= /work/shasta/sfrench/epd/current
PYTHONVERSION	= 2.7

all:	dirs pyseism_driver.py

.PHONY:	dirs
dirs:
	mkdir -p $O

pyseism_driver.py: _pyseism_driver.so

pyseism_driver_wrap.c: pyseism_driver.i pyseism_driver.h numpy.i
	swig -python pyseism_driver.i

pyseism_driver_wrap.o: pyseism_driver_wrap.c
	$(CC) -pthread -O2 -fPIC \
		-I${PYTHONROOT}/include \
		-I${PYTHONROOT}/include/python${PYTHONVERSION} \
		-I${PYTHONROOT}/lib/python${PYTHONVERSION}/site-packages/numpy/core/include \
		-c pyseism_driver_wrap.c \
		-o pyseism_driver_wrap.o

$O/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

$O/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

_pyseism_driver.so: pyseism_driver_wrap.o $(OBJS)
	$(CC) -pthread -shared \
		pyseism_driver_wrap.o $(OBJS) \
		-L${PYTHONROOT}/lib \
		-lpython${PYTHONVERSION} -lgfortran \
		-o _pyseism_driver.so

.PHONY: clean
clean:
	rm -rf $O pyseism_driver.py _pyseism_driver.so pyseism_driver_wrap.c pyseism_driver_wrap.o
	@echo
	@echo "Reminder: files possibly installed to $(PREFIX) have not been removed!"
	@echo

.PHONY: install
install: pyseism.py pyseism_driver.py _pyseism_driver.so
	cp pyseism.py pyseism_driver.py _pyseism_driver.so $(PREFIX)/

.PHONY: uninstall
uninstall:
	rm -f $(PREFIX)/pyseism.py
	rm -f $(PREFIX)/pyseism.pyc
	rm -f $(PREFIX)/pyseism_driver.py
	rm -f $(PREFIX)/pyseism_driver.pyc
	rm -f $(PREFIX)/_pyseism_driver.so
