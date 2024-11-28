CC	= gcc
FC	= gfortran
LD	= gcc

PREFIX	= ..

# example PYTHONROOT on OS X:   /Library/Frameworks/EPD64.framework/Versions/Current (2.7)
# example PYTHONROOT on shasta: /work/shasta/sfrench/epd/current (2.7)
#PYTHONROOT	=  /opt/rh/rh-python35/root/usr
PYTHONROOT      = /home/mw1685/anaconda3
PYTHONVERSION	=  3.8

all:	_pyminos_driver.so

pyminos_driver.py: pyminos_driver_wrap.c

pyminos_driver_wrap.c: pyminos_driver.i pyminos_driver.h numpy.i
	swig -python pyminos_driver.i

pyminos_driver_wrap.o: pyminos_driver_wrap.c
	$(CC) -pthread -O2 -fPIC \
		-I${PYTHONROOT}/include \
		-I${PYTHONROOT}/include/python${PYTHONVERSION} \
		-I${PYTHONROOT}/lib/python${PYTHONVERSION}/site-packages/numpy/core/include \
		-c pyminos_driver_wrap.c \
		-o pyminos_driver_wrap.o

cminos.o: cminos.c cminos.h
	$(CC) -O2 -fPIC \
		-c cminos.c \
		-o cminos.o

minos_bran_headless.o: minos_bran_headless.f
	$(FC) -O2 -fPIC \
		-c minos_bran_headless.f \
		-o minos_bran_headless.o

_pyminos_driver.so: pyminos_driver_wrap.o cminos.o minos_bran_headless.o
	$(CC) -pthread -shared \
		pyminos_driver_wrap.o cminos.o minos_bran_headless.o \
	        -L${PYTHONROOT}/lib -lpython3 -lm \
		-lgfortran \
		-o _pyminos_driver.so

.PHONY: clean
clean:
	rm -f cminos.o minos_bran_headless.o pyminos_driver.py _pyminos_driver.so pyminos_driver_wrap.c pyminos_driver_wrap.o
	@echo
	@echo "Reminder: files possibly installed to $(PREFIX) have not been removed!"
	@echo

.PHONY: install
install: pyminos.py pyminos_driver.py _pyminos_driver.so
	cp pyminos.py pyminos_driver.py _pyminos_driver.so $(PREFIX)/

.PHONY: uninstall
uninstall:
	rm -f $(PREFIX)/pyminos.py
	rm -f $(PREFIX)/pyminos_driver.py
	rm -f $(PREFIX)/_pyminos_driver.so

