FC = mpif90-openmpi-mp
#FFLAGS = -O3 -fdefault-real-8 -fdefault-double-8 -Wall -finit-real=NaN -finit-integer=-2147483648 -g -fbacktrace -fimplicit-none -fcheck=all -ffpe-trap=invalid,zero,overflow,denormal 
FFLAGS = -O3 -fdefault-real-8 -fdefault-double-8 
CODEHOME = ~/Dropbox/papers/UltimaThule/yoshikozai

main: $(CODEHOME)/yoshikozai.f90
	$(FC) $(FFLAGS) $(CODEHOME)/yoshikozai.f90 -o yoshikozai.x

clean:
	rm -f yoshikozai.x