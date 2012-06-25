all: exe

exe:
	(cd obj; make -j3; cd ..)
	(cd bin; make -j3; cd ..)

debug:
	(cd obj; make debug -j3; cd ..)
	(cd bin; make debug -j3; cd ..)

prof:
	(cd obj; make prof -j3; cd ..)
	(cd bin; make prof -j3; cd ..)

doc:
	cd doc; doxygen Doxyfile

clean :
	(cd obj; make clean; cd ..)
	(cd bin; make clean; cd ..)
