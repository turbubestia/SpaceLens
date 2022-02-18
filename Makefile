.PHONY: qmake release

qmake:
	cd ./build/eris; qmake -makefile -o Makefile ../../eris.pro
	cd ./build/ceres; qmake -makefile -o Makefile ../../ceres.pro
	cd ./build/makemake; qmake -makefile -o Makefile ../../makemake.pro
	cd ./build/haumea; qmake -makefile -o Makefile ../../haumea.pro
	
eris-release:
	cd ./build/eris; mingw32-make release
	
eris-clean:
	cd ./build/eris; mingw32-make clean	
	
	
ceres-release:
	cd ./build/ceres; mingw32-make release
	
ceres-clean:
	cd ./build/ceres; mingw32-make clean	
	
	
makemake-release:
	cd ./build/makemake; mingw32-make release

makemake-debug:
	cd ./build/makemake; mingw32-make debug	
	
makemake-clean:
	cd ./build/makemake; mingw32-make clean	
	
haumea-release:
	cd ./build/haumea; mingw32-make release
	
haumea-clean:
	cd ./build/haumea; mingw32-make clean	