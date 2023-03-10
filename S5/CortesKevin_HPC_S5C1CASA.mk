CortesKevin_HPC_S5C1CASA.pdf : CortesKevin_HPC_S5C1CASA.tex histogramas.pdf
	pdflatex CortesKevin_HPC_S5C1CASA.tex
histogramas.pdf : CortesKevin_HPC_S5C1CASA.py CortesKevin_HPC_S5C1.dat
	python3 CortesKevin_HPC_S5C1CASA.py
CortesKevin_HPC_S5C1.dat : CortesKevin_HPC_S5C1CASA.exe
	./CortesKevin_HPC_S5C1CASA.exe
CortesKevin_HPC_S5C1CASA.exe : CortesKevin_HPC_S5C1CASA.cpp
	g++ CortesKevin_HPC_S5C1CASA.cpp -o CortesKevin_HPC_S5C1CASA.exe
