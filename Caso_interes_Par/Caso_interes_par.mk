Malla_par.pdf Dist_T_par.pdf Residuales_par.pdf : Temperatura.dat Residuo_norm.dat Posprocesamiento.py
	python3 Posprocesamiento.py
Temperatura.dat Residuo_norm.dat : Procesamiento.exe
	for thr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; do \
		OMP_NUM_THREADS=$$thr ./Procesamiento.exe $$thr; \
	done
Procesamiento.exe : Fronteras_ctes.dat Fronteras_vrbls.dat Propiedades.dat Term_src.dat Procesamiento.cpp
	g++ -std=c++17 -fopenmp -fsanitize=address,undefined Procesamiento.cpp -o Procesamiento.exe
Fronteras_ctes.dat Fronteras_vrbls.dat : Frontera_t0.exe
	./Frontera_t0.exe
Frontera_t0.exe : Carasx_malla.dat Carasy_malla.dat Celdasx_malla.dat Celdasy_malla.dat Frontera_t0.cpp 
	g++ -std=c++17 -fsanitize=address,undefined Frontera_t0.cpp -o Frontera_t0.exe
Propiedades.dat Term_src.dat : Propiedades.exe
	./Propiedades.exe 
Propiedades.exe : Carasx_malla.dat Carasy_malla.dat Celdasx_malla.dat Celdasy_malla.dat Propiedades.cpp	
	g++ -std=c++17 -fsanitize=address,undefined Propiedades.cpp -o Propiedades.exe
Carasx_malla.dat Carasy_malla.dat Celdasx_malla.dat Celdasy_malla.dat : Malla.exe
	for thr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; do \
		OMP_NUM_THREADS=$$thr ./Malla.exe $$thr; \
	done
Malla.exe : Malla.cpp
	g++ -std=c++17 -fopenmp -fsanitize=address,undefined Malla.cpp -o Malla.exe
