Malla.eps Dist_T.eps Dist_T_An Residuales.eps : Temperatura.dat Residuo_norm.dat x_analitica.dat y_analitica.dat Temp_analitica.dat Posprocesamiento.py
	python3 Posprocesamiento.py
x_analitica.dat y_analitica.dat Temp_analitica.dat : Sol_analitica.exe
	./Sol_analitica.exe
Sol_analitica.exe : Sol_analitica.cpp
	g++ Sol_analitica.cpp -o Sol_analitica.exe
Temperatura.dat Residuo_norm.dat : Procesamiento.exe
	./Procesamiento.exe
Procesamiento.exe : Fronteras_ctes.dat Propiedades.dat Term_src.dat Procesamiento.cpp
	g++ Procesamiento.cpp -o Procesamiento.exe
Fronteras_ctes.dat : Frontera_t0.exe
	./Frontera_t0.exe
Frontera_t0.exe : Carasx_malla.dat Carasy_malla.dat Celdasx_malla.dat Celdasy_malla.dat Frontera_t0.cpp 
	g++ Frontera_t0.cpp -o Frontera_t0.exe
Propiedades.dat Term_src.dat : Propiedades.exe
	./Propiedades.exe 
Propiedades.exe : Carasx_malla.dat Carasy_malla.dat Celdasx_malla.dat Celdasy_malla.dat Propiedades.cpp	
	g++ Propiedades.cpp -o Propiedades.exe
Carasx_malla.dat Carasy_malla.dat Celdasx_malla.dat Celdasy_malla.dat : Malla.exe
	./Malla.exe
Malla.exe : Malla.cpp
	g++ Malla.cpp -o Malla.exe
