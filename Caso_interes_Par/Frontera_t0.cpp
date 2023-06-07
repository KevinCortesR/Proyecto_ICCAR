#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <vector>
#include <math.h>

int main(){
	// Recopilación de información:
        std::ifstream fin;
	fin.open("Celdasx_malla.dat");
		double valor;
		std::vector<double> xcr;
		while (fin >> valor){
			xcr.push_back(valor);
		}
	fin.close();
	// Condiciones de frontera:
	double BC_E = 15; // [°C] Condición Dirichlet en la frontera este (Fila 1 Fronteras_ctes.dat)
	double BC_N = 0; // [°C / m] Condición Neumann en la frontera norte (Fila 2 Fronteras_ctes.dat)
	double BC_W = -5; // [°C] Condición Dirichlet en la frontera oeste (Fila 3 Fronteras_ctes.dat)
	std::vector<double> BC_S(xcr.size()); // [°C] Condición Dirichlet en la frontera sur
	for (int i = 0; i < BC_S.size(); i++){
		BC_S[i] = -5 + 20 * pow((xcr[i] / xcr[xcr.size() - 1]), 2) + 20 * sin(4 * M_PI * xcr[i] / xcr[xcr.size() - 1]);
	}
	// Escritura de la información:
	std::ofstream fout;
	fout.open("Fronteras_ctes.dat");
		fout << BC_E << "\n";
		fout << BC_N << "\n";
		fout << BC_W << "\n";
	fout.close();
	fout.open("Fronteras_vrbls.dat");
		for (int i = 0; i < BC_S.size(); i++){
			fout << BC_S[i] << "\n";
		}
	fout.close();
	return 0;
}
