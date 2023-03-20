#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <vector>
#include <math.h>

using namespace std;

int main(){
	// Recopilación de información:
	ifstream fin;
	fin.open("Celdasx_malla.dat");
		float valor;
		vector<float> xcr;
		while (fin >> valor){
			xcr.push_back(valor);
		}
	fin.close();
	// Condiciones de frontera:
	float BC_E = 15; // [°C] Condición Dirichlet en la frontera este (Fila 1 Fronteras_ctes.dat)
	float BC_N = 0; // [°C / m] Condición Neumann en la frontera norte (Fila 2 Fronteras_ctes.dat)
	float BC_W = -5; // [°C] Condición Dirichlet en la frontera oeste (Fila 3 Fronteras_ctes.dat)
	vector<float> BC_S(xcr.size()); // [°C] Condición Dirichlet en la frontera sur
	for (int i = 0; i < BC_S.size(); i++){
		BC_S[i] = -5 + 20 * pow((xcr[i] / xcr[xcr.size() - 1]), 2) + 20 * sin(4 * M_PI * xcr[i] / xcr[xcr.size() - 1]);
	}
	// Escritura de la información:
	ofstream fout;
	fout.open("Fronteras_ctes.dat");
		fout << BC_E << endl;
		fout << BC_N << endl;
		fout << BC_W << endl;
	fout.close();
	fout.open("Fronteras_vrbls.dat");
		for (int i = 0; i < BC_S.size(); i++){
			fout << BC_S[i] <<endl;
		}
	fout.close();
	return 0;
}
