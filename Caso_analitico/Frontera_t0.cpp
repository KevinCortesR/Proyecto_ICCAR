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
		double valor;
		vector<double> xcr;
		while (fin >> valor){
			xcr.push_back(valor);
		}
	fin.close();
	// Condiciones de frontera:
	double BC_E = 150; // [°C] Condición Dirichlet en la frontera este (Fila 1 Fronteras_ctes.dat)
	double BC_N = 15; // [°C] Condición Dirichlet en la frontera norte (Fila 2 Fronteras_ctes.dat)
	double BC_W = 15; // [°C] Condición Dirichlet en la frontera oeste (Fila 3 Fronteras_ctes.dat)
	double BC_S = 15; // [°C] Condición Dirichlet en la frontera sur (Fila 4 Frotneras_ctes.dat)
	// Escritura de la información:
	ofstream fout;
	fout.open("Fronteras_ctes.dat");
		fout << BC_E << endl;
		fout << BC_N << endl;
		fout << BC_W << endl;
		fout << BC_S << endl;
	fout.close();
	return 0;
}
