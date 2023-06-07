#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <vector>

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
	// Propiedades del material:
	std::vector<double> kappa(xcr.size());	// [W/(m*s)]
	for (int i = 0; i < kappa.size(); i++){
		kappa[i] = 2 * (2.5 - xcr[i] / xcr[xcr.size() - 1]);
	}
	// Término fuente:
	double B = 100;
	// Escritura de la información:
	std::ofstream fout;
	fout.open("Propiedades.dat");
		for (int i = 0; i < kappa.size(); i++){
			fout << kappa[i] << "\n";
		}
	fout.close();
	fout.open("Term_src.dat");
		fout << B << "\n";
	fout.close();
	return 0;
}
