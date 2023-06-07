#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <vector>

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
	// Propiedades del material:
	vector<double> kappa(xcr.size());	// [W/(m*s)]
	for (int i = 0; i < kappa.size(); i++){
		kappa[i] = 1;
	}
	// Término fuente:
	double B = 0;
	// Escritura de la información:
	ofstream fout;
	fout.open("Propiedades.dat");
		for (int i = 0; i < kappa.size(); i++){
			fout << kappa[i] << endl;
		}
	fout.close();
	fout.open("Term_src.dat");
		fout << B << endl;
	fout.close();
	return 0;
}
