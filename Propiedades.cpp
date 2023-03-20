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
		float valor;
		vector<float> xcr;
		while (fin >> valor){
			xcr.push_back(valor);
		}
	fin.close();
	// Propiedades del material:
	vector<float> kappa(xcr.size());	// [W/(m*s)]
	for (int i = 0; i < kappa.size(); i++){
		kappa[i] = 2 * (2.5 - xcr[i] / xcr[xcr.size() - 1]);
	}
	// Término fuente:
	float B = 100;
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
