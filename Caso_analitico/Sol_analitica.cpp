#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

int main(){
	// Recopilación de información:
	ifstream fin;
	fin.open("Fronteras_ctes.dat");
		double valor = 0;
		vector<double> BC;
		while (fin >> valor){
			BC.push_back(valor);
		}
		double T_max = BC[0]; // Temperatura de la frontera más caliente.
		double T_min = BC[1]; // Temperatura de la frontera menos caliente.
	fin.close();
	// Parámetros geométricos:
	double L = 3; // [m] Tamaño en x
	double H = 2; // [m] Tamaño en y
	// Dominio para la solución analítica:
	int nx = 180; // [-] Divisiones en x
	int ny = 120; // [-] Divisiones en y
	double dx = L / double(nx); // [m] Paso en x
	double dy = H / double(ny); // [m] Paso en y
	vector<double> x(nx + 1); // [m] Vector de longitud en x
	vector<double> y(ny + 1); // [m] Vector de longitud en y
	for (int i = 0; i < x.size(); i++){
		if (i == 0){
			x[i] = 0;
		}
		else{
			x[i] = x[i-1] + dx;
		}
	}
	for (int i = 0; i < y.size(); i++){
		if (i == 0){
			y[i] = 0;
		}
		else{
			y[i] = y[i-1] + dy;
		}
	}
	// Creación de la matriz de Temperatura de la solución analítica:
	vector<double> T(x.size() * y.size()); // [°C] Matriz de temperatura
	double T_sum = 0; // [-] Valor puntual de la temperatura adimensional
	int n = 75; // [-] Número de términos para la solución analítica
	for (int i = 0; i < y.size(); i++){
		for (int j = 0; j < x.size(); j++){
			int p = i * x.size() + j;
			double c = 0;
			double Term1 = 0;
			double Term2 = 0;
			for (int k = 0; k < n; k++){
				c = double(2 * k + 1);
				Term1 = (2 * sinh(c * M_PI * x[j] / H)) / (c * sinh(c * M_PI * L / H));
				Term2 = sin(c * M_PI * y[i] / H);
				T_sum = T_sum + Term1 * Term2;
			}
			T[p] = (2 / M_PI) * T_sum;
			T_sum = 0;
		}
	}
	double T_ad_max = *max_element(T.begin(), T.end());
	for (int i = 0; i < y.size(); i++){
		for (int j = 0; j < x.size(); j++){
			int p = i * x.size() + j;
			T[p] = ((T_max - T_min) / T_ad_max) * T[p] + T_min;
		}
	}
	// Escritura de la información:
	ofstream fout;
	fout.open("x_analitica.dat");
		for (int i = 0; i < x.size(); i++){
			fout << x[i] << endl;
		}
	fout.close();
	fout.open("y_analitica.dat");
		for (int i = 0; i < y.size(); i++){
			fout << y[i] << endl;
		}
	fout.close();
	fout.open("Temp_analitica.dat");
		for (int i = 0; i < y.size(); i++){
			for (int j = 0; j < x.size(); j++){
				int p = i * x.size() + j;
				fout << T[p] << "\t"; 
			}
			fout << "\n";
		}
	fout.close();
}
