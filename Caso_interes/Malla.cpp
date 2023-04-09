#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

// Para realizar un enmallado uniforme haga Delta = Delta_mean.
vector<double> refinamiento(double lim_inf, double lim_sup, vector<double> &p_int, double Delta_mean, double Delta){
	vector<double> x_ref(1);
	int nmc = ceil((lim_sup - lim_inf) / Delta_mean);
	int cont = 1;
	x_ref[0] = lim_inf;
	for (int i = 0; i < p_int.size(); i++){
		if (p_int[i] >= (Delta / 2)){
			p_int[i] = p_int[i] - Delta / 2;
		}
	}
	while (x_ref[cont - 1] < lim_sup){
		x_ref.push_back(0);
		vector<double> v_dist(p_int.size());
		for (int i = 0; i < p_int.size(); i++){
			v_dist[i] = abs(p_int[i] - x_ref[cont - 1]); 
		}
		double dist = *min_element(v_dist.begin(), v_dist.end());
		int dist_cell = ceil(dist / Delta);
		if (dist_cell > nmc){
			x_ref[cont] = x_ref[cont - 1] + Delta_mean;
			cont++;
		}
		else if (dist_cell <= nmc || dist_cell > 2){
			double Delta_cell = (Delta_mean - Delta) / (nmc) * dist_cell + Delta;
			x_ref[cont] = x_ref[cont - 1] + Delta_cell;
      			cont++;
		}
		else if (dist_cell <= 2){
			x_ref[cont] = x_ref[cont - 1] + Delta;
      			cont++;
		}
	}
	double fin = round(x_ref[x_ref.size() - 1] * 1e6) / 1e6;
	double antefin = round(x_ref[x_ref.size() - 2] * 1e6) / 1e6;
	if (((fin - antefin) < Delta / 2) || antefin >= lim_sup){
		x_ref.pop_back();
	}
	x_ref[x_ref.size() - 1] = lim_sup;
	return x_ref;
}

int main() {
	// Parámetros geométricos:
	double L = 3; // [m] Tamaño en x
	double H = 2; // [m] Tamaño en y
	// Parámetros de la malla:
	int ncx = 60; // [-] Número de celdas en x
	int ncy = 40; // [-] Número de celdas en y
	double Deltax_mean = L / double(ncx); // [m] Tamaño promedio de las celdas en x
	double Deltay_mean = H / double(ncy); // [m] Tamaño promedio de las celdas en y
	double Deltax = 0.25 * Deltax_mean; // [m] Tamaño mínimo de las celdas en x
	double Deltay = 0.25 * Deltay_mean; // [m] Tamaño mínimo de las celdas en y
	int nfx = ncx + 1; // [-] Número de caras en x
	int nfy = ncy + 1; // [-] Número de caras en y
	vector<double> px_int = {0, L}; // [m] Puntos de interés para refinar en dirección x
	vector<double> py_int = {0, H}; // [m] Puntos de interés para refinar en dirección y
	// Creación de la malla:
	vector<double> xfr = refinamiento(0, L, px_int, Deltax_mean, Deltax); // [m] Vector de caras en x con refinamiento
	vector<double> yfr = refinamiento(0, H, py_int, Deltay_mean, Deltay); // [m] Vector de caras en y con refinamiento
	vector<double> xcr(xfr.size() - 1); // [m] Vector de celdas (centroides) en x con refinamiento
	vector<double> ycr(yfr.size() - 1); // [m] Vector de celdas (centroides) en y con refinamiento
	for (int i = 0; i < xcr.size(); i++){
		xcr[i] = 0.5 * (xfr[i] + xfr[i + 1]);
	}
	for (int i = 0; i < ycr.size(); i++){
		ycr[i] = 0.5 * (yfr[i] + yfr[i + 1]);
	}
	// Escritura de la información:
	ofstream file;
	file.open("Carasx_malla.dat");
		for (int i = 0; i < xfr.size(); i++){
			file << xfr[i] << endl;
		}
	file.close();
	file.open("Carasy_malla.dat");
		for (int i = 0; i < yfr.size(); i++){
			file << yfr[i] << endl;
		}
	file.close();
	file.open("Celdasx_malla.dat");
		for (int i = 0; i < xcr.size(); i++){
			file << xcr[i] << endl;
		}
	file.close();
	file.open("Celdasy_malla.dat");
		for (int i = 0; i < ycr.size(); i++){
			file << ycr[i] << endl;
		}
	file.close();
	return 0;
}
