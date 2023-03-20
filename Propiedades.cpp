#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>

using namespace std;

int main(){
	// Propiedades del material:
	float kappa = 1;	// [W/(m*s)]
	float C = 1;		// [J/(kg*m)]
	float rho = 1;		// [kg/(m³)]
	float alfa = k/(rho*C);// [m²/s]
	ofstream file;
	file.open("Propiedades.dat");
		file << "kappa [W/(m*s)]" << "\t" << "C [J/(kg*m)]" << "\t" << "rho [kg/(m³)]" << "\t" << "alfa [m²/s]" << endl;
		file << kappa << "\t" << C << "\t" << rho << "\t" << alfa << endl;
	file.close();
	
	return 0;
}
