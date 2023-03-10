#include <iostream>
#include <random>
#include <fstream>

using namespace std;

int main (){
	// Generar números aleatorios con distribución uniforme:
	const int n = 50;
	random_device rd_1{};
	mt19937 gen_1{rd_1()};
	uniform_int_distribution<int> d_1(0,10);
	int Numeros[n]{};
	for (int i = 0; i < n; i++) {
    		Numeros[i] = d_1(gen_1);
    		cout << Numeros[i] << " ";
  	}
  	cout << endl;
  	// Generar números aleatorios con distribución Gaussiana:
  	random_device rd_2{};
  	mt19937 gen_2{rd_2()};
  	normal_distribution<float> d_2{-10, 17};
  	float Num_flo[n]{};
  	for (int i = 0; i < n; i++) {
    		Num_flo[i] = d_2(gen_2);
    		cout << Num_flo[i] << " ";
  	}
  	cout << endl;
  	// Condicionales:
  	int uni = 0;
  	int gauss = 0;
  	for (int i = 0; i < n; i++){
  		if (Numeros[i] >= 7 && uni == 0){
  			uni = 1;
  		}
  		if (Num_flo[i] >= 7 && gauss == 0){
  			gauss = 1;
  		}
  		if (uni == 1 && gauss == 1){
  			break;
  		}
  		else if (uni == 1 && gauss == 0 && i == (n-1)){
  			cout << "Distribución gaussiana: Todos menores que 7" << endl;
  		}
  		else if (uni == 0 && gauss == 1 && i == (n-1)){
  			cout << "Distribución uniforme: Todos menores que 7" << endl;
  		}
  		else if (uni == 0 && gauss == 0 && i == (n-1)){
  			cout << "Distribución uniforme: Todos menores que 7" << endl;
  			cout << "Distribución gaussiana: Todos menores que 7" << endl;
  		}
  	}
  	uni = 0;
  	gauss = 0;
  	for (int i = 0; i < n; i++){
  		if (Numeros[i] > 9 && uni == 0){
  			uni = 1;
  		}
  		if (Num_flo[i] > 9 && gauss == 0){
  			gauss = 1;
  		}
  		if (uni == 1 && gauss == 1){
  			cout << "Distribución uniforme: Al menos uno mayor que 9" << endl;
  			cout << "Distribución gaussiana: Al menos uno mayor que 9" << endl;
  			break;
  		}
  		else if (uni == 1 && gauss == 0 && i == (n-1)){
  			cout << "Distribución uniforme: Al menos uno mayor que 9" << endl;
  		}
  		else if (uni == 0 && gauss == 1 && i == (n-1)){
  			cout << "Distribución gaussiana: Al menos uno mayor que 9" << endl;
  		}
  	}
  	uni = 0;
  	gauss = 0;
  	for (int i = 0; i < n; i++){
  		if (Numeros[i] > 3 && uni == 0){
  			cout << "Distribución uniforme > 3: " << Numeros[i] << endl;
  			uni = 1;
  		}
  		if (Num_flo[i] > 3 && gauss == 0){
  			cout << "Distribución gaussiana > 3: " << Num_flo[i] << endl;
  			gauss = 1;
  		}
  		if (uni == 1 && gauss == 1){
  			break;
  		}
  	}
  	// Generar números aleatorios con distribución uniforme:
	const int m = 1000;
	random_device rd_3{};
	mt19937 gen_3{rd_3()};
	uniform_int_distribution<int> d_3(0,10);
	int Numeros_1[m]{};
	for (int i = 0; i < m; i++) {
    		Numeros_1[i] = d_3(gen_3);
  	}
  	// Generar números aleatorios con distribución Gaussiana:
  	random_device rd_4{};
  	mt19937 gen_4{rd_4()};
  	normal_distribution<float> d_4{-10, 17};
  	float Num_flo_1[m]{};
  	for (int i = 0; i < m; i++) {
    		Num_flo_1[i] = d_4(gen_4);
  	}
  	ofstream file;
	file.open("CortesKevin_HPC_S5C1.dat");
		file << "Dist. Uniforme" << "\t" << "Dist. Gaussiana" << endl;
		for (int i = 0; i < m; i++){
			file << Numeros_1[i] << "\t" << Num_flo_1[i] << endl;
		}
	file.close();
	return 0;
}
