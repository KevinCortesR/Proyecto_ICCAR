#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

// Función generadora de la matriz de coeficientes de los vecinos y de la celda y del vector de términos fuente:
vector<vector<float>> Coeficientes(vector<float> &xcr, vector<float> &ycr, vector<float> &xfr, vector<float> &yfr, vector<float> &kappa, float B){
	int nxcr = xcr.size();
	int nycr = ycr.size();
	vector<vector<float>> Mat_A(nxcr * nycr, vector<float> (6)); // 0: A_P, 1: A_E, 2: A_N, 3: A_W, 4: A_S, 5: Q_P
	for (int j = 0; j < nycr; j++){
		for (int i = 0; i < nxcr; i++){
			int nc = j * nxcr + i;
			if ((i > 0 && i < (nxcr - 1)) && j == (nycr - 1)){ // Celdas N (Tipo 1)
				Mat_A[nc][1] = kappa[i] * (yfr[j + 1] - yfr[j]) / (xcr[i + 1] - xcr[i]);
				Mat_A[nc][2] = 0;
				Mat_A[nc][3] = kappa[i] * (yfr[j + 1] - yfr[j]) / (xcr[i] - xcr[i - 1]);
				Mat_A[nc][4] = kappa[i] * (xfr[i + 1] - xfr[i]) / (ycr[j] - ycr[j - 1]);
				Mat_A[nc][0] = - (Mat_A[nc][1] + Mat_A[nc][2] + Mat_A[nc][3] + Mat_A[nc][4]);
				Mat_A[nc][5] = - B * (xfr[i + 1] - xfr[i]) * (yfr[j + 1] - yfr[j]);
			}
			else if ((i > 0 && i < (nxcr - 1)) && j == 0){ // Celdas S (Tipo 2)
				Mat_A[nc][1] = kappa[i] * (yfr[j + 1] - yfr[j]) / (xcr[i + 1] - xcr[i]);
				Mat_A[nc][2] = kappa[i] * (xfr[i + 1] - xfr[i]) / (ycr[j + 1] - ycr[j]);
				Mat_A[nc][3] = kappa[i] * (yfr[j + 1] - yfr[j]) / (xcr[i] - xcr[i - 1]);
				Mat_A[nc][4] = kappa[i] * (xfr[i + 1] - xfr[i]) / (ycr[j] - yfr[j]);
				Mat_A[nc][0] = - (Mat_A[nc][1] + Mat_A[nc][2] + Mat_A[nc][3] + Mat_A[nc][4]);
				Mat_A[nc][5] = - B * (xfr[i + 1] - xfr[i]) * (yfr[j + 1] - yfr[j]);
			}
			else if (i == 0 && (j > 0 && j < (nycr - 1))){ // Celdas W (Tipo 3)
				Mat_A[nc][1] = kappa[i] * (yfr[j + 1] - yfr[j]) / (xcr[i + 1] - xcr[i]);
				Mat_A[nc][2] = kappa[i] * (xfr[i + 1] - xfr[i]) / (ycr[j + 1] - ycr[j]);
				Mat_A[nc][3] = kappa[i] * (yfr[j + 1] - yfr[j]) / (xcr[i] - xfr[i]);
				Mat_A[nc][4] = kappa[i] * (xfr[i + 1] - xfr[i]) / (ycr[j] - ycr[j - 1]);
				Mat_A[nc][0] = - (Mat_A[nc][1] + Mat_A[nc][2] + Mat_A[nc][3] + Mat_A[nc][4]);
				Mat_A[nc][5] = - B * (xfr[i + 1] - xfr[i]) * (yfr[j + 1] - yfr[j]);
			}
			else if (i == (nxcr - 1) && (j > 0 && j < (nycr - 1))){ // Celdas E (Tipo 4)
				Mat_A[nc][1] = kappa[i] * (yfr[j + 1] - yfr[j]) / (xfr[i + 1] - xcr[i]);
				Mat_A[nc][2] = kappa[i] * (xfr[i + 1] - xfr[i]) / (ycr[j + 1] - ycr[j]);
				Mat_A[nc][3] = kappa[i] * (yfr[j + 1] - yfr[j]) / (xcr[i] - xcr[i - 1]);
				Mat_A[nc][4] = kappa[i] * (xfr[i + 1] - xfr[i]) / (ycr[j] - ycr[j - 1]);
				Mat_A[nc][0] = - (Mat_A[nc][1] + Mat_A[nc][2] + Mat_A[nc][3] + Mat_A[nc][4]);
				Mat_A[nc][5] = - B * (xfr[i + 1] - xfr[i]) * (yfr[j + 1] - yfr[j]);
			}
			else if (i == (nxcr - 1) && j == (nycr - 1)){ // Celda NE (Tipo 5)
				Mat_A[nc][1] = kappa[i] * (yfr[j + 1] - yfr[j]) / (xfr[i + 1] - xcr[i]);
				Mat_A[nc][2] = 0;
				Mat_A[nc][3] = kappa[i] * (yfr[j + 1] - yfr[j]) / (xcr[i] - xcr[i - 1]);
				Mat_A[nc][4] = kappa[i] * (xfr[i + 1] - xfr[i]) / (ycr[j] - ycr[j - 1]);
				Mat_A[nc][0] = - (Mat_A[nc][1] + Mat_A[nc][2] + Mat_A[nc][3] + Mat_A[nc][4]);
				Mat_A[nc][5] = - B * (xfr[i + 1] - xfr[i]) * (yfr[j + 1] - yfr[j]);
			}
			else if (i == 0 && j == (nycr - 1)){ // Celda NW (Tipo 6)
				Mat_A[nc][1] = kappa[i] * (yfr[j + 1] - yfr[j]) / (xfr[i + 1] - xcr[i]);
				Mat_A[nc][2] = 0;
				Mat_A[nc][3] = kappa[i] * (yfr[j + 1] - yfr[j]) / (xcr[i] - xfr[i]);
				Mat_A[nc][4] = kappa[i] * (xfr[i + 1] - xfr[i]) / (ycr[j] - ycr[j - 1]);
				Mat_A[nc][0] = - (Mat_A[nc][1] + Mat_A[nc][2] + Mat_A[nc][3] + Mat_A[nc][4]);
				Mat_A[nc][5] = - B * (xfr[i + 1] - xfr[i]) * (yfr[j + 1] - yfr[j]);
			}
			else if (i == 0 && j == 0){ // Celda SW (Tipo 7)
				Mat_A[nc][1] = kappa[i] * (yfr[j + 1] - yfr[j]) / (xcr[i + 1] - xcr[i]);
				Mat_A[nc][2] = kappa[i] * (xfr[i + 1] - xfr[i]) / (ycr[j + 1] - ycr[j]);
				Mat_A[nc][3] = kappa[i] * (yfr[j + 1] - yfr[j]) / (xcr[i] - xfr[i]);
				Mat_A[nc][4] = kappa[i] * (xfr[i + 1] - xfr[i]) / (ycr[j] - yfr[j]);
				Mat_A[nc][0] = - (Mat_A[nc][1] + Mat_A[nc][2] + Mat_A[nc][3] + Mat_A[nc][4]);
				Mat_A[nc][5] = - B * (xfr[i + 1] - xfr[i]) * (yfr[j + 1] - yfr[j]);
			}
			else if (i == (nxcr - 1) && j == 0){ // Celda SE (Tipo 8)
				Mat_A[nc][1] = kappa[i] * (yfr[j + 1] - yfr[j]) / (xfr[i + 1] - xcr[i]);
				Mat_A[nc][2] = kappa[i] * (xfr[i + 1] - xfr[i]) / (ycr[j + 1] - ycr[j]);
				Mat_A[nc][3] = kappa[i] * (yfr[j + 1] - yfr[j]) / (xcr[i] - xcr[i - 1]);
				Mat_A[nc][4] = kappa[i] * (xfr[i + 1] - xfr[i]) / (ycr[j] - yfr[j]);
				Mat_A[nc][0] = - (Mat_A[nc][1] + Mat_A[nc][2] + Mat_A[nc][3] + Mat_A[nc][4]);
				Mat_A[nc][5] = - B * (xfr[i + 1] - xfr[i]) * (yfr[j + 1] - yfr[j]);
			}
			else if ((i > 0 && i < (nxcr - 1)) && (j > 0 && j < (nycr - 1))){ // Celdas interiores (Tipo 9)
				Mat_A[nc][1] = kappa[i] * (yfr[j + 1] - yfr[j]) / (xcr[i + 1] - xcr[i]);
				Mat_A[nc][2] = kappa[i] * (xfr[i + 1] - xfr[i]) / (ycr[j + 1] - ycr[j]);
				Mat_A[nc][3] = kappa[i] * (yfr[j + 1] - yfr[j]) / (xcr[i] - xcr[i - 1]);
				Mat_A[nc][4] = kappa[i] * (xfr[i + 1] - xfr[i]) / (ycr[j] - ycr[j - 1]);
				Mat_A[nc][0] = - (Mat_A[nc][1] + Mat_A[nc][2] + Mat_A[nc][3] + Mat_A[nc][4]);
				Mat_A[nc][5] = - B * (xfr[i + 1] - xfr[i]) * (yfr[j + 1] - yfr[j]);
			}
		}
	}
	return Mat_A;
}

// La función Gauss-Seidel barre la malla de izquierda a derecha yendo de abajo hacia arriba
vector<vector<float>> Gauss_Seidel(float BC_N, vector<float> &BC_S, float BC_W, float BC_E, vector<float> &xcr, vector<float> &ycr, vector<vector<float>> &Mat_A, vector<float> &Vec_Q, vector<float> &T, float epsilon){
	int it = 0;
	float Res_norm = epsilon;
	float Res_norm_ant = 0;
	vector<float> it_V;
	vector<float> Rn_V;
	int nxcr = xcr.size();
	int nycr = ycr.size();
	while (Res_norm >= epsilon){
		float R = 0;
		float F = 0;
		for (int j = 0; j < nycr; j++){
			for (int i = 0; i < nxcr; i++){
				int nc = j * nxcr + i;
				int E = nc + 1;
				int N = nc + nxcr;
				int W = nc - 1;
				int S = nc - nxcr;
				if ((i > 0 && i < (nxcr - 1)) && j == (nycr - 1)){ // Celdas N (Tipo 1)
					T[nc] = (Vec_Q[nc] - (Mat_A[nc][1] * T[E] + Mat_A[nc][3] * T[W] + Mat_A[nc][4] * T[S])) / Mat_A[nc][0];
				}
				else if ((i > 0 && i < (nxcr - 1)) && j == 0){ // Celdas S (Tipo 2)
					T[nc] = (Vec_Q[nc] - (Mat_A[nc][1] * T[E] + Mat_A[nc][2] * T[N] + Mat_A[nc][3] * T[W] + Mat_A[nc][4] * BC_S[i])) / Mat_A[nc][0];
				}
				else if (i == 0 && (j > 0 && j < (nycr - 1))){ // Celdas W (Tipo 3)
					T[nc] = (Vec_Q[nc] - (Mat_A[nc][1] * T[E] + Mat_A[nc][2] * T[N] + Mat_A[nc][3] * BC_W + Mat_A[nc][4] * T[S])) / Mat_A[nc][0];
				}
				else if (i == (nxcr - 1) && (j > 0 && j < (nycr - 1))){ // Celdas E (Tipo 4)
					T[nc] = (Vec_Q[nc] - (Mat_A[nc][1] * BC_E + Mat_A[nc][2] * T[N] + Mat_A[nc][3] * T[W] + Mat_A[nc][4] * T[S])) / Mat_A[nc][0];
				}
				else if (i == (nxcr - 1) && j == (nycr - 1)){ // Celda NE (Tipo 5)
					T[nc] = (Vec_Q[nc] - (Mat_A[nc][1] * BC_E + Mat_A[nc][3] * T[W] + Mat_A[nc][4] * T[S])) / Mat_A[nc][0];
				}
				else if (i == 0 && j == (nycr - 1)){ // Celda NW (Tipo 6)
					T[nc] = (Vec_Q[nc] - (Mat_A[nc][1] * T[E] + Mat_A[nc][3] * BC_W + Mat_A[nc][4] * T[S])) / Mat_A[nc][0];
				}
				else if (i == 0 && j == 0){ // Celda SW (Tipo 7)
					T[nc] = (Vec_Q[nc] - (Mat_A[nc][1] * T[E] + Mat_A[nc][2] * T[N] + Mat_A[nc][3] * BC_W + Mat_A[nc][4] * BC_S[i])) / Mat_A[nc][0];
				}
				else if (i == (nxcr - 1) && j == 0){ // Celda SE (Tipo 8)
					T[nc] = (Vec_Q[nc] - (Mat_A[nc][1] * BC_E + Mat_A[nc][2] * T[N] + Mat_A[nc][3] * T[W] + Mat_A[nc][4] * BC_S[i])) / Mat_A[nc][0];
				}
				else if ((i > 0 && i < (nxcr - 1)) && (j > 0 && j < (nycr - 1))){ // Celdas interiores (Tipo 9)
					T[nc] = (Vec_Q[nc] - (Mat_A[nc][1] * T[E] + Mat_A[nc][2] * T[N] + Mat_A[nc][3] * T[W] + Mat_A[nc][4] * T[S])) / Mat_A[nc][0];
				}
			}
		}
		for (int j = 0; j < nycr; j++){
			for (int i = 0; i < nxcr; i++){
				int nc = j * nxcr + i;
				int E = nc + 1;
				int N = nc + nxcr;
				int W = nc - 1;
				int S = nc - nxcr;
				if ((i > 0 && i < (nxcr - 1)) && j == (nycr - 1)){ // Celdas N (Tipo 1)
					R = R + abs(Vec_Q[nc] - (Mat_A[nc][1] * T[E] + Mat_A[nc][3] * T[W] + Mat_A[nc][4] * T[S] + Mat_A[nc][0] * T[nc]));
				}
				else if ((i > 0 && i < (nxcr - 1)) && j == 0){ // Celdas S (Tipo 2)
					R = R + abs(Vec_Q[nc] - (Mat_A[nc][1] * T[E] + Mat_A[nc][2] * T[N] + Mat_A[nc][3] * T[W] + Mat_A[nc][4] * BC_S[i] + Mat_A[nc][0] * T[nc]));
				}
				else if (i == 0 && (j > 0 && j < (nycr - 1))){ // Celdas W (Tipo 3)
					R = R + abs(Vec_Q[nc] - (Mat_A[nc][1] * T[E] + Mat_A[nc][2] * T[N] + Mat_A[nc][3] * BC_W + Mat_A[nc][4] * T[S] + Mat_A[nc][0] * T[nc]));
				}
				else if (i == (nxcr - 1) && (j > 0 && j < (nycr - 1))){ // Celdas E (Tipo 4)
					R = R + abs(Vec_Q[nc] - (Mat_A[nc][1] * BC_E + Mat_A[nc][2] * T[N] + Mat_A[nc][3] * T[W] + Mat_A[nc][4] * T[S] + Mat_A[nc][0] * T[nc]));
				}
				else if (i == (nxcr - 1) && j == (nycr - 1)){ // Celda NE (Tipo 5)
					R = R + abs(Vec_Q[nc] - (Mat_A[nc][1] * BC_E + Mat_A[nc][3] * T[W] + Mat_A[nc][4] * T[S] + Mat_A[nc][0] * T[nc]));
				}
				else if (i == 0 && j == (nycr - 1)){ // Celda NW (Tipo 6)
					R = R + abs(Vec_Q[nc] - (Mat_A[nc][1] * T[E] + Mat_A[nc][3] * BC_W + Mat_A[nc][4] * T[S] + Mat_A[nc][0] * T[nc]));
				}
				else if (i == 0 && j == 0){ // Celda SW (Tipo 7)
					R = R + abs(Vec_Q[nc] - (Mat_A[nc][1] * BC_E + Mat_A[nc][2] * T[N] + Mat_A[nc][3] * T[W] + Mat_A[nc][4] * BC_S[i] + Mat_A[nc][0] * T[nc]));
				}
				else if (i == (nxcr - 1) && j == 0){ // Celda SE (Tipo 8)
					R = R + abs(Vec_Q[nc] - (Mat_A[nc][1] * T[E] + Mat_A[nc][2] * T[N] + Mat_A[nc][3] * BC_W + Mat_A[nc][4] * BC_S[i] + Mat_A[nc][0] * T[nc]));
				}
				else if ((i > 0 && i < (nxcr - 1)) && (j > 0 && j < (nycr - 1))){ // Celdas interiores (Tipo 9)
					R = R + abs(Vec_Q[nc] - (Mat_A[nc][1] * T[E] + Mat_A[nc][2] * T[N] + Mat_A[nc][3] * T[W] + Mat_A[nc][4] * T[S] + Mat_A[nc][0] * T[nc]));
				}
				F = F + abs(Mat_A[nc][0] * T[nc]);
			}
		}
		Res_norm_ant = Res_norm;
		Res_norm = R / F;
		it++;
		it_V.push_back(0);
		Rn_V.push_back(0);
		it_V[it - 1] = it;
		Rn_V[it - 1] = Res_norm;
		cout << "Iteración " << it << " con Residual Normalizado: " << Res_norm << endl;
		if ((Res_norm_ant - Res_norm) < 1e-8 && (Res_norm_ant - Res_norm) >= 0){
			break;
		}
		else if ((Res_norm_ant - Res_norm) < 0 && it > 1000){
			break;
		}
	}
	int tam_T = T.size();
	int tam_it_V = it_V.size();
	int tam_max;
	if (tam_T > tam_it_V){
		tam_max = tam_T;
	}
	else{
		tam_max = tam_it_V;
	}
	vector<vector<float>> Result(tam_max, vector<float> (3)); // 0: T, 1: it_V, 2: Rn_V
	for (int i = 0; i < T.size(); i++){
		Result[i][0] = T[i];
	}
	for (int i = 0; i < it_V.size(); i++){
		Result[i][1] = it_V[i];
		Result[i][2] = Rn_V[i];
	}
	return Result;
}

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
	fin.open("Celdasy_malla.dat");
		valor = 0;
		vector<float> ycr;
		while (fin >> valor){
			ycr.push_back(valor);
		}
	fin.close();
	fin.open("Carasx_malla.dat");
		valor = 0;
		vector<float> xfr;
		while (fin >> valor){
			xfr.push_back(valor);
		}
	fin.close();
	fin.open("Carasy_malla.dat");
		valor = 0;
		vector<float> yfr;
		while (fin >> valor){
			yfr.push_back(valor);
		}
	fin.close();
	fin.open("Propiedades.dat");
		valor = 0;
		vector<float> kappa;
		while (fin >> valor){
			kappa.push_back(valor);
		}
	fin.close();
	fin.open("Term_src.dat");
		float B;
		fin >> B;
	fin.close();
	fin.open("Fronteras_vrbls.dat");
		valor = 0;
		vector<float> BC_S;
		while (fin >> valor){
			BC_S.push_back(valor);
		}
	fin.close();
	fin.open("Fronteras_ctes.dat");
		valor= 0;
		vector<float> BC;
		while (fin >> valor){
			BC.push_back(valor);
		}
		float BC_E = BC[0];
		float BC_N = BC[1];
		float BC_W = BC[2];
	fin.close();
	// Definición de variables:
	float epsilon = 5e-4; // [-] Criterio de convergencia
	vector<vector<float>> Coef(xcr.size() * ycr.size(), vector<float> (6)); // Matriz de coeficientes y términos independientes
	vector<vector<float>> A(xcr.size() * ycr.size(), vector<float> (5)); // Matriz de coeficientes de las celdas y sus vecinos
	vector<float> Q(xcr.size() * ycr.size()); // Vector de términos independientes para cada celda
	vector<vector<float>> Result; // Matriz de temperatura, iteraciones y residuos normalizados
	vector<float> T(xcr.size() * ycr.size(),10); // [°C] Vector de temperaturas con initial guess = 10
	vector<float> it_V; // Vector de iteraciones
	vector<float> Rn_V; // Vector de residuos normalizados
	// Llamado de las funciones:
	Coef = Coeficientes(xcr, ycr, xfr, yfr, kappa, B);
	for (int i = 0; i < Q.size(); i++){
		A[i][0] = Coef[i][0];
		A[i][1] = Coef[i][1];
		A[i][2] = Coef[i][2];
		A[i][3] = Coef[i][3];
		A[i][4] = Coef[i][4];
		Q[i] = Coef[i][5];
	}
	Result = Gauss_Seidel(BC_N, BC_S, BC_W, BC_E, xcr, ycr, A, Q, T, epsilon);
	for (int i = 0; i < T.size(); i++){
		T[i] = Result[i][0];
	}
	for (int i = 0; i < Result.size(); i++){
		if(i == 0 && (Result[i + 1][1] - Result[i][1]) > 0){
			it_V.push_back(0);
			Rn_V.push_back(0);
			it_V[i] = Result[i][1];
			Rn_V[i] = Result[i][2];
		}
		else if (i > 0 && (Result[i][1] - Result[i - 1][1]) > 0){
			it_V.push_back(0);
			Rn_V.push_back(0);
			it_V[i] = Result[i][1];
			Rn_V[i] = Result[i][2];
		}
	}
	// Escritura de la información:
	ofstream fout;
	fout.open("Temperatura.dat");
		for (int i = 0; i < T.size(); i++){
			fout << T[i] << endl;
		}
	fout.close();
	fout.open("Residuo_norm.dat");
		for (int i = 0; i < it_V.size(); i++){
			fout << it_V[i] << "\t" << Rn_V[i] << endl;
		}
	fout.close();
	return 0;
}
