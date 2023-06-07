#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <chrono>

// Función generadora de la matriz de coeficientes de los vecinos y de la celda y del vector de términos fuente:
Eigen::MatrixXd Coeficientes(Eigen::VectorXd &xcr, Eigen::VectorXd &ycr, Eigen::VectorXd &xfr, Eigen::VectorXd &yfr, Eigen::VectorXd &kappa, double B){
  // Declaración de las variables
  int nxcr = xcr.size();
  int nycr = ycr.size();
  Eigen::MatrixXd Mat_A(nxcr * nycr,6); // 0: A_P, 1: A_E, 2: A_N, 3: A_W, 4: A_S, 5: Q_P
  for (int j = 0; j < nycr; j++){
    for (int i = 0; i < nxcr; i++){
      int nc = j * nxcr + i;
      if ((i > 0 && i < (nxcr - 1)) && j == (nycr - 1)){ // Celdas N (Tipo 1)
	Mat_A(nc,1) = kappa(i) * (yfr(j + 1) - yfr(j)) / (xcr(i + 1) - xcr(i));
	Mat_A(nc,2) = 0;
	Mat_A(nc,3) = kappa(i) * (yfr(j + 1) - yfr(j)) / (xcr(i) - xcr(i - 1));
	Mat_A(nc,4) = kappa(i) * (xfr(i + 1) - xfr(i)) / (ycr(j) - ycr(j - 1));
	Mat_A(nc,0) = - (Mat_A(nc,1) + Mat_A(nc,2) + Mat_A(nc,3) + Mat_A(nc,4));
	Mat_A(nc,5) = - B * (xfr(i + 1) - xfr(i)) * (yfr(j + 1) - yfr(j));
      }
      else if ((i > 0 && i < (nxcr - 1)) && j == 0){ // Celdas S (Tipo 2)
	Mat_A(nc,1) = kappa(i) * (yfr(j + 1) - yfr(j)) / (xcr(i + 1) - xcr(i));
	Mat_A(nc,2) = kappa(i) * (xfr(i + 1) - xfr(i)) / (ycr(j + 1) - ycr(j));
	Mat_A(nc,3) = kappa(i) * (yfr(j + 1) - yfr(j)) / (xcr(i) - xcr(i - 1));
	Mat_A(nc,4) = kappa(i) * (xfr(i + 1) - xfr(i)) / (ycr(j) - yfr(j));
	Mat_A(nc,0) = - (Mat_A(nc,1) + Mat_A(nc,2) + Mat_A(nc,3) + Mat_A(nc,4));
	Mat_A(nc,5) = - B * (xfr(i + 1) - xfr(i)) * (yfr(j + 1) - yfr(j));
      }
      else if (i == 0 && (j > 0 && j < (nycr - 1))){ // Celdas W (Tipo 3)
	Mat_A(nc,1) = kappa(i) * (yfr(j + 1) - yfr(j)) / (xcr(i + 1) - xcr(i));
	Mat_A(nc,2) = kappa(i) * (xfr(i + 1) - xfr(i)) / (ycr(j + 1) - ycr(j));
	Mat_A(nc,3) = kappa(i) * (yfr(j + 1) - yfr(j)) / (xcr(i) - xfr(i));
	Mat_A(nc,4) = kappa(i) * (xfr(i + 1) - xfr(i)) / (ycr(j) - ycr(j - 1));
	Mat_A(nc,0) = - (Mat_A(nc,1) + Mat_A(nc,2) + Mat_A(nc,3) + Mat_A(nc,4));
	Mat_A(nc,5) = - B * (xfr(i + 1) - xfr(i)) * (yfr(j + 1) - yfr(j));
      }
      else if (i == (nxcr - 1) && (j > 0 && j < (nycr - 1))){ // Celdas E (Tipo 4)
	Mat_A(nc,1) = kappa(i) * (yfr(j + 1) - yfr(j)) / (xfr(i + 1) - xcr(i));
	Mat_A(nc,2) = kappa(i) * (xfr(i + 1) - xfr(i)) / (ycr(j + 1) - ycr(j));
	Mat_A(nc,3) = kappa(i) * (yfr(j + 1) - yfr(j)) / (xcr(i) - xcr(i - 1));
	Mat_A(nc,4) = kappa(i) * (xfr(i + 1) - xfr(i)) / (ycr(j) - ycr(j - 1));
	Mat_A(nc,0) = - (Mat_A(nc,1) + Mat_A(nc,2) + Mat_A(nc,3) + Mat_A(nc,4));
	Mat_A(nc,5) = - B * (xfr(i + 1) - xfr(i)) * (yfr(j + 1) - yfr(j));
      }
      else if (i == (nxcr - 1) && j == (nycr - 1)){ // Celda NE (Tipo 5)
	Mat_A(nc,1) = kappa(i) * (yfr(j + 1) - yfr(j)) / (xfr(i + 1) - xcr(i));
	Mat_A(nc,2) = 0;
	Mat_A(nc,3) = kappa(i) * (yfr(j + 1) - yfr(j)) / (xcr(i) - xcr(i - 1));
	Mat_A(nc,4) = kappa(i) * (xfr(i + 1) - xfr(i)) / (ycr(j) - ycr(j - 1));
	Mat_A(nc,0) = - (Mat_A(nc,1) + Mat_A(nc,2) + Mat_A(nc,3) + Mat_A(nc,4));
	Mat_A(nc,5) = - B * (xfr(i + 1) - xfr(i)) * (yfr(j + 1) - yfr(j));
      }
      else if (i == 0 && j == (nycr - 1)){ // Celda NW (Tipo 6)
	Mat_A(nc,1) = kappa(i) * (yfr(j + 1) - yfr(j)) / (xfr(i + 1) - xcr(i));
	Mat_A(nc,2) = 0;
	Mat_A(nc,3) = kappa(i) * (yfr(j + 1) - yfr(j)) / (xcr(i) - xfr(i));
	Mat_A(nc,4) = kappa(i) * (xfr(i + 1) - xfr(i)) / (ycr(j) - ycr(j - 1));
	Mat_A(nc,0) = - (Mat_A(nc,1) + Mat_A(nc,2) + Mat_A(nc,3) + Mat_A(nc,4));
	Mat_A(nc,5) = - B * (xfr(i + 1) - xfr(i)) * (yfr(j + 1) - yfr(j));
      }
      else if (i == 0 && j == 0){ // Celda SW (Tipo 7)
	Mat_A(nc,1) = kappa(i) * (yfr(j + 1) - yfr(j)) / (xcr(i + 1) - xcr(i));
	Mat_A(nc,2) = kappa(i) * (xfr(i + 1) - xfr(i)) / (ycr(j + 1) - ycr(j));
	Mat_A(nc,3) = kappa(i) * (yfr(j + 1) - yfr(j)) / (xcr(i) - xfr(i));
	Mat_A(nc,4) = kappa(i) * (xfr(i + 1) - xfr(i)) / (ycr(j) - yfr(j));
	Mat_A(nc,0) = - (Mat_A(nc,1) + Mat_A(nc,2) + Mat_A(nc,3) + Mat_A(nc,4));
	Mat_A(nc,5) = - B * (xfr(i + 1) - xfr(i)) * (yfr(j + 1) - yfr(j));
      }
      else if (i == (nxcr - 1) && j == 0){ // Celda SE (Tipo 8)
	Mat_A(nc,1) = kappa(i) * (yfr(j + 1) - yfr(j)) / (xfr(i + 1) - xcr(i));
	Mat_A(nc,2) = kappa(i) * (xfr(i + 1) - xfr(i)) / (ycr(j + 1) - ycr(j));
	Mat_A(nc,3) = kappa(i) * (yfr(j + 1) - yfr(j)) / (xcr(i) - xcr(i - 1));
	Mat_A(nc,4) = kappa(i) * (xfr(i + 1) - xfr(i)) / (ycr(j) - yfr(j));
	Mat_A(nc,0) = - (Mat_A(nc,1) + Mat_A(nc,2) + Mat_A(nc,3) + Mat_A(nc,4));
	Mat_A(nc,5) = - B * (xfr(i + 1) - xfr(i)) * (yfr(j + 1) - yfr(j));
      }
      else if ((i > 0 && i < (nxcr - 1)) && (j > 0 && j < (nycr - 1))){ // Celdas interiores (Tipo 9)
	Mat_A(nc,1) = kappa(i) * (yfr(j + 1) - yfr(j)) / (xcr(i + 1) - xcr(i));
	Mat_A(nc,2) = kappa(i) * (xfr(i + 1) - xfr(i)) / (ycr(j + 1) - ycr(j));
	Mat_A(nc,3) = kappa(i) * (yfr(j + 1) - yfr(j)) / (xcr(i) - xcr(i - 1));
	Mat_A(nc,4) = kappa(i) * (xfr(i + 1) - xfr(i)) / (ycr(j) - ycr(j - 1));
	Mat_A(nc,0) = - (Mat_A(nc,1) + Mat_A(nc,2) + Mat_A(nc,3) + Mat_A(nc,4));
	Mat_A(nc,5) = - B * (xfr(i + 1) - xfr(i)) * (yfr(j + 1) - yfr(j));
      }
    }
  }
  return Mat_A;
}

// La función Gauss-Seidel barre la malla de izquierda a derecha yendo de abajo hacia arriba
Eigen::MatrixXd Gauss_Seidel(double BC_N, Eigen::VectorXd &BC_S, double BC_W, double BC_E, Eigen::VectorXd &xcr, Eigen::VectorXd &ycr, Eigen::MatrixXd &Mat_A, Eigen::VectorXd &Vec_Q, Eigen::VectorXd &T, double epsilon){
  // Declaración de las variables
  int it = 0;
  double Res_norm = 10 * epsilon;
  double Res_norm_ant = 0;
  Eigen::VectorXd it_V;
  Eigen::VectorXd Rn_V;
  int nxcr = xcr.size();
  int nycr = ycr.size();
  while (Res_norm >= epsilon){
    double R = 0;
    double F = 0;
    for (int j = 0; j < nycr; j++){
      for (int i = 0; i < nxcr; i++){
	int nc = j * nxcr + i;
	int E = nc + 1;
	int N = nc + nxcr;
	int W = nc - 1;
	int S = nc - nxcr;
	if ((i > 0 && i < (nxcr - 1)) && j == (nycr - 1)){ // Celdas N (Tipo 1)
	  T(nc) = (Vec_Q(nc) - (Mat_A(nc,1) * T(E) + Mat_A(nc,3) * T(W) + Mat_A(nc,4) * T(S))) / Mat_A(nc,0);
	}
	else if ((i > 0 && i < (nxcr - 1)) && j == 0){ // Celdas S (Tipo 2)
	  T(nc) = (Vec_Q(nc) - (Mat_A(nc,1) * T(E) + Mat_A(nc,2) * T(N) + Mat_A(nc,3) * T(W) + Mat_A(nc,4) * BC_S(i))) / Mat_A(nc,0);
	}
	else if (i == 0 && (j > 0 && j < (nycr - 1))){ // Celdas W (Tipo 3)
	  T(nc) = (Vec_Q(nc) - (Mat_A(nc,1) * T(E) + Mat_A(nc,2) * T(N) + Mat_A(nc,3) * BC_W + Mat_A(nc,4) * T(S))) / Mat_A(nc,0);
	}
	else if (i == (nxcr - 1) && (j > 0 && j < (nycr - 1))){ // Celdas E (Tipo 4)
	  T(nc) = (Vec_Q(nc) - (Mat_A(nc,1) * BC_E + Mat_A(nc,2) * T(N) + Mat_A(nc,3) * T(W) + Mat_A(nc,4) * T(S))) / Mat_A(nc,0);
	}
	else if (i == (nxcr - 1) && j == (nycr - 1)){ // Celda NE (Tipo 5)
	  T(nc) = (Vec_Q(nc) - (Mat_A(nc,1) * BC_E + Mat_A(nc,3) * T(W) + Mat_A(nc,4) * T(S))) / Mat_A(nc,0);
	}
	else if (i == 0 && j == (nycr - 1)){ // Celda NW (Tipo 6)
	  T(nc) = (Vec_Q(nc) - (Mat_A(nc,1) * T(E) + Mat_A(nc,3) * BC_W + Mat_A(nc,4) * T(S))) / Mat_A(nc,0);
	}
	else if (i == 0 && j == 0){ // Celda SW (Tipo 7)
	  T(nc) = (Vec_Q(nc) - (Mat_A(nc,1) * T(E) + Mat_A(nc,2) * T(N) + Mat_A(nc,3) * BC_W + Mat_A(nc,4) * BC_S(i))) / Mat_A(nc,0);
	}
	else if (i == (nxcr - 1) && j == 0){ // Celda SE (Tipo 8)
	  T(nc) = (Vec_Q(nc) - (Mat_A(nc,1) * BC_E + Mat_A(nc,2) * T(N) + Mat_A(nc,3) * T(W) + Mat_A(nc,4) * BC_S(i))) / Mat_A(nc,0);
	}
	else if ((i > 0 && i < (nxcr - 1)) && (j > 0 && j < (nycr - 1))){ // Celdas interiores (Tipo 9)
	  T(nc) = (Vec_Q(nc) - (Mat_A(nc,1) * T(E) + Mat_A(nc,2) * T(N) + Mat_A(nc,3) * T(W) + Mat_A(nc,4) * T(S))) / Mat_A(nc,0);
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
	  R = R + std::abs(Vec_Q(nc) - (Mat_A(nc,1) * T(E) + Mat_A(nc,3) * T(W) + Mat_A(nc,4) * T(S) + Mat_A(nc,0) * T(nc)));
	}
	else if ((i > 0 && i < (nxcr - 1)) && j == 0){ // Celdas S (Tipo 2)
	  R = R + std::abs(Vec_Q(nc) - (Mat_A(nc,1) * T(E) + Mat_A(nc,2) * T(N) + Mat_A(nc,3) * T(W) + Mat_A(nc,4) * BC_S(i) + Mat_A(nc,0) * T(nc)));
	}
	else if (i == 0 && (j > 0 && j < (nycr - 1))){ // Celdas W (Tipo 3)
	  R = R + std::abs(Vec_Q(nc) - (Mat_A(nc,1) * T(E) + Mat_A(nc,2) * T(N) + Mat_A(nc,3) * BC_W + Mat_A(nc,4) * T(S) + Mat_A(nc,0) * T(nc)));
	}
	else if (i == (nxcr - 1) && (j > 0 && j < (nycr - 1))){ // Celdas E (Tipo 4)
	  R = R + std::abs(Vec_Q(nc) - (Mat_A(nc,1) * BC_E + Mat_A(nc,2) * T(N) + Mat_A(nc,3) * T(W) + Mat_A(nc,4) * T(S) + Mat_A(nc,0) * T(nc)));
	}
	else if (i == (nxcr - 1) && j == (nycr - 1)){ // Celda NE (Tipo 5)
	  R = R + std::abs(Vec_Q(nc) - (Mat_A(nc,1) * BC_E + Mat_A(nc,3) * T(W) + Mat_A(nc,4) * T(S) + Mat_A(nc,0) * T(nc)));
	}
	else if (i == 0 && j == (nycr - 1)){ // Celda NW (Tipo 6)
	  R = R + std::abs(Vec_Q(nc) - (Mat_A(nc,1) * T(E) + Mat_A(nc,3) * BC_W + Mat_A(nc,4) * T(S) + Mat_A(nc,0) * T(nc)));
	}
	else if (i == 0 && j == 0){ // Celda SW (Tipo 7)
	  R = R + std::abs(Vec_Q(nc) - (Mat_A(nc,1) * T(E) + Mat_A(nc,2) * T(N) + Mat_A(nc,3) * BC_W + Mat_A(nc,4) * BC_S(i) + Mat_A(nc,0) * T(nc)));
	}
	else if (i == (nxcr - 1) && j == 0){ // Celda SE (Tipo 8)
	  R = R + std::abs(Vec_Q(nc) - (Mat_A(nc,1) * BC_E + Mat_A(nc,2) * T(N) + Mat_A(nc,3) * T(W) + Mat_A(nc,4) * BC_S(i) + Mat_A(nc,0) * T(nc)));
	}
	else if ((i > 0 && i < (nxcr - 1)) && (j > 0 && j < (nycr - 1))){ // Celdas interiores (Tipo 9)
	  R = R + std::abs(Vec_Q(nc) - (Mat_A(nc,1) * T(E) + Mat_A(nc,2) * T(N) + Mat_A(nc,3) * T(W) + Mat_A(nc,4) * T(S) + Mat_A(nc,0) * T(nc)));
	}
	F = F + std::abs(Mat_A(nc,0) * T(nc));
      }
    }
    Res_norm_ant = Res_norm;
    Res_norm = R / F;
    it++;
    it_V.conservativeResize(it);
    Rn_V.conservativeResize(it);
    it_V(it - 1) = it;
    Rn_V(it - 1) = Res_norm;
    std::cout << "Iteración " << it << " con Residual Normalizado: " << Res_norm << "\n";
    if ((Res_norm_ant - Res_norm) < 1e-9 && (Res_norm_ant - Res_norm) >= 0){
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
  Eigen::MatrixXd Result(tam_max,3); // 0: T, 1: it_V, 2: Rn_V
  for (int i = 0; i < T.size(); i++){
    Result(i,0) = T(i);
  }
  for (int i = 0; i < it_V.size(); i++){
    Result(i,1) = it_V(i);
    Result(i,2) = Rn_V(i);
  }
  return Result;
}

int main(int argc, char ** argv) {
  int num_thr = std::stoi(argv[1]); // [-] Número de threads que se utilizan
  // Recopilación de información:
  std::ifstream fin;
  fin.open("Celdasx_malla.dat");
  double valor;
  Eigen::VectorXd xcr;
  while (fin >> valor){
    xcr.conservativeResize(xcr.size() + 1);
    xcr(xcr.size() - 1) = valor;
  }
  fin.close();
  fin.open("Celdasy_malla.dat");
  valor = 0;
  Eigen::VectorXd ycr;
  while (fin >> valor){
    ycr.conservativeResize(ycr.size() + 1);
    ycr(ycr.size() - 1) = valor;
  }
  fin.close();
  fin.open("Carasx_malla.dat");
  valor = 0;
  Eigen::VectorXd xfr;
  while (fin >> valor){
    xfr.conservativeResize(xfr.size() + 1);
    xfr(xfr.size() - 1) = valor;
  }
  fin.close();
  fin.open("Carasy_malla.dat");
  valor = 0;
  Eigen::VectorXd yfr;
  while (fin >> valor){
    yfr.conservativeResize(yfr.size() + 1);
    yfr(yfr.size() - 1) = valor;
  }
  fin.close();
  fin.open("Propiedades.dat");
  valor = 0;
  Eigen::VectorXd kappa;
  while (fin >> valor){
    kappa.conservativeResize(kappa.size() + 1);
    kappa(kappa.size() - 1) = valor;
  }
  fin.close();
  fin.open("Term_src.dat");
  double B;
  fin >> B;
  fin.close();
  fin.open("Fronteras_vrbls.dat");
  valor = 0;
  Eigen::VectorXd BC_S;
  while (fin >> valor){
    BC_S.conservativeResize(BC_S.size() + 1);
    BC_S(BC_S.size() - 1) = valor;
  }
  fin.close();
  fin.open("Fronteras_ctes.dat");
  valor= 0;
  Eigen::VectorXd BC;
  while (fin >> valor){
    BC.conservativeResize(BC.size() + 1);
    BC(BC.size() - 1);
  }
  double BC_E = BC(0);
  double BC_N = BC(1);
  double BC_W = BC(2);
  fin.close();
  // Inicialización de los tiempos
  std::clock_t cpu_prcss_start = std::clock();
  auto wall_prcss_start = std::chrono::steady_clock::now();
  // Definición de variables:
  double epsilon = 1e-7; // [-] Criterio de convergencia
  Eigen::MatrixXd Coef(xcr.size() * ycr.size(),6); // Matriz de coeficientes y términos independientes
  Eigen::MatrixXd A(xcr.size() * ycr.size(),5); // Matriz de coeficientes de las celdas y sus vecinos
  Eigen::VectorXd Q(xcr.size() * ycr.size()); // Vector de términos independientes para cada celda
  Eigen::MatrixXd Result; // Matriz de temperatura, iteraciones y residuos normalizados
  Eigen::VectorXd T;
  T.setConstant(xcr.size() * ycr.size(),10); // [°C] Vector de temperaturas con initial guess = 10
  Eigen::VectorXd it_V; // Vector de iteraciones
  Eigen::VectorXd Rn_V; // Vector de residuos normalizados
  // Llamado de las funciones:
  Coef = Coeficientes(xcr, ycr, xfr, yfr, kappa, B);
  for (int i = 0; i < Q.size(); i++){
    A(i,0) = Coef(i,0);
    A(i,1) = Coef(i,1);
    A(i,2) = Coef(i,2);
    A(i,3) = Coef(i,3);
    A(i,4) = Coef(i,4);
    Q(i) = Coef(i,5);
  }
  Result = Gauss_Seidel(BC_N, BC_S, BC_W, BC_E, xcr, ycr, A, Q, T, epsilon);
  for (int i = 0; i < T.size(); i++){
    T(i) = Result(i,0);
  }
  for (int i = 0; i < Result.rows(); i++){
    if(i == 0 && (Result(i + 1,1) - Result(i,1)) > 0){
      it_V.conservativeResize(it_V.size() + 1);
      Rn_V.conservativeResize(Rn_V.size() + 1);
      it_V(i) = Result(i,1);
      Rn_V(i) = Result(i,2);
    }
    else if (i > 0 && (Result(i,1) - Result(i - 1,1)) > 0){
      it_V.conservativeResize(it_V.size() + 1);
      Rn_V.conservativeResize(Rn_V.size() + 1);
      it_V(i) = Result(i,1);
      Rn_V(i) = Result(i,2);
    }
  }
  // Finalización de los tiempos
  auto wall_prcss_end = std::chrono::steady_clock::now();
  std::chrono::duration<double> wall_prcss_elapsed_sec = wall_prcss_end - wall_prcss_start;
  std::clock_t cpu_prcss_end = std::clock();
  // Escritura de la información:
  std::ofstream fout;
  fout.open("Temperatura.dat");
  for (int i = 0; i < T.size(); i++){
    fout << T(i) << "\n";
  }
  fout.close();
  fout.open("Residuo_norm.dat");
  for (int i = 0; i < it_V.size(); i++){
    fout << it_V(i) << "\t" << Rn_V(i) << "\n";
  }
  fout.close();
  fout.open("Parallel_time_prcss.dat", std::ios::app);
  fout << num_thr << "\t" << wall_prcss_elapsed_sec.count() << "\t" << double (cpu_prcss_end - cpu_prcss_start) / CLOCKS_PER_SEC << "\n";
  fout.close();
  return 0;
}
