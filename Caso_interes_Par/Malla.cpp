#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <algorithm>
#include <chrono>
#include <numeric>

// Para realizar un enmallado uniforme haga Delta = Delta_mean.
Eigen::VectorXd refinamiento(double lim_inf, double lim_sup, Eigen::VectorXd &p_int, double Delta_max, double Delta){
  if (Delta_max == Delta){
    p_int.resize(1);
    p_int(0) = 0.5 * (lim_sup + lim_inf);
  }
  int n_min_c = round((lim_sup - lim_inf) / Delta_max); // [-] Número mínimo de celdas
  int n_max_c = round((lim_sup - lim_inf) / Delta); // [-] Número máximo de celdas
  double r_max = n_min_c * Delta; // [m] Distáncia máxima de influencia de los puntos de interés
  double r_min = (r_max / 10) * (double(n_max_c - n_min_c) / double(n_max_c)); // [m] Distancia alrededor de los puntos de interés donde debe darse el tamaño mínimo de las celdas
  Eigen::VectorXd zoi_d(6 * p_int.size()); // [m] Vector de zonas de influencia en términos de distancia nominal
  zoi_d.setZero();
  Eigen::VectorXi zoi_c(6 * p_int.size()); // [-] Vector de zonas de influencia en términos de celdas
  zoi_c.setZero();
  Eigen::VectorXd zoi_dc(6 * p_int.size()); // [m] Vector de zonas de influencia en términos de distancia de las celdas
  zoi_dc.setZero();
  Eigen::VectorXd scl_fac(6 * p_int.size()); // [-] Vecotr de factores de escala de las zonas de influencia.
  scl_fac.setOnes();
  double zoi_sr; // Zona de influencia sin refinamiento
  double zoi_rf; // Zona de influencia con refinamiento
  double zoi_tr; // Zona de influencia en transición
  int num_csr; // Número de celdas en zoi_sr
  int num_crf; // Número de celdas en zoi_rf
  int num_ctr; // Número de celdas en zoi_tr
  
  for (int i = 0; i < p_int.size(); i++){ // Cálculo del número de celdas y las distancias nominales de las zonas de influencia
    if (i == 0){
      double r_left = p_int(i) - lim_inf;
      if(r_left > r_max){
	zoi_sr = r_left - r_max;
	zoi_rf = r_min;
	num_csr = ceil(zoi_sr / Delta_max);
	zoi_sr = num_csr * Delta_max;
	num_crf = floor(zoi_rf / Delta);
	zoi_rf = num_crf * Delta;
	zoi_tr = r_left - zoi_sr - zoi_rf;
	num_ctr = ceil((10 * zoi_tr) / (7 * Delta_max + 3 * Delta));
      }
      else if(r_left <= r_max && r_left > r_min){
	zoi_sr = 0;
	zoi_rf = r_min;
	num_csr = ceil(zoi_sr / Delta_max);
	zoi_sr = num_csr * Delta_max;
	num_crf = floor(zoi_rf / Delta);
	zoi_rf = num_crf * Delta;
	zoi_tr = r_left - zoi_sr - zoi_rf;
	num_ctr = ceil((10 * zoi_tr) / (7 * Delta_max + 3 * Delta));
      }
      else if(r_left <= r_min){
	zoi_sr = 0;
	zoi_rf = r_left;
	num_csr = ceil(zoi_sr / Delta_max);
	zoi_sr = num_csr * Delta_max;
	num_crf = floor(zoi_rf / Delta);
	zoi_rf = num_crf * Delta;
	zoi_tr = r_left - zoi_sr - zoi_rf;
	num_ctr = ceil((10 * zoi_tr) / (7 * Delta_max + 3 * Delta));
      }
      zoi_d(6 * i) = zoi_sr;
      zoi_d(6 * i + 1) = zoi_tr;
      zoi_d(6 * i + 2) = zoi_rf;
      zoi_c(6 * i) = num_csr;
      zoi_c(6 * i + 1) = num_ctr;
      zoi_c(6 * i + 2) = num_crf;
    }
    if (i == (p_int.size() - 1)){
      double r_right = lim_sup - p_int(i);
      if(r_right > r_max){
	zoi_sr = r_right - r_max;
	zoi_rf = r_min;
	num_csr = ceil(zoi_sr / Delta_max);
	zoi_sr = num_csr * Delta_max;
	num_crf = floor(zoi_rf / Delta);
	zoi_rf = num_crf * Delta;
	zoi_tr = r_right - zoi_sr - zoi_rf;
	num_ctr = ceil((10 * zoi_tr) / (7 * Delta_max + 3 * Delta));
      }
      else if(r_right <= r_max && r_right > r_min){
	zoi_sr = 0;
	zoi_rf = r_min;
	num_csr = ceil(zoi_sr / Delta_max);
	zoi_sr = num_csr * Delta_max;
	num_crf = floor(zoi_rf / Delta);
	zoi_rf = num_crf * Delta;
	zoi_tr = r_right - zoi_sr - zoi_rf;
	num_ctr = ceil((10 * zoi_tr) / (7 * Delta_max + 3 * Delta));
      }
      else if(r_right <= r_min){
	zoi_sr = 0;
	zoi_rf = r_right;
	num_csr = ceil(zoi_sr / Delta_max);
	zoi_sr = num_csr * Delta_max;
	num_crf = floor(zoi_rf / Delta);
	zoi_rf = num_crf * Delta;
	zoi_tr = r_right - zoi_sr - zoi_rf;
	num_ctr = ceil((10 * zoi_tr) / (7 * Delta_max + 3 * Delta));
      }
      zoi_d(6 * i + 3) = zoi_rf;
      zoi_d(6 * i + 4) = zoi_tr;
      zoi_d(6 * i + 5) = zoi_sr;
      zoi_c(6 * i + 3) = num_crf;
      zoi_c(6 * i + 4) = num_ctr;
      zoi_c(6 * i + 5) = num_csr;
    }
    if (p_int.size() > 1 && i < (p_int.size() - 1)){
      double r_mid = (p_int(i+1) - p_int(i)) / 2;
      if(r_mid > r_max){
	zoi_sr = r_mid - r_max;
	zoi_rf = r_min;
	num_csr = ceil(zoi_sr / Delta_max);
	zoi_sr = num_csr * Delta_max;
	num_crf = floor(zoi_rf / Delta);
	zoi_rf = num_crf * Delta;
	zoi_tr = r_mid - zoi_sr - zoi_rf;
	num_ctr = ceil((10 * zoi_tr) / (3 * Delta_max + 7 * Delta));
      }
      else if(r_mid <= r_max && r_mid > r_min){
	zoi_sr = 0;
	zoi_rf = r_min;
	num_csr = ceil(zoi_sr / Delta_max);
	zoi_sr = num_csr * Delta_max;
	num_crf = floor(zoi_rf / Delta);
	zoi_rf = num_crf * Delta;
	zoi_tr = r_mid - zoi_sr - zoi_rf;
	num_ctr = ceil((10 * zoi_tr) / (3 * Delta_max + 7 * Delta));
      }
      else if(r_mid <= r_min){
	zoi_sr = 0;
	zoi_rf = r_mid;
	num_csr = ceil(zoi_sr / Delta_max);
	zoi_sr = num_csr * Delta_max;
	num_crf = floor(zoi_rf / Delta);
	zoi_rf = num_crf * Delta;
	zoi_tr = r_mid - zoi_sr - zoi_rf;
	num_ctr = ceil((10 * zoi_tr) / (3 * Delta_max + 7 * Delta));
      }
      zoi_d(6 * i + 3) = zoi_rf;
      zoi_d(6 * i + 4) = zoi_tr;
      zoi_d(6 * i + 5) = zoi_sr;
      zoi_d(6 * (i + 1)) = zoi_sr;
      zoi_d(6 * (i + 1) + 1) = zoi_tr;
      zoi_d(6 * (i + 1) + 2) = zoi_rf;
      zoi_c(6 * i + 3) = num_crf;
      zoi_c(6 * i + 4) = num_ctr;
      zoi_c(6 * i + 5) = num_csr;
      zoi_c(6 * (i + 1)) = num_csr;
      zoi_c(6 * (i + 1) + 1) = num_ctr;
      zoi_c(6 * (i + 1) + 2) = num_crf;
    }
  }

  if (Delta_max == Delta){ // Compensación de celdas cuando el tamaño de la celda es uniforme
    int dif = n_min_c - std::reduce(zoi_c.begin(), zoi_c.end(), 0);
    int dif_left = floor(dif / 2);
    int dif_right = ceil(dif / 2);
    zoi_c(1) += dif_left;
    zoi_c(4) += dif_right;
  }

  int num_cell = std::reduce(zoi_c.begin(), zoi_c.end(), 0);
  Eigen::VectorXd x_ref_aux(num_cell + 1);
  x_ref_aux(0) = lim_inf;
  Eigen::VectorXd x_ref(num_cell + 1);
  x_ref(0) = lim_inf;
  int cont = 0;
  for (int i = 0; i < zoi_c.size(); i++){ // Llenado del vector de caras refinado sin corrección
    if (zoi_c(i) > 0){
      for (int j = 0; j < zoi_c(i); j++){
	++cont;
	if ((i % 6) == 0 || ((i - 5) % 6) == 0){ // Sin refinamiento (sr)
	  x_ref_aux(cont) = x_ref_aux(cont - 1) + Delta_max;
	}
	else if (((i - 1) % 6) == 0){ // Transición (tr) en acercamiento
	  double Delta_celda = ((Delta_max - Delta) / (r_max - r_min)) * ((7 * Delta_max + 3 * Delta) / 10) * (zoi_c(i) - j) + Delta;
	  x_ref_aux(cont) = x_ref_aux(cont - 1) + Delta_celda;
	}
	else if (((i - 4) % 6) == 0){ // Transición (tr) en alejamiento
	  double Delta_celda = ((Delta_max - Delta) / (r_max - r_min)) * ((7 * Delta_max + 3 * Delta) / 10) * (j + 1) + Delta;
	  x_ref_aux(cont) = x_ref_aux(cont - 1) + Delta_celda;
	}
	else if (((i - 2) % 6) == 0 || ((i - 3) % 6) == 0){ // Refinamiento (rf)
	  x_ref_aux(cont) = x_ref_aux(cont - 1) + Delta;
	}
      }
    }
  }

  int c_in = 0; // Celda inicial para la escritura de los vectores x_ref_aux y x_ref
  double dist;

  for (int i = 0; i < zoi_c.size(); i++){ // Cálculo de las distancias con las celdas obtenidas
    if (zoi_c(i) > 0){
      dist = 0;
      for (int j = 0; j < zoi_c(i); j++){
	dist += x_ref_aux(c_in + j + 1) - x_ref_aux(c_in + j);
      }
      zoi_dc(i) = dist;
      c_in += zoi_c(i);
    }
  }

  if (Delta_max != Delta){
    for (int i = 0; i < scl_fac.size(); i++){ // Cálculo de los factores de escala de cada zona de influencia
      if (zoi_dc(i) > 0){
	scl_fac(i) = zoi_d(i) / zoi_dc(i);
      }
      else {
	scl_fac(i) = 0;
      }
    }
  }

  cont = 0;
  c_in = 0;
  double celda;

  for (int i = 0; i < zoi_c.size(); i++){ // Aplicación de los factores de escala en cada zona de influencia
    if (zoi_c(i) > 0){
      for (int j = 0; j < zoi_c(i); j++){
	++cont;
	celda = x_ref_aux(c_in + j + 1) - x_ref_aux(c_in + j);
	x_ref(cont) = x_ref(cont - 1) + scl_fac(i) * celda;
      }
      c_in += zoi_c(i);
    }
  }
  
  return x_ref;
}

int main(int argc, char ** argv) {
  // Inicialización de los tiempos
  std::clock_t cpu_mesh_start = std::clock();
  auto wall_mesh_start = std::chrono::steady_clock::now();
  // Parámetros geométricos:
  double L = 3; // [m] Tamaño en x
  double H = 2; // [m] Tamaño en y
  // Parámetros de la malla:
  int ncx = 180; // [-] Número de celdas en x
  int ncy = 120; // [-] Número de celdas en y
  double Deltax_max = L / double(ncx); // [m] Tamaño máximo de las celdas en x
  double Deltay_max = H / double(ncy); // [m] Tamaño máximo de las celdas en y
  double Deltax = 1 * Deltax_max; // [m] Tamaño mínimo de las celdas en x
  double Deltay = 1 * Deltay_max; // [m] Tamaño mínimo de las celdas en y
  int nfx = ncx + 1; // [-] Número de caras en x
  int nfy = ncy + 1; // [-] Número de caras en y
  Eigen::VectorXd px_int {{0, L}}; // [m] Puntos de interés para refinar en dirección x
  Eigen::VectorXd py_int {{0, 0.05 * H, H}}; // [m] Puntos de interés para refinar en dirección y
  // Creación de la malla:
  Eigen::VectorXd xfr = refinamiento(0, L, px_int, Deltax_max, Deltax); // [m] Vector de caras en x con refinamiento
  Eigen::VectorXd yfr = refinamiento(0, H, py_int, Deltay_max, Deltay); // [m] Vector de caras en y con refinamiento
  Eigen::VectorXd xcr(xfr.size() - 1); // [m] Vector de celdas (centroides) en x con refinamiento
  Eigen::VectorXd ycr(yfr.size() - 1); // [m] Vector de celdas (centroides) en y con refinamiento
  for (int i = 0; i < xcr.size(); i++){
    xcr(i) = 0.5 * (xfr(i) + xfr(i + 1));
  }
  for (int i = 0; i < ycr.size(); i++){
    ycr(i) = 0.5 * (yfr(i) + yfr(i + 1));
  }
  // Finalización de los tiempos
  auto wall_mesh_end = std::chrono::steady_clock::now();
  std::chrono::duration<double> wall_mesh_elapsed_sec = wall_mesh_end - wall_mesh_start;
  std::clock_t cpu_mesh_end = std::clock();
  // Escritura de la información:
  std::ofstream file;
  file.open("Carasx_malla.dat");
  for (int i = 0; i < xfr.size(); i++){
    file << xfr(i) << "\n";
  }
  file.close();
  file.open("Carasy_malla.dat");
  for (int i = 0; i < yfr.size(); i++){
    file << yfr(i) << "\n";
  }
  file.close();
  file.open("Celdasx_malla.dat");
  for (int i = 0; i < xcr.size(); i++){
    file << xcr(i) << "\n";
  }
  file.close();
  file.open("Celdasy_malla.dat");
  for (int i = 0; i < ycr.size(); i++){
    file << ycr(i) << "\n";
  }
  file.close();
  file.open("Time_mesh.dat", std::ios::app);
  file << wall_mesh_elapsed_sec.count() << "\t" << double (cpu_mesh_end - cpu_mesh_start) / CLOCKS_PER_SEC << "\n";
  file.close();
  return 0;
}
