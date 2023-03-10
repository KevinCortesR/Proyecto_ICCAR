#include <iostream>

int main() {
  std::cout << "¡Hola, mundo!" << std::endl;
  std::cout << "¿Qué tal si jugamos un juego?" << std::endl;
  std::cout << "Vamos a contar hasta 10: ";
  
  for(int i = 1; i <= 10; i++) {
    std::cout << i << ", ";
  }
  
  std::cout << "¡Listo o no, ahí voy!" << std::endl;
  
  return 0;
}
