#!/usr/bin/bash
cd Caso_analitico/
make -f Caso_analitico.mk
cd ../Caso_interes/
make -f Caso_interes.mk
cd ../
