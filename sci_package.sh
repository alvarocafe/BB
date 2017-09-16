#Bash script para instalar programas cientificos para o Linux Mint.
#Esse shell script instala as dependencias e os programas necessarios
#para analise numerica utilizando o metodo dos elementos de contorno.
#Tambem pretende-se instalar e clonar o repositorio base do metodo dos elementos 
#de contorno.
#Autor: Alvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com

echo Instalar dependencias:
echo **********************Atualizar repositorios**********************
sudo apt-get update
echo **********************Instalar pacotes**********************
sudo apt-get install texlive-latex-base
sudo apt-get install texstudio
sudo apt-get install julia
sudo apt-get install git
sudo apt-get install freecad
sudo apt-get install gmsh
#echo **********************Fim da instalacao dos pacotes**********************
echo **********************Atualizar os softwares**********************
sudo apt-get update
sudo apt-get upgrade-y
echo **********************Clonar programas de BEM pelo GitHub.**********************
git clone https://github.com/alvarocafe/BEM_base.git
