#Bash script para instalar programas cientificos para o Linux Mint.
#Esse shell script instala as dependencias e os programas necessarios
#para analise numerica utilizando o metodo dos elementos de contorno.
#Tambem pretende-se instalar e clonar o repositorio base do metodo dos elementos 
#de contorno.
#Autor: Alvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
#Tela de inicio
echo '              
  __   ____  ____  __ _                                           
 /  \ (  _ \(  __)(  ( \                                          
(  O ) ) __/ ) _) /    /                                          
 \__/ (__)  (____)\_)__)                                          
 ____   __   _  _  _  _  __     __   ____   __    __   __ _  ____ 
/ ___) (  ) ( \/ )/ )( \(  )   / _\ (_  _) (  )  /  \ (  ( \/ ___)
\___ \  )(  / \/ \) \/ (/ (_/\/    \  )(    )(  (  O )/    /\___ \
(____/ (__) \_)(_/\____/\____/\_/\_/ (__)  (__)  \__/ \_)__)(____/
                                                                  
'
#echo Instalar dependencias:
#sudo add-apt-repository ppa:webupd8team/atom #retirado para a eficiencia.
echo **********************Atualizar repositorios**********************
sudo apt-get update
sudo apt-get upgrade-y
echo **********************Instalar pacotes**********************
sudo apt-get install texlive-latex-base
#sudo apt-get install atom #retirado para a eficiencia.Gedit já faz o necessario.
sudo apt-get install texstudio
sudo apt-get install julia
sudo apt-get install octave
sudo apt-get install git
sudo apt-get install freecad
sudo apt-get install gmsh
#echo **********************Fim da instalacao dos pacotes**********************
echo **********************Atualizar os softwares**********************
sudo apt-get update
sudo apt-get upgrade-y
echo **********************Clonar programas de BEM pelo GitHub.**********************
git clone https://github.com/alvarocafe/BEM_base.git
