function mapa_de_cor(NOS,ELEM,T)
% geometria da estrutura
figure
%patch('faces',ELEM,'vertices',NOS, ...
%    'FaceVertexCData',T,'CDataMapping','scaled','FaceColor','flat');
patch('faces',ELEM,'vertices',NOS, ...
    'FaceVertexCData',T,'CDataMapping','scaled', 'facecolor','interp', ...
    'EdgeColor','none');
view(2)
%view(3)
colorbar
axis image
xlabel('x')
ylabel('y')
%zlabel('z')
colorbar

% Apagar linhas das arestas

% For FaceVertexCData to work properly, we have to use an m-by-1 column
% vector for indexed colors or m-by-3 array of RBG triplets, 
% where m is the number of rows in the "faces" property

% 1ª Tentativa - Erro: Não plota figura
%size(NOS) = 256 3
%size(ELEM) = 68 2
%size(T) = 68 1

% 2ª Tentativa - Erro: Plota apenas o contorno da figura
%size(NOS) = 68 3
%size(ELEM) = 68 2
%size(T) = 68 1

% Teste bracket - Funciona OK. 3D
%size(NOS) = 25266 3
%size(ELEM) = 8422 3
%size(T) = 8422 1
