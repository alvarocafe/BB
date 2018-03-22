% Mostrar numeração dos nós
  for no  = 1 : length(NOS(:,1))
    n_no = num2str(NOS(no,1));
    x_no = NOS(no,2);
    y_no = NOS(no,3);
    text(x_no,y_no,n_no);
  end;
