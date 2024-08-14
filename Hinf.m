% Gera os dados do sistema
DataSysLin

% Sistema 1 
H1 = decoupled_system(1,:);
[num_d, den_d] = tfdata(H1);

[A,B,C,D] = tf2ss(num_d{1}, den_d{1});

res = Lema31Finsler(A,B,eye(2,1),C,D,zeros(1,1));
