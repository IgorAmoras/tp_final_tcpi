s = tf('s');
%% Dados do sistema
% Sistema linearizado
h11 = -0.1806*exp(-89.5*s)/(150*s + 1);
h21 = -0.8357*exp(-220*s)/(580*s + 1);
h12 = -0.05705*exp(-101*s)/(140*s + 1);
h22 = 0.134*exp(-180*s)/(610*s + 1);

H = [h11 h12; h21 h22];

% Sistema sem atraso
h11_nd = -0.1806/(150*s + 1);
h21_nd = -0.8357/(580*s + 1);
h12_nd = -0.05705/(140*s + 1);
h22_nd = 0.134/(610*s + 1);

H_nd = [h11_nd h12_nd; h21_nd h22_nd];

%% Controladores
% Controlador por desacoplamento
[num_tf, dent_f] = tfdata(H_nd);
num = cell2mat(num_tf);
inv_num = inv(num(:,[2,4])); 

decoupler = cell(2,2);

for i = 1:2
    for j = 1:2
        decoupler{i, j} = tf(inv_num(i, j), dent_f{i, j});
    end
end

%% Teste de desacoplamento

decoupled_system = MultiplyMIMO(decoupler, ToCellArray(H_nd));
decoupled_system = ToTfMatrix(decoupled_system);
step(decoupled_system)
function result = MultiplyMIMO(A, B)
    [m, n] = size(A);
    [~, q] = size(B);
    result = cell(m, q);
    for i = 1:m
        for j = 1:q
            result_tf = tf(0);
            for k = 1:n
                result_tf = result_tf + A{i, k} * B{k, j};
            end
            result{i, j} = result_tf;
        end
    end
end

function cell_array = ToCellArray(M)
    [m, n] = size(M);
    cell_array = cell(m, n);

    for i = 1:m
        for j = 1:n
            cell_array{i, j} = M(i, j);
        end
    end
end

function tf_matrix = ToTfMatrix(cell_array)
    [m, n] = size(cell_array);
    for i = 1:m
        for j = 1:n
            tf_matrix(i, j) = cell_array{i, j};
        end
    end
end
