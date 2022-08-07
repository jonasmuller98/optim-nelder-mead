%% Engenharia Aeroespacial - Algoritmos de Otimização
% Nelder Mead - 2-Simplex
% Jonas Müller Gonçalves
clear all; clc;

iusuario = input('Deseja inserir valores manuais (=1) ou teste cadastrado(~=1)? ');
if iusuario == 1
    iusuario2 = 0;
    fteste = inputdlg('Digite a função (a(x(1)) + b(x(2)) ...): '); s = fteste{:};f = str2func(['@(x)' s]);
    alfa = input('Insira o valor de alpha (recomendado - 1): ');
    beta = input('Insira o valor de beta (recomendado - 2): ');
    gama =  input('Insira o valor de gamma (recomendado - 0.5): ');
    delta = input('Insira o valor de delta (recomendado - 0.5): ');
    P1 = input('Insira o ponto 1 do chute inicial i.e. [x0 y0]: ');
    P2 = input('Insira o ponto 2 do chute inicial i.e. [x1 y1]: ');
    P3 = input('Insira o ponto 3 do chute inicial i.e. [x2 y2]: ');
    tolerancia = input('Insira a área desejada para o simplex final (condição de parada): '); %Condição de parada - Área < tolerancia
    
else
    disp('Testes Cadastrados: ');
    disp('(1) Booth Function');
    disp('(2) Matyas Function');
    disp('(3) Easom Function');
    disp('(4) Ackley Function');
    disp('(5) Função vista em aula');
    iusuario2 = input('Qual o teste desejado? ');
    
    if iusuario2 == 1
        f = @(x) (x(1)+2*x(2)-7)^2+(2*x(1)+x(2)-5)^2;  % Booth Function
        P1 = [8 5]; P2 = [6 2]; P3 = [7 5];
    elseif iusuario2 == 2
        f = @(x) 0.26*(x(1)^2 + x(2)^2) - 0.48*x(1)*x(2); % Matyas Function
        P1 = [0.2 5]; P2 = [1.5 2]; P3 = [3 0.4];
    elseif iusuario2 == 3
        f = @(x) -cos(x(1))*cos(x(2))*exp(-((x(1) - pi)^2 + (x(2) - pi)^2)); % Easom Function
        P1 = [0.2 5]; P2 = [1.5 2]; P3 = [3 0.4];
    elseif iusuario2 == 4
        f = @(x) -20*exp(-0.2*sqrt(0.5*(x(1)^2 + x(2)^2))) - exp(0.5*(cos(2*pi*x(1)) + cos(2*pi*x(2)))) + exp(1) + 20; % Ackley Function
        P1 = [0.2 5]; P2 = [1.5 2]; P3 = [3 0.4];
    elseif iusuario2 == 5
        f = @(x) x(1)^2 + x(2)^2 -2*x(1) - 2*x(2) + 2;  % Função Utilizada em aula
        P1 = [0.2 5]; P2 = [1.5 2]; P3 = [3 0.4];
    end
    
    %% Definição das variáveis do método
    alfa = 1; % input('Insira o valor de alpha (recomendado - 1): ')
    beta = 2; % input('Insira o valor de beta (recomendado - 2): ')
    gama = 0.5; % input('Insira o valor de gamma (recomendado - 0.5): ')
    delta = 0.5;% input('Insira o valor de delta (recomendado - 0.5): ')
    
    %% Condição de Erro
    tolerancia = 1e-5;  % input('Insira a área objetivo: ') Condição de parada - Área < tolerancia
end

Points = [P1;P2;P3];
%% Variáveis Auxiliares
ierro = 1;
iteracoes = 0;
triang_area = 0;
tarea = 0;

%% Validação dos Resultados
options = optimset('Tolfun',tolerancia,'MaxFunEvals',10000);
[x,fval,exitflag,output] = fminsearch(f,P1,options);

%% Laço Externo
while ierro == 1
    iteracoes = iteracoes + 1;
    %% Identifica os valores do Simplex
    P1 = Points(1,:);
    P2 = Points(2,:);
    P3 = Points(3,:);
    
    %% Organização dos vetores
    [Xf,iordem] = sort([f(P1);f(P2);f(P3)]);    % Coloca em ordem crescente
    x_novo = [Points(iordem(1),:);Points(iordem(2),:);Points(iordem(3),:)]; % Cria o vetor a ser atualizado
    N = length(Xf) - 1; % Obtém o n do simplex
    
    Xb = (1/N)*sum(x_novo(1:N,:),1);    % Obtém o centróide
    Xr = Xb + alfa*(Xb - x_novo(end,:));   % Reflete o pior ponto
    
    %% Comparações
    if f(x_novo(1,:))<=f(Xr) && f(Xr)<f(x_novo(2,:))
        x_novo(3,:) = Xr;
    elseif f(Xr)<f(x_novo(1,:))
        Xe = Xr + beta*(Xr - Xb);
        if f(Xe)<f(Xr)
            x_novo(3,:) = Xe;
        else
            x_novo(3,:) = Xr;
        end
    elseif f(Xr)>f(x_novo(2,:))
        Xc = x_novo(3,:) + gama*(Xb - x_novo(3,:));
        if f(Xc)<f(x_novo(3,:))
            x_novo(3,:) = Xc;
        else
            for i=1:N+1
                x_novo(i,:) = x_novo(1,:) + delta*(x_novo(i,:) - x_novo(1,:));
            end
        end
    end
    %% Verificação do Erro
    determ = [x_novo(1,:) 1;x_novo(2,:) 1;x_novo(3,:) 1];   % Monta a matriz para determinante
    triang_area = abs(0.5*det(determ)); % Cálcula A = 0.5*|D|
    tarea(iteracoes) = triang_area; % Armazena os valores da área
    
    if triang_area < tolerancia    % Se a área do triângulo é menor que a tolerância
        ierro = 0;  % Sai do While e encontrou o melhor ponto
    else    % Se não
        Points = x_novo;    % Atribui ao vetor os valores para reiniciar
    end
    
    %% Armazena os valores das iterações dos Symplex em matrizes
    symplex_x(:,iteracoes) = x_novo(:,1);
    symplex_y(:,iteracoes) = x_novo(:,2);
    symplex_z(:,iteracoes) = [f(x_novo(1,:));f(x_novo(2,:));f(x_novo(3,:))];
    
end

%% Comparação dos Resultados
disp('COMPARAÇÃO DE RESULTADOS');
disp(' ');
disp('--- Resultados da fminsearch ---');
disp(['x0 = [',num2str(x),']']);
disp(['f(x0) = ',num2str(fval)]);
disp(['Iterações: ',num2str(output.iterations)]);
disp(' ');
disp('--- Resultados do Algoritmo ---');
disp(['x0 = [',num2str(x_novo(1,:)),']']);
disp(['f(x0) = ',num2str(f(x_novo(1,:)))]);
disp(['Iterações: ',num2str(iteracoes)]);

%% Plot dos Gráficos
figure(1)
if iusuario == 1
    fa = strrep(s,'x(1)','x');
    fa = strrep(fa,'x(2)','y');
    fb = str2func(['@(x,y) ' fa]);
    fsurf(fb);title('Função Inserida'); hold on;
end

if iusuario2 == 1
    xx = -10:.2:10;
    yy = -10:.2:10;
    [X,Y] = meshgrid(xx,yy);
    Z =(X+2.*Y-7).^2+(2.*X+Y-5).^2;
elseif iusuario2 == 2
    xx = -10:.2:10;
    yy = -10:.2:10;
    [X,Y] = meshgrid(xx,yy);
    Z = 0.26.*(X.^2 + Y.^2) - 0.48.*X.*Y;
elseif iusuario2 == 3
    xx = -100:0.8:100;
    yy = -100:0.8:100;
    [X,Y] = meshgrid(xx,yy);
    Z = -cos(X).*cos(Y).*exp(-((X - pi).^2 + (Y - pi).^2));
elseif iusuario2 == 4
    xx = -5:.2:5;
    yy = -5:.2:5;
    [X,Y]=meshgrid(xx,yy);
    Z = -20.*exp(-0.2.*sqrt(0.5.*(X.^2 + Y.^2))) - exp(0.5.*(cos(2*pi.*X) + cos(2*pi.*Y))) + exp(1) + 20; % Ackley Function
elseif iusuario2 == 5
    xx=-10:.2:10;
    yy=-10:.2:10;
    [X,Y]=meshgrid(xx,yy);
    Z=X.^2+Y.^2-2.*X-2.*Y+2;
end

if iusuario ~= 1
    surf(X,Y,Z); hold on;
end

if iusuario2 == 1, title('Booth Function');
elseif iusuario2 == 2, title('Matyas Function');
elseif iusuario2 == 3, title('Easom Function');
elseif iusuario2 == 4, title('Ackley Function');
elseif iusuario2 == 5, title('Função da aula');
end

for i = 1:length(symplex_x)
    plot3(symplex_x(:,i),symplex_y(:,i),symplex_z(:,i),'c-','linewidth',1.4);hold on
end

figure(2);
% suptitle('Reprodução da Otimização');
sub1=subplot(2,1,1);
grid on; grid minor;
xlabel('Iteração');ylabel('Área do Triângulo');title('Variação da área do triângulo');
ani1 = animatedline('Color','r','linewidth',1.4);
subplot(2,1,2)
grid on;
ani2 = animatedline('Color','b');
xlabel('Coordenada x');ylabel('Coordenada y'); zlabel('Coordenada z');title('Variação dos Simplex');
view([-114 15]);
iteract = linspace(1,iteracoes,iteracoes);  % Cria um vetor de iterações (usado no plot somente)

for k=1:length(iteract)
    addpoints(ani1,iteract(k),tarea(k));
    addpoints(ani2,symplex_x(:,k),symplex_y(:,k),symplex_z(:,k));
    drawnow
    pause(0.1);
end