function [] = GRAPHE()
%% GRAPHE DE LA DISTRIBUTION D'ELECTRON
    A = load('COORDONNEES.dat');
    x(1:2:2*size(A,1)) = A(:,1);
    x(2:2:2*size(A,1)) = A(:,4);
    y(1:2:2*size(A,1)) = A(:,2);
    y(2:2:2*size(A,1)) = A(:,5);
    z(1:2:2*size(A,1)) = A(:,3);
    z(2:2:2*size(A,1)) = A(:,6);
    C = repmat([255,0,0;0,0,255],size(A,1),1);
    scatter3(x',y',z',1,C,'filled')
    hold on
    grid on
    xlim([floor(min(x)), ceil(max(x))]);
    ylim([floor(min(y)), ceil(max(y))]);
    zlim([floor(min(z)), ceil(max(z))]);
    set(gca,'FontName','Times New Roman','FontSize',16);
    xlabel('x (Å)','FontName','Times New Roman','FontSize',16)
    ylabel('y (Å)','FontName','Times New Roman','FontSize',16)
    zlabel('z (Å)','FontName','Times New Roman','FontSize',16)

%% GRAPHE DE L'ENERGIE EN FONCTION DE LA DISTANCE
    figure
    A = load('MORSE.dat');
    plot(A(:,1),A(:,2))
    hold on
    set(gca,'FontName','Times New Roman','FontSize',16);
    xlabel('S (Å)','FontName','Times New Roman','FontSize',16)
    ylabel('$U(S)$ (eV)','Interpreter','LaTeX','FontSize',16)
    ylim([-50 300]);

%% GRAPHE DE LE SOLUTION DE L'EQUATION TRANSCENDANTE
    figure
    A = load('SOLUTION.dat');
    plot(A(:,1),A(:,2))
    hold on
    set(gca,'FontName','Times New Roman','FontSize',16);
    xlabel('S (Å)','FontName','Times New Roman','FontSize',16)
    ylabel('$a(S)$','Interpreter','LaTeX','FontSize',16)
end