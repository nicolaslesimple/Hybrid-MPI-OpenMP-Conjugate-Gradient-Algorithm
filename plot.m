%% Strong scaling analysis with relative speed up
clc;
clear all;

x = [1,2,4];
list_1 = [3812.08,933.78,681.258,];
list_2 = [3816.17,765.292,687.069];
list_4 = [3818.68,744.149, 635.851];
% Relative speed up calculus
list_1 = list_1(1)./list_1;
list_2 = list_2(1)./list_2;
list_4 = list_4(1)./list_4;
x_big = linspace(1,16,100);

figure(1)
semilogy(x_big, interp1(x,list_1,x_big), '--r', 'LineWidth', 0.1)
hold on
semilogy(x_big, interp1(x,list_2,x_big'),'--b', 'LineWidth', 0.1)
hold on
semilogy(x_big, interp1(x,list_4,x_big), '--g', 'LineWidth', 0.1)
hold on
semilogy(x, list_1, '*r', 'LineWidth', 0.1)
hold on
semilogy(x, list_2, '*b', 'LineWidth', 0.1)
hold on
semilogy(x, list_4, '*g', 'LineWidth', 0.1)
hold on

leg1 = legend('1 OpenMP Thread','2 OpenMP Threads','4 OpenMP Threads','Location','SouthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title('Strong Scaling Analysis of the CG method : ');
set(title1,'FontName','Arial','FontSize',12)
xlabel('Number of MPI ranks','FontName','Arial','FontSize',10);
ylabel('Relative Speedup','FontName','Arial','FontSize',10);
xlim([0.5,5])
ylim([0,8])
grid on;
hold off;

% Save plot
filename='./plot/strong.eps';
print(gcf,'-depsc',filename)

%% Strong scaling analysis
clc;
clear all;

x = [1,2,4];
list_1 = [3812.08,933.78,681.258];
list_2 = [3816.17,765.292,687.069];
list_4 = [3818.68,744.149,635.851];
x_big = linspace(1,16,100);

figure(2)
semilogy(x_big, interp1(x,list_1,x_big), '--r', 'LineWidth', 0.1)
hold on
semilogy(x_big, interp1(x,list_2,x_big'),'--b', 'LineWidth', 0.1)
hold on
semilogy(x_big, interp1(x,list_4,x_big), '--g', 'LineWidth', 0.1)
hold on
semilogy(x, list_1, '*r', 'LineWidth', 0.1)
hold on
semilogy(x, list_2, '*b', 'LineWidth', 0.1)
hold on
semilogy(x, list_4, '*g', 'LineWidth', 0.1)
hold on

leg1 = legend('1 OpenMP Thread','2 OpenMP Threads','4 OpenMP Threads','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title('Strong Scaling Analysis of the CG method : ');
set(title1,'FontName','Arial','FontSize',12)
xlabel('Number of MPI ranks','FontName','Arial','FontSize',10);
ylabel('Time per CG [s]','FontName','Arial','FontSize',10);
grid on;
xlim([0.5,4.5])
hold off;

% Save the plot
filename='./plot/strong_2.eps';
print(gcf,'-depsc',filename)

%% Weak scaling analysis with relative speed up
clc;
clear all;

x = [1,2,4];
list_1 = [0.00360577,0.00274242,0.00430072];
list_2 = [0.00360708,0.00271577,0.0047284];
list_4 = [0.00360943,0.00271577,0.00487931];
list_8 = [0.00360601,0.00271577,0.00451811];
% Relative speed up calculus
list_1 = list_1(1).*x./list_1;
list_2 = list_2(1).*x./list_2;
list_4 = list_4(1).*x./list_4;
list_8 = list_8(1).*x./list_8;
x_big = linspace(1,15,1000);

figure(3)
semilogy(x_big, interp1(x,list_1,x_big), '--r', 'LineWidth', 0.1)
hold on
semilogy(x_big, interp1(x,list_2,x_big), '--b', 'LineWidth', 0.1)
hold on
semilogy(x_big, interp1(x,list_4,x_big), '--g', 'LineWidth', 0.1)
hold on
semilogy(x_big, interp1(x,list_8,x_big), '--c', 'LineWidth', 0.1)
hold on
semilogy(x, list_1, '*r', 'LineWidth', 0.1)
hold on
semilogy(x, list_2, '*b', 'LineWidth', 0.1)
hold on
semilogy(x, list_4, '*g', 'LineWidth', 0.1)
hold on
semilogy(x, list_8, '*c', 'LineWidth', 0.1)
hold on

leg1 = legend('1 OpenMP Threads','2 OpenMP Threads','4 OpenMP Threads','8 OpenMP Threads','Location','SouthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title('Weak Scaling Analysis of the CG method : ');
set(title1,'FontName','Arial','FontSize',12)
xlabel('Number of MPI Ranks','FontName','Arial','FontSize',10);
ylabel('Relative Speedup','FontName','Arial','FontSize',10);
xlim([0,5]);
grid on;
hold off;

% Save plot
filename='./plot/weak.eps';
print(gcf,'-depsc',filename)

%% Weak scaling analysis
clc;
clear all;

x = [1,2,4];
list_1 = [0.00360577,0.00274242,0.00430072];
list_2 = [0.00360708,0.00271577,0.0047284];
list_4 = [0.00360943,0.00271577,0.00487931];
list_8 = [0.00360601,0.00271577,0.00451811];
x_big = linspace(1,15,1000);

figure(4)
semilogy(x_big, interp1(x,list_1,x_big), '--r', 'LineWidth', 0.1)
hold on
semilogy(x_big, interp1(x,list_2,x_big), '--b', 'LineWidth', 0.1)
hold on
semilogy(x_big, interp1(x,list_4,x_big), '--g', 'LineWidth', 0.1)
hold on
semilogy(x_big, interp1(x,list_8,x_big), '--c', 'LineWidth', 0.1)
hold on
semilogy(x, list_1, '*r', 'LineWidth', 0.1)
hold on
semilogy(x, list_2, '*b', 'LineWidth', 0.1)
hold on
semilogy(x, list_4, '*g', 'LineWidth', 0.1)
hold on
semilogy(x, list_8, '*c', 'LineWidth', 0.1)
hold on

leg1 = legend('1 OpenMP Threads','2 OpenMP Threads','4 OpenMP Threads','8 OpenMP Threads','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title('Weak Scaling Analysis of the CG method : ');
set(title1,'FontName','Arial','FontSize',12)
xlabel('Number of MPI Ranks','FontName','Arial','FontSize',10);
ylabel('CPU time per CG iter [s]','FontName','Arial','FontSize',10);
xlim([0,5]);
grid on;
hold off;

% Save plot
filename='./plot/weak_2.eps';
print(gcf,'-depsc',filename)