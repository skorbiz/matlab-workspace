function plot_sol(T, u)
% Generate 2 plots vertically

subplot(2,1,1);
% Plot x(1) and x(2)
    plot(T,u(:,1),'r', T,u(:,2),'g')

    xlabel('t'); ylabel('x');
    title('Coordinates')
    legend('x(1)', 'x(2)')

subplot (2,1,2);
% Plot xd(1) and xd(2)
    plot(T,u(:,3),'r', T,u(:,4),'g')

    xlabel('t'); ylabel('xd');
    title ('Velocities')
    legend('xd(1)', 'xd(2)')
