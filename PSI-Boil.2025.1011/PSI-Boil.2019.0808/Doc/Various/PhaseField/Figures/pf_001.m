
W=2.0;

x=-10:0.1:10.0;

y1 = -tanh(x/sqrt(2)/W);

y2 = -0.5*tanh(x/sqrt(2)/W)+0.5;

figure(1);
plot(x,y1,'r','LineWidth',2); 
axis([-10 10 -1.1 1.1]);
xlabel('n',      'FontSize',16,'Interpreter','latex');
ylabel('$\phi$', 'FontSize',16,'Interpreter','latex');
title ('$\phi=-\tanh \frac{n}{\sqrt{2} W}$',            'FontSize',20,'Interpreter','latex');
print -depsc 'y1.eps'

figure(2);
plot(x,y2,'g','LineWidth',2);
axis([-10 10 -1.1 1.1]);
xlabel('n',      'FontSize',16,'Interpreter','latex');
ylabel('$\phi$', 'FontSize',16,'Interpreter','latex');
title ('$\phi=-0.5\tanh \frac{n}{\sqrt{2} W} + 0.5$',   'FontSize',20,'Interpreter','latex');
print -depsc 'y2.eps'

xy1 = [x',y1'];
xy2 = [x',y2'];

save 'xy1.dat' -ascii xy1
save 'xy2.dat' -ascii xy2


