%%
x = 3.002;
y = linspace(-12*pi,+6*pi,1000);
figure();
subplot(2,1,1);
plot(y/pi, angle(Faddeeva_w(x+1i*y/2/abs(x)))); grid on;
xlabel('$y / \pi$', 'interpreter', 'latex');
title(sprintf('$w(x + i \\frac{y}{2|x|}), x = %g$', x), 'interpreter', 'latex');
ylabel('Phase Angle');

subplot(2,1,2);
semilogy(y/pi, abs(Faddeeva_w(x+1i*y/2/abs(x)))); grid on;
xlabel('$y / \pi$', 'interpreter', 'latex');
ylabel('Abs. Value');