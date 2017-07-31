using PyPlot
theta=-pi:0.01:pi
n_antenna =24
x0 = ones(n_antenna)

function g(x, theta, d=0.5)
  n = length(x);
  af = 0
  for i = 1 : n
    af = af + x[i] .*   abs(sin(pi*d*n_antenna*(sin(theta)))./sin(pi*d*sin(theta)))/n_antenna
  end
  return 10 * log10(abs(af).^2)
end
function heaviside(theta)
   0.5 * (sign(theta) + 1)
end

function gd(theta, a=0, b=20, theta0=-3.5 * pi / 180, theta1=3.5 * pi / 180)
  -b*( heaviside(-theta+theta0) + heaviside(theta-theta1) ) - a
end

PyPlot.plot(theta,g(x0,theta))
PyPlot.plot(theta,gd(theta))
