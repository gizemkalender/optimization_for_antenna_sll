# include("./optimize.jl")
# include("./constants.jl")
import optimize
using PyPlot

theta_i=-pi;
theta_f=pi;
n=100;
delta_theta=(theta_f-theta_i)/n
theta=-pi:0.01:pi


function heaviside(theta)
   0.5 * (sign(theta) + 1)
end

# PyPlot.plot(t,y)

function gd(theta, a=0, b=20, theta0=-3.5 * pi / 180, theta1=3.5 * pi / 180)
  return -b*( heaviside(-theta+theta0) + heaviside(theta-theta1) ) - a
end
# (cos(2*i-1)*pi*d.*cos(theta))
# * abs(sin(pi*d*n_antenna*(sind(theta)))./sin(pi*d*sind(theta)))/n_antenna

# ./sin(pi*d*sind(theta)))/n_antenna

function g(x, theta, d=1.)
  n = length(x);
  k = 0:n-1
  af = sum(x.*exp(complex(0, k .* pi * d .* sin(theta))))
  out = 20 * log10(abs(af))
  return out
end


# x = ones(8)
# domain = linspace(-pi/2, pi/2, 1001)
# res = zeros(length(domain))
# for k = 1 : length(domain)
#   res[k] = g(x, domain[k])
# end
# PyPlot.plot(domain * 180 / pi, res)

function integrate(f, x, initial, final, n=1000)
  delta_theta = (final - initial) / n
  res = 0
  for i = 1 : n
    res = res + f(x, i * delta_theta) * delta_theta
  end
  return res
 end

 function error_power(x, theta)
   pow =  abs(g(x, theta) - gd(theta)).^2
   return pow
 end

 function cost(x)
   return integrate(error_power, x, theta_i, theta_f) + 1e6*(norm(x, 1)-1)^2
 end


n_antenna = 6
x0 = ones(n_antenna)
xmin, fmin, num_iter = optimize.gradient_descent(cost, optimize.estimate_gradient, x0)
