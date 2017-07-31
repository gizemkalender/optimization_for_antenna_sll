
import optimize

function rosenbrock(x)   # It has a global minimum of 0 at the point (1, 1).
    return (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
end

function drosenbrock(x)
    return [-2 + 2 * x[1] - 400 * x[1] * (x[2] - x[1]^2);
            200 * (x[2] - x[1]^2)]
end

function quadratic(x)  # It has a global minimum of -25/8 at the point (-5/4, 1)
    2 * x[1]^2 + x[2]^2 + 5*x[1]
end

function dquadratic(x)
    [4 * x[1] + 5;
     2 * x[2]]
end

# function print_current_state(k, x, xold)
#   if k % 100 == 0;
#     println("k: ", k, " x: ", x, " f(x):",
#     rosenbrock(x), " Δx:", sqrt(optimize.norm2(x-xold)), " Δf:", sqrt(optimize.norm2(rosenbrock(x)-rosenbrock(xold))))
#   end
# end
# function plot(k, x, xold)


x0 = [5.; 5.]
# xmin, fmin, niter = optimize.gradient_descent(rosenbrock, optimize.estimate_gradient, x0)
xmin, fmin, niter = optimize.gradient_descent(rosenbrock, optimize.estimate_gradient, x0)
