module optimize
import constants

function estimate_gradient(f, x, h=constants.DERIV_TOL)
    n = length(x)
    df = zeros(n);
    for k = 1:n
        delta = zeros(n)
        delta[k] = h
        df[k] = (f(x+delta)-f(x-delta))[1]./2./h
    end
    return df
end

function norm2(x)
    return sum(x.^2)
end

function stopping_criteria(f, k, x, df, x_old, df_old, N=constants.MAX_ITER, eps=[constants.GRAD_TOL, constants.STEP_TOL, constants.STEP_TOL])
    output = false
    if k > N
      println("No convergence in N steps")
      output=true
    elseif norm2(df) / norm2(df_old) < eps[1]^2
      output=true
    elseif norm2(x-x_old)< eps[2]^2
      println("Slow converge or no convergence. x does not change much.")
      output=true
    elseif norm2(f(x)-f(x_old))<eps[3]^2
      println("Slow converge or no convergence. No significiant cost decrease.")
      output=true
    end
    return output
end

function line_along(f, x, d)
    g(α) = (f(x + α * d))[1]
    dg(α) = ((estimate_gradient(f, x+ α*d))'*d)[1]
    return g, dg
end

function cubic_interpolation(g, dg, cubic_tol=constants.CUBIC_TOl, β=constants.BETA, s=constants.CUBIC_INT_S)
    # Determine initial interval
    a = 0
    b = s
    while ~(dg(b) >= 0 && g(b) >= g(a))
        a = b
        b = 2 * β * s
    end

    # Update  current interval
    alphabar = nothing
    while abs(a - b) > cubic_tol
        ga = g(a)
        gb = g(b)
        dga = dg(a)
        dgb = dg(b)
        z = 3 * (ga - gb) / (b - a) + dga + dgb
        w = sqrt(z^2 - dga + dgb)
        alphabar = b - (dgb + w - z) / (dgb - dga + 2 * w) * (b - a)

        if dg(alphabar) >= 0 || g(alphabar) >= g(a)
            b = alphabar
        end
        if dg(alphabar) < 0 && g(alphabar) < g(a)
            a = alphabar
        end
    end
    return alphabar
end

function gradient_descent(f, grad_f, x0)
    # Initialization
    x_old = x0
    k = 0
    x = x0
    while true
        # Current point
        dfx_old = grad_f(f, x_old)
        d = -dfx_old
        g, dg = line_along(f, x_old, d)
        alpha = cubic_interpolation(g, dg)

        # Update current point
        x = x_old + alpha * d
        dfx = grad_f(f, x)

        # Check stopping criteria
        if stopping_criteria(f, k, x, dfx, x_old, dfx_old)
            break
        end
        # callback(k, x, x_old)
        x_old = x
        k = k + 1
    end
    # callback(0, x, x_old)
    # Return results
    return x, f(x), k
end


end
