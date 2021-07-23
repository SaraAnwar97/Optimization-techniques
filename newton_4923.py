import math
from scipy.misc import derivative


def func(x):
    return 2 * math.sin(x) - 0.1 * (math.pow(x, 2))


print("NEWTON'S METHOD")


def Newton(x_old, iteration):
    for i in range(iteration):
        first_derivative = derivative(func, x_old, dx=1e-9)
        second_derivative = derivative(func, x_old, dx=1e-6, n=2)
        x_new = x_old - (first_derivative / second_derivative)
        x_new = round(x_new, 3)
        err = abs((x_new - x_old) / x_new)
        err = round(err, 3)
        print("Iteration number =", i + 1, "\nX_old =", x_old, "\tX_new =", x_new, "\tRelative Error =", err)
        x_old = x_new
    return x_old


x_optimal = Newton(2.5, 3)
fmax = round(func(x_optimal), 3)
print("Fmax = f(", x_optimal, ")=", fmax)

print("=======================================================")
print("GOLDEN'S SECTION METHOD")


def Golden_Section(xl, xu, iteration):
    global x2
    for i in range(iteration):
        x1 = round(xl + 0.618 * (xu - xl), 4)
        x2 = round(xu - 0.618 * (xu - xl), 4)
        f1 = round(func(x1), 4)
        f2 = round(func(x2), 4)
        len = abs(round(xu - xl, 4))
        err = round((0.618 * 0.618 * len), 4)
        print(
            f'Iteration Number = {i}\nX1 = {x1},\tX2 = {x2},\tXU = {xu},\tXL = {xl},\tLength = {len},\tF(x1) = {f1},\tF(x2) = {f2},\tRelative Error = {err}')
        if f1 > f2:
            xl = x2
            x2 = x1
        else:
            xu = x1
            x1 = x2
    return x2


x2_optimal = Golden_Section(0, 4, 8)
fmax2 = round(func(x2_optimal), 4)
print("Fmax = f(", x2_optimal, ")=", fmax2)
