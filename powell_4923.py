import numpy as np
from sympy import *

def powell_method(X1,x,y,lamda,func):
    S = np.array([1, 0])
    X = np.zeros((5 + 1, 2))
    e = 0.01
    X[0] = X1
    l = np.zeros(5)
    for i in range(0,5):
            print("Iteration #", i+1 , ":")
            if i != 3:
                if S[1] == 0:
                    S = np.array([0, 1])
                else:
                    S = np.array([1, 0])
            else:
                S = X[3] - X[1]
            print("S", S)
            f = func.subs({x: X[i, 0], y: X[i, 1]})
            fp = func.subs({x: X[i, 0] + e * S[0], y: X[i, 1] + e * S[1]})
            fn = func.subs({x: X[i, 0] - e * S[0], y: X[i, 1] - e * S[1]})
            print(f'f_current_value = {f}\nf_positive = {fp}\nf_negative = {fn}')
            if fn < f:
                S *= -1
            fl = func.subs({x: X[i, 0] + S[0] * lamda, y: X[i, 1] + S[1] * lamda})
            dfdl = Derivative(fl, lamda).doit()
            print("Derivative of lambda =", dfdl)
            l[i] = solve(dfdl, lamda)[0]
            print("lambda =", l[i])
            X[i + 1] = X[i] + l[i] * S
            print("X(",i + 2, ")=", X[i + 1], "\n----------------------------------------------")

def steepest_descent_method(X1,x,y,lamda,func):
    dfdx = Derivative(func, x).doit()
    dfdy = Derivative(func, y).doit()
    grad_f = np.zeros((4, 2))
    X = X1
    l = np.zeros(4)
    for i in range(0,4):
        print("Iteration #", i + 1 ,":")
        grad_f[i, 0] = dfdx.subs({x: X[0], y: X[1]})
        grad_f[i, 1] = dfdy.subs({x: X[0], y: X[1]})
        print("grad_f", i + 1, " = ", grad_f[i])
        fl = func.subs({x: X[0] - grad_f[i, 0] * lamda, y: X[1] - grad_f[i, 1] * lamda})
        dfdl = Derivative(fl, lamda).doit()
        print("Derivative of lamdba =", dfdl)
        l[i] = solve(dfdl, lamda)[0]
        print("lambda =", l[i])
        X = X - l[i] * grad_f[i]
        print("X(", i + 1, ")=", X, "\n----------------------------------------------")

def main():
    x = symbols('x')
    y = symbols('y')
    lamda = symbols('lamda')
    func = x - y + 2 * pow(x, 2) + 2 * x * y + pow(y, 2)
    X1 = [0, 0]
    print("----------------------------------------------POWELL'S METHOD----------------------------------------------")
    powell_method(X1,x,y,lamda,func)
    print("----------------------------------------------STEEPEST DESCENT'S METHOD----------------------------------------------")
    steepest_descent_method(X1,x,y,lamda,func)



if __name__ == "__main__":
    main()
