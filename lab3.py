import numpy as np
import scipy.integrate as spi
from sympy import limit, symbols, sqrt, ln
from sympy.calculus.singularities import singularities


def rectangles(func, a, b, n, left=False, right=False, middle=False):
    h = (b - a) / n
    x = a
    integral = 0

    if left:
        for i in range(n):
            integral += func(x)
            x += h
    elif right:
        x = a + h
        for i in range(n):
            integral += func(x)
            x += h
    elif middle:
        x = a + h / 2
        for i in range(n - 1):
            integral += func(x)
            x += h

    return integral * h


def trapezoids(func, a, b, n, tr=True):
    h = (b - a) / n
    x = a
    integral = 0

    integral += func(x) / 2
    x += h

    for i in range(n - 1):
        integral += func(x)
        x += h

    integral += func(x) / 2

    return integral * h


def simpson(func, a, b, n, tr=True):
    h = (b - a) / n
    x = a
    integral = 0

    integral += func(x)
    x += h

    for i in range(1, n, 2):
        integral += 4 * func(x)
        x += 2 * h

    x = a + 2 * h
    for i in range(2, n - 1, 2):
        integral += 2 * func(x)
        x += 2 * h

    integral += func(b)

    return integral * h / 3


def rung_rule(func, method, a, b, eps, k, n=4, **kwargs):
    integral_n = method(func, a, b, n, kwargs)
    integral_2n = method(func, a, b, 2 * n, kwargs)
    print(integral_n, integral_2n, n)
    if abs(integral_n - integral_2n) / (2 ** k - 1) < eps:
        return integral_n, n

    while abs(integral_n - integral_2n) / (2 ** k - 1) >= eps:
        n *= 2
        integral_n = integral_2n
        integral_2n = method(func, a, b, 2 * n, kwargs)
        print(integral_n, integral_2n, n, abs(integral_n - integral_2n) / (2 ** k - 1))

    return integral_2n, n


def func_int(func, a, b):
    return spi.quad(func, a, b)


def int_convergence(func, dot, a, b):
    x = symbols('x')
    if dot < a or dot > b:
        return True
    elif dot == a:
        if (limit(func, x, a, "+") is not None
                and limit(func, x, a, "+") != np.inf and limit(func, x, a, "+") != -1 * np.inf):
            return True
    elif dot == b:
        if (limit(func, x, b, "-") is not None
                and limit(func, x, b, "-") != np.inf and limit(func, x, b, "-") != -1 * np.inf):
            return True
    elif (limit(func, x, dot, "-") is not None and limit(func, x, dot, "+") is not None
          and limit(func, x, dot, "-") != np.inf and limit(func, x, dot, "-") != -1 * np.inf
          and limit(func, x, dot, "+") != np.inf and limit(func, x, dot, "+") != -1 * np.inf):
        return True
    else:
        return False


def main():
    print("Выберите функцию для интегрирования:")
    print("1. x^2")
    print("2. x^3")
    print("3. sin(x)")
    print("4. x^3 - 5x^2 + 3x - 16")
    print("Несобственные интегралы второго рода:")
    print("5. 1/sqrt(x), 0")
    print("6. 1/(1 - x), 1")
    print("7. 1 / sqrt(2 * x - x^2)")
    choice = int(input())

    if choice == 1:
        func = lambda x: x ** 2
    elif choice == 2:
        func = lambda x: x ** 3
    elif choice == 3:
        func = np.sin
    elif choice == 4:
        func = lambda x: x ** 3 - 5 * x ** 2 + 3 * x - 16
    elif choice == 5:
        func = lambda x: 1 / sqrt(x)
        # dot = 0
        x = symbols('x')
        F = 2 * sqrt(x)
        f = 1 / sqrt(x)
    elif choice == 6:
        func = lambda x: 1 / (1 - x)
        # dot = 1
        x = symbols('x')
        F = -1 * ln(abs(1 - x))
        f = 1 / (1 - x)
    elif choice == 7:
        func = lambda x: 1 / sqrt(2 * x - x ** 2)
        # dot = 1
        x = symbols('x')
        F = -1 * ln(abs(1 - x))
        f = 1 / sqrt(2 * x - x ** 2)
    else:
        print("Неверный выбор функции.")
        return None
    if func is None:
        return

    a = float(input("Введите нижний предел интегрирования: "))
    b = float(input("Введите верхний предел интегрирования: "))
    eps = float(input("Введите точность вычислений: "))

    if choice == 5 or choice == 6:
        x = symbols('x')
        dots = singularities(f, x)
        for dot in dots:
            if not int_convergence(F, dot, a, b):
                print("Интеграл не существует")
                exit()
            else:
                print("Интеграл сходится")

    print("\nВыберите метод интегрирования:")
    print("1. Метод прямоугольников (левые)")
    print("2. Метод прямоугольников (правые)")
    print("3. Метод прямоугольников (средние)")
    print("4. Метод трапеций")
    print("5. Метод Симпсона")
    choice_m = int(input())

    if choice != 5 and choice != 6 and choice != 7:
        if choice_m == 1:
            integral, n = rung_rule(func, rectangles, a, b, eps, 1, left=True)
        elif choice_m == 2:
            integral, n = rung_rule(func, rectangles, a, b, eps, 1, right=True)
        elif choice_m == 3:
            integral, n = rung_rule(func, rectangles, a, b, eps, 2, middle=True)
        elif choice_m == 4:
            integral, n = rung_rule(func, trapezoids, a, b, eps, 2, tr=True)
        elif choice_m == 5:
            integral, n = rung_rule(func, simpson, a, b, eps, 4, tr=True)
        else:
            print("Неверный выбор метода.")
            return

        # int_count, error = func_int(func, a, b)

        print(f"\nЗначение интеграла: {integral}")
        print(f"Число разбиения интервала интегрирования для достижения требуемой точности: {n}")
        # print(f"Погрешность: {abs(int_count - integral)}")

    else:
        dots_new = []
        dots = singularities(f, x)
        for d in dots:
            if d > a and d < b:
                dots_new.append(d)

        if len(dots_new) > 1:
            intervals = [[a, dots_new[0]]]
            for i in range(1, len(dots_new)):
                intervals.append([dots_new[i - 1], dots_new[i]])
            intervals.append([dots_new[-1], b])
        elif len(dots_new) == 0:
            intervals = [[a, b]]
        elif len(dots_new) == 1:
            intervals = [[a, dots_new[0]], [dots_new[0], b]]
        integrall = 0
        n = 0
        for i in intervals:
            a_i = i[0]
            b_i = i[1]
            if i[0] in dots:
                a_i += 0.000001
            if i[1] in dots:
                b_i -= 0.000001
            if choice_m == 1:
                integral_i, n_i = rung_rule(func, rectangles, a_i, b_i, eps, 1, left=True)
            elif choice_m == 2:
                integral_i, n_i = rung_rule(func, rectangles, a_i, b_i, eps, 1, right=True)
            elif choice_m == 3:
                integral_i, n_i = rung_rule(func, rectangles, a_i, b_i, eps, 2, middle=True)
            elif choice_m == 4:
                integral_i, n_i = rung_rule(func, trapezoids, a_i, b_i, eps, 2, tr=True)
            elif choice_m == 5:
                integral_i, n_i = rung_rule(func, simpson, a_i, b_i, eps, 4, tr=True)
            integrall += integral_i
            n += n_i

        print(f"\nЗначение интеграла: {integrall}")
        print(f"Число разбиения интервала интегрирования для достижения требуемой точности: {n}")





if __name__ == "__main__":
    main()
