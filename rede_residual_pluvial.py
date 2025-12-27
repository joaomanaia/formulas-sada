import math


def dx(fn, x, delta=0.001):
    return (fn(x + delta) - fn(x)) / delta


def solve(fn, value, x=0.5, maxtries=1000, maxerr=0.00001):
    for tries in range(maxtries):
        err = fn(x) - value
        if abs(err) < maxerr:
            return x
        slope = dx(fn, x)
        if abs(slope) < 1e-10:
            if err > 0:
                x = x + 0.01
            else:
                x = x - 0.01
        else:
            x = x - err / slope
    return x


while True:
    print("1: Q = U*A")
    print("2: Cap. Vazao")
    print("3: h/D (i)")
    print("4: Vmin (i)")
    print("5: Vmax (i)")
    formula = input("Formula: ")

    if formula == "1":
        Qmax = float(input("Qmax (m3/s): "))
        Vmax = float(input("Vmax (m/s): "))
        hDivD = float(input("h/D: "))

        Theta = 2 * math.acos(1 - 2 * hDivD)
        ThetaMinusSin = Theta - math.sin(Theta)

        D = pow(((8 * Qmax) / (Vmax * ThetaMinusSin)), 1 / 2)

        print("theta = {:.2f}rad".format(Theta))
        print("theta-sin(theta) = {:.2f} rad".format(ThetaMinusSin))
        print("D >= {:.4f}m ({:.2f}mm)".format(D, D * 1000))

    elif formula == "2":
        Qmax = float(input("Qmax (m3/s): "))
        Ks = float(input("Ks: "))
        hDivD = float(input("h/D: "))
        imax = float(input("imax: "))

        Theta = 2 * math.acos(1 - 2 * hDivD)
        ThetaMinusSin = Theta - math.sin(Theta)

        D = pow(
            (pow(8, 5 / 3) * Qmax)
            / (
                Ks
                * pow(ThetaMinusSin, 5 / 3)
                * math.sqrt(imax)
                * pow(Theta / 2, -2 / 3)
            ),
            3 / 8,
        )

        print("theta = {:.2f} rad".format(Theta))
        print("theta-sin(theta) = {:.2f} rad".format(ThetaMinusSin))
        print("D >= {:.4f} m ({:.2f} mm)".format(D, D * 1000))
    elif formula == "3":
        Qmax = float(input("Qmax (m3/s): "))
        Ks = float(input("Ks: "))
        D = float(input("D (m): "))
        hDivD = float(input("h/D: "))

        Theta = 2 * math.acos(1 - 2 * hDivD)
        ThetaMinusSin = Theta - math.sin(Theta)

        a = Qmax / (Ks * pow(D, 8 / 3))
        b = pow(Theta, 2 / 3) / pow(ThetaMinusSin, 5 / 3)

        i = pow(8 * pow(4, 2 / 3) * a * b, 2)

        print("theta = {:.2f} rad".format(Theta))
        print("theta-sin(theta) = {:.2f} rad".format(ThetaMinusSin))
        print("i >= {:.4f} m/m\ni = {:.2f} %".format(i, i * 100))
    elif formula == "4":
        Qal = float(input("Qal (m3/s): "))
        Ks = float(input("Ks: "))
        D = float(input("D (m): "))
        Vmin = float(input("Vmin (m/s): "))

        ThetaMinusSin = (8 * Qal) / (pow(D, 2) * Vmin)

        def fn(Theta):
            return Theta - math.sin(Theta)

        Theta = solve(fn, ThetaMinusSin)

        a = Qal / (Ks * pow(D, 8 / 3))
        b = pow(Theta, 2 / 3) / pow(ThetaMinusSin, 5 / 3)

        i = pow(8 * pow(4, 2 / 3) * a * b, 2)

        print("theta = {:.2f} rad".format(Theta))
        print("theta-sin(theta) = {:.2f} rad".format(ThetaMinusSin))
        print("i >= {:.4f} m/m\ni = {:.2f} %".format(i, i * 100))

    elif formula == "5":
        Qmax = float(input("Qmax (m3/s): "))
        Ks = float(input("Ks: "))
        D = float(input("D (m): "))
        Vmax = float(input("Vmax (m/s): "))

        ThetaMinusSin = (8 * Qmax) / (pow(D, 2) * Vmax)

        def fn(Theta):
            return Theta - math.sin(Theta)

        Theta = solve(fn, ThetaMinusSin)

        a = Qmax / (Ks * pow(D, 8 / 3))
        b = pow(Theta, 2 / 3) / pow(ThetaMinusSin, 5 / 3)

        i = pow(8 * pow(4, 2 / 3) * a * b, 2)

        print("theta = {:.2f} rad".format(Theta))
        print("theta-sin(theta) = {:.2f} rad".format(ThetaMinusSin))
        print("i <= {:.4f} m/m\ni = {:.2f} %".format(i, i * 100))
    else:
        print("Formula nao reconhecida")

    r = input("Voltar? (s/n): ")
    if r.lower() != "s":
        break
