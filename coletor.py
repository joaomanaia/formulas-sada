# Coletor RGSPPDADAR
from math import acos, sin, sqrt, ceil, pi

# Regulamento
IMIN = 0.003
IMAX = 0.15
PMIN = 1.0
DMIN = 0.2
DCOM = [200, 300, 400, 500, 600, 800, 1000]


def slv(v, x=0.5):
    for _ in range(50):
        f = x - sin(x) - v
        d = 1 - sin(x + 1.57)  # cos aprox
        if abs(d) < 0.01:
            d = 0.01
        x -= f / d
        x = max(0.01, min(pi - 0.01, x))
        if abs(f) < 1e-5:
            break
    return x


# i por Vmin/Vmax (formula 4/5)
def iV(Q, K, D, V):
    ts = (8 * Q) / (D * D * V)
    t = slv(ts)
    a = Q / (K * D ** (8 / 3))
    b = t ** (2 / 3) / ts ** (5 / 3)
    return (8 * 4 ** (2 / 3) * a * b) ** 2


# i por h/D (formula 3)
def ihD(Q, K, D, hD):
    t = 2 * acos(1 - 2 * hD)
    ts = t - sin(t)
    a = Q / (K * D ** (8 / 3))
    b = t ** (2 / 3) / ts ** (5 / 3)
    return (8 * 4 ** (2 / 3) * a * b) ** 2


def arD(D):
    for d in DCOM:
        if d >= D * 1000:
            return d / 1000
    return DCOM[-1] / 1000


def dimD(K, Qm, Qa, Vm, Vi, hD):
    t = 2 * acos(1 - 2 * hD)
    ts = t - sin(t)
    # Diametro por Vmax (Q=U*A)
    Dv = sqrt((8 * Qm) / (Vm * ts))
    # Diametro por capacidade de vazao
    Dc = (
        (8 ** (5 / 3) * Qm) / (K * ts ** (5 / 3) * sqrt(IMAX) * (t / 2) ** (-2 / 3))
    ) ** (3 / 8)
    D = arD(max(Dv, Dc, DMIN))
    # imin = max(i por h/D com Qmax, i por Vmin com Qal, IMIN reg)
    i_hD = ihD(Qm, K, D, hD)  # Qmax para h/D
    i_Vm = iV(Qa, K, D, Vi)  # Qal para Vmin
    im = max(i_hD, i_Vm, IMIN)
    # imax = min(i por Vmax, IMAX reg)
    ix = min(iV(Qm, K, D, Vm), IMAX)
    return D, im, ix


def impl(nc, Ls, its, p0, im, ix):
    ps = [p0]
    ics = []
    for k in range(nc - 1):
        L = Ls[k]
        it = its[k]  # + = desce, - = sobe
        dt = -it * L  # desce = cota diminui = delta negativo
        pi_ = ps[-1] + im * L + dt
        if pi_ >= PMIN:
            ps.append(pi_)
            ics.append(im)
        else:
            in_ = (PMIN - ps[-1] - dt) / L
            if in_ <= ix:
                ps.append(PMIN)
                ics.append(max(in_, im))
            else:
                px = ps[-1] + ix * L + dt
                ps.append(max(px, PMIN))
                ics.append(ix)
    return ps, ics


# === MAIN ===
K = float(input("Ks="))
Qm = float(input("Qmax(l/s)=")) / 1000
Qa = float(input("Qal(l/s)=")) / 1000
Vm = float(input("Vmax(m/s)="))
Vi = float(input("Vmin(m/s)="))
hD = float(input("h/D="))
p0 = float(input("Prof CV1(m)="))
n = int(input("N de CVs="))

Ls = []
its = []
for k in range(n - 1):
    print("-" * 14)
    print("Trecho " + str(k + 1) + ":")
    Ls.append(float(input(" L(m)=")))
    its.append(float(input(" i terr(%)=")))

Ls = [L for L in Ls]
its = [i / 100 for i in its]

D, im, ix = dimD(K, Qm, Qa, Vm, Vi, hD)
ps, ics = impl(n, Ls, its, p0, im, ix)

print()
print("=" * 14)
print("* RESULTADO *")
print("=" * 14)
print("D = {:.2f} mm".format(D * 1000))
print("i min = {:.2f}%".format(im * 100))
print("i max = {:.2f}%".format(ix * 100))
print()
print("-- TRECHOS --")
for k in range(n - 1):
    print("T{}: i={:.2f}%".format(k + 1, ics[k] * 100))
print()
print("-- PROF CVs --")
for k in range(len(ps)):
    print("CV{}: {:.2f}m".format(k + 1, ps[k]))
