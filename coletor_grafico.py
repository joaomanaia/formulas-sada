# Coletor RGSPPDADAR
from math import acos, sin, sqrt, ceil, pi
from casioplot import set_pixel, draw_string, show_screen, clear_screen

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


# === DESENHO COM CASIOPLOT ===
def draw_line(x1, y1, x2, y2, color=(0, 0, 0)):
    """Desenha uma linha entre dois pontos usando Bresenham."""
    dx = abs(x2 - x1)
    dy = abs(y2 - y1)
    sx = 1 if x1 < x2 else -1
    sy = 1 if y1 < y2 else -1
    err = dx - dy
    while True:
        set_pixel(x1, y1, color)
        if x1 == x2 and y1 == y2:
            break
        e2 = 2 * err
        if e2 > -dy:
            err -= dy
            x1 += sx
        if e2 < dx:
            err += dx
            y1 += sy


def draw_result():
    clear_screen()

    # Dimensoes da tela Casio fx-CG50: 384x192
    W, H = 384, 192
    MARGIN = 30

    # Cores
    BLACK = (0, 0, 0)
    BLUE = (0, 0, 255)
    RED = (255, 0, 0)
    GREEN = (0, 128, 0)

    # Titulo
    draw_string(5, 2, "Coletor RGSPPDADAR", BLACK, "small")
    draw_string(5, 14, "D={}mm".format(int(D * 1000)), BLUE, "small")

    # Calcular posicoes X para cada CV
    total_L = sum(Ls) if Ls else 1
    xs = [MARGIN]
    for L in Ls:
        xs.append(xs[-1] + int((W - 2 * MARGIN) * L / total_L))

    # Calcular cotas do terreno (i positivo = desce = cota diminui)
    # Assumimos cota inicial do terreno = 100 (valor arbitrario)
    cotas_terr = [100.0]
    for k in range(len(Ls)):
        # i positivo = terreno desce = delta negativo na cota
        delta = -its[k] * Ls[k]
        cotas_terr.append(cotas_terr[-1] + delta)

    # Calcular cotas do fundo do coletor (cota terreno - profundidade)
    cotas_col = [cotas_terr[k] - ps[k] for k in range(len(ps))]

    # Encontrar min/max de todas as cotas para escala Y
    all_cotas = cotas_terr + cotas_col
    c_max = max(all_cotas)
    c_min = min(all_cotas)
    c_range = c_max - c_min if c_max > c_min else 1

    # Area de desenho
    Y_TOP = 35
    Y_BOT = H - 25
    Y_RANGE = Y_BOT - Y_TOP

    def cota_to_y(c):
        # Cota maior = mais alto na tela (Y menor)
        return int(Y_BOT - ((c - c_min) / c_range) * Y_RANGE)

    # Desenhar linha do terreno
    for k in range(len(cotas_terr) - 1):
        x1, y1 = xs[k], cota_to_y(cotas_terr[k])
        x2, y2 = xs[k + 1], cota_to_y(cotas_terr[k + 1])
        draw_line(x1, y1, x2, y2, GREEN)
    draw_string(W - MARGIN + 2, cota_to_y(cotas_terr[-1]) - 5, "Terr", GREEN, "small")

    # Desenhar perfil do coletor e CVs
    for k in range(len(ps)):
        x = xs[k]
        yt = cota_to_y(cotas_terr[k])  # topo (terreno)
        yc = cota_to_y(cotas_col[k])  # fundo (coletor)

        # Desenhar linha vertical da CV (poco)
        draw_line(x, yt, x, yc, BLUE)

        # Desenhar ponto da CV
        for dx in range(-2, 3):
            for dy in range(-2, 3):
                set_pixel(x + dx, yc + dy, RED)

        # Label da CV
        draw_string(x - 8, yc + 5, "CV{}".format(k + 1), BLACK, "small")
        draw_string(x - 10, yc + 15, "{:.1f}m".format(ps[k]), BLACK, "small")

    # Desenhar trechos do coletor (linhas entre CVs)
    for k in range(len(ps) - 1):
        x1 = xs[k]
        x2 = xs[k + 1]
        y1 = cota_to_y(cotas_col[k])
        y2 = cota_to_y(cotas_col[k + 1])
        draw_line(x1, y1, x2, y2, BLUE)

        # Mostrar inclinacao no meio do trecho
        xm = (x1 + x2) // 2
        ym = (y1 + y2) // 2 - 8
        draw_string(xm - 15, ym, "i={:.1f}%".format(ics[k] * 100), RED, "small")

    show_screen()


draw_result()
