def t_xy_xy(l, m, n, Vdds, Vddp, Vddd):
    return (3 * l**2 * m**2 * Vdds
            + (l**2 + m**2 - 4 * l**2 * m**2) * Vddp
            + (n**2 + l**2 * m**2) * Vddd)

def t_xy_yz(l, m, n, Vdds, Vddp, Vddd):
    return (3 * l * m**2 * n * Vdds
            + l * n * (1 - 4 * m**2) * Vddp
            + l * n * (m**2 - 1) * Vddd)

def t_xy_zx(l, m, n, Vdds, Vddp, Vddd):
    return (3 * l**2 * m * n * Vdds
            + m * n * (1 - 4 * l**2) * Vddp
            + m * n * (l**2 - 1) * Vddd)

# yz band: one cyclic permutation l->m, m->n, n->l
def t_yz_yz(l, m, n, Vdds, Vddp, Vddd):
    return t_xy_xy(m, n, l, Vdds, Vddp, Vddd)

def t_yz_zx(l, m, n, Vdds, Vddp, Vddd):
    return t_xy_yz(m, n, l, Vdds, Vddp, Vddd)

def t_yz_xy(l, m, n, Vdds, Vddp, Vddd):
    return t_xy_zx(m, n, l, Vdds, Vddp, Vddd)

# xz band: two cyclic permutations l->n, m->l, n->m
def t_xz_xz(l, m, n, Vdds, Vddp, Vddd):
    return t_xy_xy(n, l, m, Vdds, Vddp, Vddd)

def t_xz_xy(l, m, n, Vdds, Vddp, Vddd):
    return t_xy_yz(n, l, m, Vdds, Vddp, Vddd)

def t_xz_yz(l, m, n, Vdds, Vddp, Vddd):
    return t_xy_zx(n, l, m, Vdds, Vddp, Vddd)


if __name__ == '__main__':
    from sympy import symbols, simplify, sqrt

    Vdds, Vddp, Vddd = symbols('V_dds V_ddp V_ddd')

    bands = [
        ('xy', [(t_xy_xy, 'xy-xy'), (t_xy_yz, 'xy-yz'), (t_xy_zx, 'xy-zx')]),
        ('yz', [(t_yz_yz, 'yz-yz'), (t_yz_zx, 'yz-zx'), (t_yz_xy, 'yz-xy')]),
        ('xz', [(t_xz_xz, 'xz-xz'), (t_xz_xy, 'xz-xy'), (t_xz_yz, 'xz-yz')]),
    ]

    for l, m, n, label in [
        (1, 0, 0, 'x'),
        (0, 1, 0, 'y'),
        (1/sqrt(2), 1/sqrt(2), 0, 'x+y'),
        (1/sqrt(2), -1/sqrt(2), 0, 'x-y'),
    ]:
        print(f"\n=== direction {label} ===")
        for band_name, funcs in bands:
            print(f"  {band_name} band:")
            for fn, name in funcs:
                val = simplify(fn(l, m, n, Vdds, Vddp, Vddd))
                print(f"    t_{name} = {val}")