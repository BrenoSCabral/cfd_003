import numpy as np
import matplotlib.pyplot as plt

# resolve a equacao de crescimento exponencial dU/dt = g*U

ttot = 300
g = 0.05 # taxa de crescimento
ci = 2 # condicao inicial
nu = 0.2 # coeficiente do filtro de Asselin (leap-frog)

dts = [1, 10]

def sol_analitica(t):
    return ci * np.exp(g * t)

for dt in dts:

    nn = int(round(ttot / dt))

    # listas com a solucao de cada esquema (comecam na condicao inicial)
    prog = [ci]
    rk2 = [ci]
    rk4 = [ci]
    mat = [ci]
    heun = [ci]
    tt = [0.0]

    # valores "correntes" de cada esquema
    u_c = ci # progressivo
    v_c = ci # rk2
    w_c = ci # rk4
    m_c = ci # matsuno
    h_c = ci # heun

    for n in range(1, nn + 1):

        # progressivo
        u_p = u_c + dt * g * u_c

        # runge-kutta 2
        kv1 = dt * g * v_c
        kv2 = dt * g * (v_c + kv1 / 2)
        v_p = v_c + kv2

        # runge-kutta 4
        kw1 = dt * g * w_c
        kw2 = dt * g * (w_c + kw1 / 2)
        kw3 = dt * g * (w_c + kw2 / 2)
        kw4 = dt * g * (w_c + kw3)
        w_p = w_c + (kw1 + 2 * kw2 + 2 * kw3 + kw4) / 6

        # matsuno: estimativa progressiva, depois usa a tendencia da estimativa
        m_est = m_c + dt * g * m_c
        m_p = m_c + dt * g * m_est

        # heun: estimativa progressiva, depois usa a media das tendencias
        h_est = h_c + dt * g * h_c
        h_p = h_c + dt / 2 * (g * h_c + g * h_est)

        # salva
        prog.append(u_p)
        rk2.append(v_p)
        rk4.append(w_p)
        mat.append(m_p)
        heun.append(h_p)
        tt.append(n * dt)

        # o futuro vira o presente
        u_c, v_c, w_c, m_c, h_c = u_p, v_p, w_p, m_p, h_p

    # leap-frog: precisa de dois tempos anteriores (n-1 e n),

    lf = [ci]
    y_ant = ci # tempo n-1 (ja filtrado)
    y_c = ci # tempo n -> primeiro passo repete a condicao inicial (como no Fortran)

    for n in range(1, nn):
        y_p = y_ant + 2 * dt * g * y_c # leap-frog
        y_c = (1 - nu) * y_c + nu * (y_ant + y_p) / 2 # filtro de Asselin no tempo n
        lf.append(y_c) # salva o valor filtrado
        y_ant = y_c
        y_c = y_p

    lf.append(y_c) # ultimo ponto (nao tem como filtrar, falta o n+1)

    # ------------------------------------------------------------------
    # graficos
    # ------------------------------------------------------------------
    tt = np.array(tt)
    ana = sol_analitica(tt)

    esquemas = {'Progressivo': prog, 'RK2': rk2, 'RK4': rk4,
                'Leap-frog': lf, 'Matsuno': mat, 'Heun': heun}

    fig, (ax_lin, ax_log, ax_err) = plt.subplots(1, 3, figsize=(12, 6))

    for ax in (ax_lin, ax_log):
        ax.plot(tt, ana, 'k:', lw=2.5, label='Analítica')

    for nome, u in esquemas.items():
        u = np.array(u)
        estilo = '--' if nome == 'Heun' else '-' # heun = rk2 neste problema
        ax_lin.plot(tt, u, estilo, label=nome)
        ax_log.plot(tt, u, estilo, label=nome)
        erro = np.abs(u - ana) / ana
        ax_err.plot(tt[1:], erro[1:], estilo, label=nome)

    # ax_lin.set_title(f'Solução (escala linear) - dt = {dt}')
    ax_lin.set_ylabel('U')

    ax_log.set_yscale('log')
    # ax_log.set_title(f'Solução (escala log) - dt = {dt}')
    ax_log.set_ylabel('U')

    ax_err.set_yscale('log')
    # ax_err.set_title(f'Erro relativo - dt = {dt}')
    ax_err.set_ylabel('|U - Ua| / Ua')

    for ax in (ax_lin, ax_log, ax_err):
        ax.set_xlabel('tempo')
        ax.grid(True, alpha=0.3)

    ax_lin.legend(fontsize=9)    # mesma legenda vale para os tres eixos

    plt.tight_layout()
    plt.savefig(f'resultado_dt{dt}.png', dpi=300)

    # tabela do erro relativo
    print(f'\ndt = {dt}')
    print(f'{"Esquema":<12} {"erro t=300":>12} {"erro medio":>12}')
    for nome, u in esquemas.items():
        erro = np.abs(np.array(u) - ana) / ana
        print(f'{nome:<12} {erro[-1]:12.2e} {erro[1:].mean():12.2e}')

# plt.show()