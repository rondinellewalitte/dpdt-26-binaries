# -*- coding: utf-8 -*-
"""TIC 232634196: a unica deteccao carregada pela alavanca de 20 anos. Segunda
rodada de revisao (item 7): em vez de descrever a ambiguidade, PUBLICAR uma
previsao falsificavel - o instante do minimo primario perto de 2028,0 sob a
efemeride quadratica e sob a linear, com as incertezas de cada uma, para que
um unico minimo medido decida.

Entrada: os 7 pontos do JSON do estado B' (3 temporadas SuperWASP com barras
de 19,8 min, 4 setores TESS com barras de 0,7-1,5 min), o P da escada e os
ajustes `cadeia_oc.ajustar` grau 1 e grau 2 (os mesmos da tabela; o JSON so
guarda chi2/dof do linear, por isso o ajuste e refeito aqui, e a sanidade
confere Q e sigma_Q contra o JSON). Referencia BTJD = BJD_TDB - 2457000.

Previsao: o ciclo E* cujo minimo linear cai mais perto de 2028,0 (ano
decimal, TDB); t_lin(E*) +- sigma pela covariancia do linear; t_quad(E*) +-
sigma pela covariancia COMPLETA do quadratico (o termo cruzado entre o
coeficiente linear e o quadratico nao e desprezivel em E* ~ 2400). A
separacao entre as duas previsoes e dada em minutos e em unidades de cada
sigma. O sigma formal do linear e reportado com o chi2_red do linear ao lado
(108/5): ele nao descreve os dados e a barra dele nao significa nada alem de
"onde a reta cairia".

EXPECTATIVA (escrita e commitada antes de rodar): E* ~ 2390-2400 (8,2 anos
de 2019,9 em 1,2511 d). Separacao quadratica - linear em 2028,0: 20-40 min
(Q E*^2 ~ -35 min, parcialmente compensado pela diferenca dos termos
lineares). sigma da previsao quadratica: 4-10 min (sigma_Q E*^2 ~ 3,4 min
mais a covariancia). A separacao e >= 3 sigma_quad: um minimo medido a 1 min
(um setor do TESS, uma noite de CCD) decide. Se a separacao sair < 2
sigma_quad, a previsao nao e um teste e a nota tem de dizer isso em vez de
publica-la como tal.
"""
import json
import sys
from pathlib import Path

import numpy as np
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cadeia_oc as co  # noqa: E402
import config  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
TIC = 232634196
ANO_ALVO = 2028.0
BTJD0 = 2457000.0


def prever(coef, cov, E):
    x = np.array([E ** k for k in range(len(coef) - 1, -1, -1)], float)   # ordem de np.vander: grau maior primeiro
    return float(x @ coef), float(np.sqrt(x @ cov @ x))


if __name__ == "__main__":
    j = json.loads((BASE / f"oc__{TIC}.json").read_text(encoding="utf-8"))
    pts = sorted(j["pontos"], key=lambda q: q["t0"])
    E = np.array([q["E"] for q in pts], float)
    t = np.array([q["t0"] for q in pts], float)
    sig = np.array([q["sig_min"] for q in pts], float) / 1440.0
    P = float(j["P_escada_d"])
    cl, rl, chi2l, covl = co.ajustar(E, t, sig, 1)
    cq, rq, chi2q, covq = co.ajustar(E, t, sig, 2)
    # sanidade: o quadratico refeito e o da tabela
    assert abs(cq[0] - j["parabola"]["Q_d_por_ciclo2"]) < 1e-15 and abs(np.sqrt(covq[0, 0]) - j["parabola"]["sQ"]) < 1e-15, (cq[0], j["parabola"])
    assert abs(chi2l - j["linear"]["chi2"]) < 1e-6 and abs(chi2q - j["parabola"]["chi2"]) < 1e-6
    print(f"TIC {TIC}: {len(E)} epocas, E de {E.min():.0f} a {E.max():.0f}; P_escada {P:.9f} d")
    print(f"  linear:     chi2 {chi2l:.1f} / {len(E) - 2} gl;  P_lin = {cl[0]:.9f} +- {np.sqrt(covl[0, 0]):.1e} d")
    print(f"  quadratico: chi2 {chi2q:.2f} / {len(E) - 3} gl;  Q = {cq[0]:.4e} +- {np.sqrt(covq[0, 0]):.2e} d/ciclo^2; dP/dt = {j['parabola']['dPdt_s_por_ano']:+.4f} +- {j['parabola']['sdPdt_s_por_ano']:.4f} s/ano")
    jd_alvo = Time(ANO_ALVO, format="decimalyear", scale="tdb").jd
    E_star = int(np.round((jd_alvo - BTJD0 - cl[1]) / cl[0]))
    t_lin, s_lin = prever(cl, covl, E_star)
    t_quad, s_quad = prever(cq, covq, E_star)
    sep_min = (t_quad - t_lin) * 1440
    s_lin_min, s_quad_min = s_lin * 1440, s_quad * 1440
    # a covariancia completa contra a diagonal, para registrar quanto o termo cruzado pesa
    s_quad_diag = float(np.sqrt(sum((E_star ** k) ** 2 * covq[i, i] for i, k in enumerate((2, 1, 0))))) * 1440
    def cal(btjd):
        return Time(btjd + BTJD0, format="jd", scale="tdb").iso[:16]
    print(f"\n== previsao para o ciclo E* = {E_star} (o minimo linear mais proximo de {ANO_ALVO} TDB)")
    print(f"  linear:     BJD_TDB {t_lin + BTJD0:.5f} ({cal(t_lin)} TDB) +- {s_lin_min:.2f} min  [sigma formal; chi2_red do linear {chi2l / (len(E) - 2):.1f}]")
    print(f"  quadratico: BJD_TDB {t_quad + BTJD0:.5f} ({cal(t_quad)} TDB) +- {s_quad_min:.2f} min  [covariancia completa; so a diagonal daria {s_quad_diag:.2f} min]")
    print(f"  separacao quadratico - linear: {sep_min:+.1f} min = {abs(sep_min) / s_quad_min:.1f} sigma_quad = {abs(sep_min) / s_lin_min:.1f} sigma_lin")
    print(f"  expectativa: E* 2390-2400; separacao 20-40 min; sigma_quad 4-10 min; separacao >= 3 sigma_quad")
    out = {"tic": TIC, "ano_alvo_tdb": ANO_ALVO, "E_star": E_star, "P_escada_d": P,
           "linear": {"BJD_TDB": t_lin + BTJD0, "sigma_min": s_lin_min, "chi2": chi2l, "dof": len(E) - 2, "iso_tdb": cal(t_lin)},
           "quadratico": {"BJD_TDB": t_quad + BTJD0, "sigma_min": s_quad_min, "sigma_min_so_diagonal": s_quad_diag, "chi2": chi2q, "dof": len(E) - 3, "iso_tdb": cal(t_quad),
                          "Q_d_por_ciclo2": cq[0], "sQ": float(np.sqrt(covq[0, 0]))},
           "separacao_min": sep_min, "separacao_em_sigma_quad": abs(sep_min) / s_quad_min, "separacao_em_sigma_lin": abs(sep_min) / s_lin_min}
    (BASE / f"previsao_{TIC}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"  -> {BASE / f'previsao_{TIC}.json'}")
