# -*- coding: utf-8 -*-
"""Rodada seis, Passo D: sensibilidade dos 26 ao vies arquival da Seccao 3.4.

O desenho nao muda: as epocas ARMAZENADAS dos `pontos` de cada registro. A
cadeia gravou cada temporada SuperWASP com t0 = t0_bruto - VIES_CADEIA_MIN/1440
(VIES_CADEIA_MIN = -2,14 min: a epoca e adiantada em 2,14 min) e barra
hypot(max(formal, dispersao), VIES_CADEIA_SIG = 0,78). Aqui a correcao e
refeita com delta em {0, -2,14, -4,28} min (t0 = bruto - delta/1440) e a
barra do vies em {0,78, 2,14} min (barra = hypot(interna, sigma_v), interna =
sqrt(barra^2 - 0,78^2)); as epocas TESS nao mudam. Seis combinacoes; a atual
e (delta -2,14, sigma_v 0,78). Cada combinacao passa por (a) a cadeia
DIAGONAL (`reinflar_tess_26.refazer`: `cadeia_oc.ajustar`, p_curv por
Delta chi2, `cadeia_oc.adversarial`, p_gof, vao) e (b) o GP AGRUPADO
Matern-3/2 do Passo B (`deriva_gp_26.refazer_gp`, A e tau_c de
deriva_gp_26.json: C = diag + K nas TESS, GLS, adversarial marginal na
metrica de C). Por alvo e combinacao: vereditos (curvatura, adversarial,
p_gof), dP/dt, sigma, e Delta(dP/dt)/sigma contra a combinacao atual do
mesmo metodo. Resumo: mudancas no D11 e no nucleo {198408416, 232634196,
329246824, 392536812}. Saidas: sensib_vies_26.parquet, sensib_vies_26.json.

EXPECTATIVA (escrita e commitada ANTES de rodar):
  - Um deslocamento comum delta do bloco SuperWASP muda Q por ~delta/DeltaE^2:
    Delta(dP/dt)/sigma por passo de 2,14 min - diagonal: mediana 0,3-0,8, max
    <= 1,5 (232634196 ~0,1: barras de 19,9 min); GP: menor (sigma maior).
  - Vereditos com sigma_v = 0,78: delta = 0 e delta = -4,28 mudam no maximo 1
    alvo do D11 cada (diagonal 10-11 de 11); o nucleo de 4 nao muda em
    nenhuma das seis combinacoes (curvatura a >= 8 sigma, barras SW 3-20 min).
  - sigma_v = 2,14: as barras SW de 1,1-1,8 min sobem x1,3-2 (230386284,
    229914020, 377253090, 390021728, 229461186, 237116051, 199688409,
    126945917, 424461577); D11 diagonal cai para 7-10; GP agrupado (7 hoje)
    para 5-7; nenhuma deteccao nova; p_gof < 0,05 de 3 para 1-2.
  - A mediana de |Delta(dP/dt)| entre delta = 0 e -4,28 e menor que a barra em
    todos os 26 exceto 0-3 (os de barra SW < 1,5 min).
Desvio em qualquer direcao vai ao relatorio.

RESULTADO (2026-09-17, contra a expectativa de a372dec): diagonal, por passo
de 2,14 min: |Delta dP/dt|/sigma mediana 0,68, max 2,32 (230386284; esperado
<= 1,5), > 1 sigma em 7 de 26; GP: mediana 0,36, max 0,56. Vereditos,
sigma_v 0,78: delta 0 -> D11 11 de 11; delta -4,28 -> 9 (caem 377253090 e
390021728, entra 115244268: 2 mudancas, esperado <= 1). sigma_v 2,14: 9 de
11 (os mesmos dois caem), sigma(dP/dt) so x1,09, p_gof flags iguais (3;
esperado 1-2). Nucleo intacto nas seis combinacoes diagonais. GP agrupado:
delta 0 -> 6 (198408416 CAI, p_adv passa de 0,03 a > 0,05: o nucleo QUEBRA
no GP com delta 0), -2,14 -> 7, -4,28 -> 8 (entra 230386284); sigma_v 2,14
nao muda o GP. Excursao de dP/dt entre delta 0 e -4,28 (diagonal): mediana
1,36 sigma, > 1 sigma em 16 de 26 (esperado 0-3): a incerteza do vies e um
termo de primeira ordem no VALOR de dP/dt da maioria, nao nos vereditos.

ESTADO E (2026-09-17, rodada sete): a constante da cadeia passou a -1,26 +- 0,75
(vies_hj_indep). ATUAL e lido da cadeia; a grade ganha -1,26 (atual) e -1,56
(ExoClock III, so como sensibilidade dos vereditos); sigma_v em {0,75, 2,14}.
EXPECTATIVA (antes de rodar; detalhe em comparar_estados_DE.py): ATUAL reproduz
a Tabela 3 do estado E (assercao); as linhas com sigma_v 2,14 nos deltas 0 /
-2,14 / -4,28 reproduzem o estado D ao digito (mesmos pontos brutos); ExoClock
-1,56: D11 11/11, nucleo intacto na diagonal, GP 7 (198408416 p_adv ~0,038);
delta 0: 11/11; -4,28: 9/11; GP 6 / 7 / 8 em 0 / -2,14 / -4,28.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cadeia_oc as co  # noqa: E402
import config  # noqa: E402
import deriva_gp_26 as GP  # noqa: E402
import reinflar_tess_26 as RT  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
# ESTADO E: a combinacao ATUAL vem da cadeia (-1,26 / 0,75); -1,56 (ExoClock III) entra so como sensibilidade dos
# vereditos; 0 / -2,14 / -4,28 ficam para comparar com o estado D (as linhas com sigma_v 2,14 tem de reproduzir o
# estado D ao digito: os pontos brutos sao os mesmos).
ATUAL = (co.VIES_CADEIA_MIN, co.VIES_CADEIA_SIG)
DELTAS = tuple(sorted({0.0, -1.26, -1.56, -2.14, -4.28, ATUAL[0]}, reverse=True))
SIGS = (ATUAL[1], 2.14)
NUCLEO = (198408416, 232634196, 329246824, 392536812)
P_DET = 0.05
EXPECTATIVA = "a372dec (estado D); estado E: ver comparar_estados_DE.py"


def pontos_com(pontos, delta, sig_v):
    p = pontos.copy()
    sw = ~p.fonte.str.startswith("TESS").values
    bruto = p.t0.values[sw] + co.VIES_CADEIA_MIN / 1440.0
    p.loc[sw, "t0"] = bruto - delta / 1440.0
    interna = np.sqrt(np.maximum(p.sig_min.values[sw] ** 2 - co.VIES_CADEIA_SIG ** 2, 1e-6))
    p.loc[sw, "sig_min"] = np.hypot(interna, sig_v)
    return p


if __name__ == "__main__":
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    tics = sorted(int(t) for t in t26.index)
    D11 = {t for t in tics if t26.loc[t, "curvatura"]}
    hip = json.loads((BASE / "deriva_gp_26.json").read_text(encoding="utf-8"))["hiper_agrupados"]["matern32"]
    A, tc = hip["A_min"], hip["tau_c_d"]
    linhas = []
    for tic in tics:
        j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
        pontos = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
        P = float(j["P_escada_d"])
        for delta in DELTAS:
            for sig_v in SIGS:
                p = pontos_com(pontos, delta, sig_v)
                if (delta, sig_v) == ATUAL:
                    assert np.allclose(p.t0.values, pontos.t0.values) and np.allclose(p.sig_min.values, pontos.sig_min.values, atol=1e-6), tic
                for metodo in ("diagonal", "gp"):
                    r = RT.refazer(p, P) if metodo == "diagonal" else GP.refazer_gp(p, P, "matern32", A, tc)
                    det = r["p_curv"] < P_DET and r["p_adv"] < P_DET and r["vao_ok"]
                    linhas.append({"tic": tic, "delta_min": delta, "sigma_v_min": sig_v, "metodo": metodo, "em_D11": tic in D11, "nucleo": tic in NUCLEO,
                                   "dPdt": r["dPdt"], "sdPdt": r["sdPdt"], "p_curv": r["p_curv"], "p_adv": r["p_adv"], "p_gof": r["p_gof"], "vao_ok": r["vao_ok"],
                                   "curvatura": r["p_curv"] < P_DET, "adversarial": r["p_adv"] < P_DET, "inadequada": r["p_gof"] < P_DET, "detecta": det,
                                   "barra_sw_mediana_min": float(np.median(p.sig_min[~p.fonte.str.startswith("TESS")]))})
        print(f"TIC {tic} ok", flush=True)
    R = pd.DataFrame(linhas)
    # controle: a combinacao atual reproduz a tabela (diagonal)
    base = R[(R.delta_min == ATUAL[0]) & (R.sigma_v_min == ATUAL[1]) & (R.metodo == "diagonal")].set_index("tic")
    assert all(abs(base.loc[t, "dPdt"] - t26.loc[t, "dPdt"]) < 1e-6 and abs(base.loc[t, "sdPdt"] - t26.loc[t, "s"]) < 1e-6 for t in tics), "a combinacao atual nao reproduz a Tabela 3"
    assert all(bool(base.loc[t, "detecta"]) == bool(t26.loc[t, "curvatura"]) for t in tics), "veredito atual difere da tabela"
    # Delta(dP/dt)/sigma contra a combinacao atual do mesmo metodo
    for metodo in ("diagonal", "gp"):
        b = R[(R.delta_min == ATUAL[0]) & (R.sigma_v_min == ATUAL[1]) & (R.metodo == metodo)].set_index("tic")
        m = R.metodo == metodo
        R.loc[m, "delta_dPdt_sobre_sigma"] = [(r.dPdt - b.loc[r.tic, "dPdt"]) / b.loc[r.tic, "sdPdt"] for r in R[m].itertuples()]
        R.loc[m, "razao_sigma"] = [r.sdPdt / b.loc[r.tic, "sdPdt"] for r in R[m].itertuples()]
    R.to_parquet(BASE / "sensib_vies_26.parquet", index=False)
    out = {"expectativa_commit": EXPECTATIVA, "combinacoes": {}}
    print("\n== vereditos por combinacao (delta, sigma_v) e metodo:")
    for metodo in ("diagonal", "gp"):
        b = R[(R.delta_min == ATUAL[0]) & (R.sigma_v_min == ATUAL[1]) & (R.metodo == metodo)].set_index("tic")
        base_det = set(b.index[b.detecta])
        for delta in DELTAS:
            for sig_v in SIGS:
                g = R[(R.delta_min == delta) & (R.sigma_v_min == sig_v) & (R.metodo == metodo)].set_index("tic")
                det = set(g.index[g.detecta]); d11 = det & D11
                entram = sorted(det - base_det); saem = sorted(base_det - det)
                nuc = {t: bool(g.loc[t, "detecta"]) for t in NUCLEO}
                dd = g.delta_dPdt_sobre_sigma.abs()
                key = f"{metodo} delta {delta:+.2f} sigma_v {sig_v:.2f}"
                out["combinacoes"][key] = {"detecta": sorted(int(t) for t in det), "D11_sobrevivem": len(d11), "entram": [int(t) for t in entram], "saem": [int(t) for t in saem],
                                           "nucleo": {str(k): v for k, v in nuc.items()}, "nucleo_intacto": all(nuc.values()), "p_gof_flags": sorted(int(t) for t in g.index[g.inadequada]),
                                           "dDpdt_sigma_mediana": float(dd.median()), "dDpdt_sigma_max": float(dd.max()), "dDpdt_sigma_max_tic": int(dd.idxmax()),
                                           "n_acima_1sigma": int((dd > 1).sum()), "razao_sigma_mediana": float(g.razao_sigma.median())}
                print(f"  {key}: D11 sobrevivem {len(d11)} de 11 | detecta {len(det)} | entram {entram} | saem {saem} | nucleo {'intacto' if all(nuc.values()) else 'QUEBRA ' + str(nuc)} | "
                      f"p_gof < 0,05: {sorted(g.index[g.inadequada])} | |Delta dP/dt|/sigma mediana {dd.median():.2f}, max {dd.max():.2f} ({int(dd.idxmax())}), > 1 em {int((dd > 1).sum())} | sigma x{g.razao_sigma.median():.2f}")
    # por alvo: a excursao de dP/dt entre delta 0 e -4,28 (sigma_v 0,78) contra a barra
    print(f"\n== por alvo, diagonal, sigma_v {ATUAL[1]:.2f}: dP/dt por delta (colunas) e a excursao 0 -> -4,28 em sigmas")
    d = R[(R.metodo == "diagonal") & (R.sigma_v_min == ATUAL[1])].pivot(index="tic", columns="delta_min", values="dPdt")
    s = R[(R.metodo == "diagonal") & (R.sigma_v_min == ATUAL[1]) & (R.delta_min == ATUAL[0])].set_index("tic").sdPdt
    exc = ((d[0.0] - d[-4.28]).abs() / s).rename("excursao_sigma")
    T = pd.concat([d.round(4), s.round(4).rename("sigma"), exc.round(2)], axis=1)
    T["em_D11"] = [t in D11 for t in T.index]; T["nucleo"] = [t in NUCLEO for t in T.index]
    print(T.to_string())
    out["excursao_0_a_-4.28_diag"] = {"mediana_sigma": float(exc.median()), "max_sigma": float(exc.max()), "max_tic": int(exc.idxmax()), "n_acima_1": int((exc > 1).sum()),
                                      "acima_1": sorted(int(t) for t in exc.index[exc > 1])}
    print(f"  excursao 0 -> -4,28 em sigmas: mediana {exc.median():.2f}, max {exc.max():.2f} ({int(exc.idxmax())}), > 1 sigma em {int((exc > 1).sum())} de 26: {sorted(int(t) for t in exc.index[exc > 1])}")
    (BASE / "sensib_vies_26.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(f"  -> {BASE / 'sensib_vies_26.parquet'}, sensib_vies_26.json")
