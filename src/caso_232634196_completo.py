# -*- coding: utf-8 -*-
"""Rodada seis, Passo F: TIC 232634196 com TODAS as epocas primarias TESS do
cache (43 setores SPOC 2-min, s14-s86, 2019,58-2024,93; nenhum > 2025,0 no
cache - o item 3, teste fora da amostra com setores > 2025,0, NAO pode ser
feito e e dito assim) mais as 3 temporadas SuperWASP do registro, com o GP
agrupado Matern-3/2 do Passo B (deriva_gp_26.json). O desenho das tabelas
nao muda: isto e validacao.

Dados: epocas primarias de jitter_primarios_26.parquet (estimador do estado
D), E e criterios por `jitter_primarios_26.classificar_setores` (P, T0 da
escada do registro); temporadas SuperWASP dos `pontos` (com a correcao
-2,14 e a barra hypot(max(formal, disp), 0,78) atuais). C = diag(sigma^2) +
K nas TESS (`deriva_gp_26.cov_desenho`), GLS (`gls_vies_26.ajustar_gls`).
1. Quadratica so no bloco TESS (43 epocas, 2019,6-2024,9): dP/dt +- sigma,
   Delta chi2 lin - quad, p_gof; contra -0,215 +- 0,021 das 7 epocas.
2. Quadratica em tudo (46): dP/dt +- sigma, chi2/nu, p_gof; linear para
   Delta chi2.
3. Nao ha setores > 2025,0 no cache: nao feito. A efemeride local 2024 (linear
   aos 11 setores 76-86, GP entre eles) e construida para o item 4.
4. Previsao do minimo primario mais proximo de 2027,5 (E = 2206, o mesmo de
   previsao_232634196 e do Passo B) com tres efemerides - linear em tudo,
   quadratica em tudo, local 2024 - cada uma com janela sqrt(var modelo +
   var posterior do GP no instante); tambem as janelas das 7 e das 18 epocas
   do Passo B para comparacao. Saidas: caso_232634196_completo.json / .parquet.

EXPECTATIVA (escrita e commitada ANTES de rodar):
  1. So TESS, 43 epocas com GP (A 0,84 min, tau_c 68 d): a curvatura interna
     ao bloco e ~3,6 min pico a pico se dP/dt = -0,215 for uniforme; sigma
     esperado 0,03-0,08 s/ano; valor -0,30 a -0,10; compativel com -0,215 a
     2 sigma em ~60% da minha probabilidade (o sigma_j quadratico de 1,6 min
     do Passo 1 diz que o bloco nao e uma parabola limpa); p_gof < 0,05.
  2. Tudo (46): dP/dt -0,21 +- 0,012-0,018; chi2/nu 1,5-2,5, p_gof < 0,05
     (as 18 epocas com GP ja davam 0,025); Delta chi2 lin - quad > 80.
  4. E = 2206: quadratica-tudo 2461588,80-0,81 +- 1,3-1,8 min; linear-tudo
     +10 a +16 min depois, +- 1,0-1,5; local 2024 entre as duas ou alem da
     linear (a inclinacao +3,4-4,6 min/ano dos 11 setores extrapolada 2,7
     anos = +9 a +12 min), janela +- 4-8 min (sigma_P de 10 meses x 720
     ciclos + GP). As tres janelas nao se sobrepoem a 1 sigma entre
     quadratica e linear; a local sobrepoe uma delas.
Desvio em qualquer direcao vai ao relatorio.

RESULTADO (2026-09-17, contra a expectativa de f5cdac7): 43 setores, nenhum
> 2025,0 (item 3 nao feito). 1. So TESS: dP/dt = -0,702 +- 0,047 s/ano,
-9,4 sigma de -0,215 (esperado -0,30 a -0,10: DESVIO GRANDE); chi2 86/40,
p_gof 0,000. 2. Tudo (46): -0,293 +- 0,019, chi2/nu 4,11, p_gof ~0 (esperado
-0,21 e chi2/nu 1,5-2,5). O O-C das 43 contra a quadratica-tudo, por grupo:
2019,6-2020,5 -2,1 min; 2021,5 +4,9; 2022-23 +2,9; 2024 -1,9 - uma
EXCURSAO de ~7 min com pico em ~2021,5 sobre uma escala de ~3 anos, que
nenhuma parabola descreve; dentro de cada ano a dispersao (1,4-2,2 min) ja
esta acima das barras. 4. E = 2206: quadratica-tudo 2461588,8064 +- 1,80
min; linear-tudo +20,6 min (+- 1,24; 9,4 sigma); local 2024 +21,7 +- 6,5
(3,2 sigma da quadratica) - a efemeride local coincide com a LINEAR global,
nao com a quadratica: P_local - P_escada = +0,46 +- 0,44 s, quando a
quadratica exige P(2024,5) ~2,3 s abaixo da media (6 sigma). Leitura: o
O-C de 2019-2025 nao e o de um dP/dt secular; a curvatura de 20 anos
(-0,215 a -0,29) e a soma da alavanca de 2008 com uma excursao de anos
dentro do TESS; o minimo de 2027,5 decide entre a linear (+20 min) e a
quadratica.
POS-HOC (sugestao do revisor interno, sem expectativa previa): (i) cubica
so-TESS: chi2 47,4/39 contra 86,4/40 da quadratica (Delta chi2 39 por um
parametro, p ~4e-10): o termo cubico e exigido; curvatura local no centro do
bloco -0,62 +- 0,05 (a quadratica da -0,70 +- 0,05: move 1,7 sigma). A
curvatura NAO e constante ao longo de 2019-2025: o -0,70 e a parabola
absorvendo a excursao, nao uma taxa. (ii) Janelas do item 4 escaladas por
sqrt(chi2/nu) do ajuste correspondente: linear-tudo +- 1,24 -> 3,86 min
(x3,11); quadratica-tudo +- 1,80 -> 3,65 (x2,03); local 2024 +- 6,47 -> 9,46
(x1,46); a separacao linear-quadratica de +20,6 min passa de 9,4 para 3,9
sigma escalada - ainda decisiva.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402
import deriva_gp_26 as GP  # noqa: E402
import gls_vies_26 as G  # noqa: E402
import jitter_primarios_26 as J  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
TIC = 232634196
S_POR_ANO = 365.25 * 86400.0
E_PRED = 2206
ANO_CORTE = 2025.0
LOCAIS = list(range(76, 87))
EXPECTATIVA = "f5cdac7"


def ano(btjd):
    return 2000.0 + (np.asarray(btjd) + 2457000.0 - 2451544.5) / 365.25


def ajustar(pp, grau, nome, A, tc, P):
    E = pp.E.values.astype(float); t = pp.t0.values.astype(float)
    C, it = GP.cov_desenho(pp, nome, A, tc)
    coef, r, chi2, cov = G.ajustar_gls(E, t, C, grau)
    dof = len(E) - grau - 1
    out = {"n": int(len(E)), "grau": grau, "chi2": float(chi2), "dof": dof, "chi2r": float(chi2 / dof), "p_gof": float(stats.chi2.sf(chi2, dof)), "coef": coef.tolist(), "cov": cov.tolist()}
    if grau == 2:
        out["dPdt"] = 2 * coef[0] / P * S_POR_ANO; out["sdPdt"] = 2 * np.sqrt(cov[0, 0]) / P * S_POR_ANO
    # predicao em E_PRED: media do modelo + posterior do GP condicionado aos residuos TESS
    x = np.array([E_PRED ** k for k in range(grau, -1, -1)], float)
    m = float(x @ coef); var_m = float(x @ cov @ x)
    tt = t[it]; Ct = C[np.ix_(it, it)]
    kst = GP.kernel(nome, np.abs(m - tt), A, tc) / 1440.0 ** 2
    mu = float(kst @ np.linalg.solve(Ct, r[it])); var_gp = float(A ** 2 / 1440.0 ** 2 - kst @ np.linalg.solve(Ct, kst))
    out["previsao"] = {"E": E_PRED, "t_bjd": 2457000.0 + m + mu, "t_modelo_bjd": 2457000.0 + m, "mu_gp_min": mu * 1440.0, "sig_modelo_min": np.sqrt(var_m) * 1440.0,
                       "sig_gp_min": np.sqrt(max(var_gp, 0.0)) * 1440.0, "sig_total_min": np.sqrt(var_m + max(var_gp, 0.0)) * 1440.0}
    return out


if __name__ == "__main__":
    j = json.loads((BASE / f"oc__{TIC}.json").read_text(encoding="utf-8"))
    P, sP, T0, E_ult = J.escada_do_registro(j)
    S = pd.read_parquet(BASE / "jitter_primarios_26.parquet")
    d = J.classificar_setores(S[(S.tic == TIC) & (S.minimo == "primario")], P, sP, T0, E_ult, 0.0)
    assert d.ok.all(), d[~d.ok]
    tess = pd.DataFrame({"fonte": ["TESS s%d" % s for s in d.setor], "t0": d.t0_btjd.values, "sig_min": d.sig_min.values, "E": d.E.values, "setor": d.setor.values, "ano": ano(d.t0_btjd.values)}).sort_values("t0").reset_index(drop=True)
    pts = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
    sw = pts[pts.fonte.str.startswith("SuperWASP")][["fonte", "t0", "sig_min", "E"]].copy(); sw["setor"] = -1; sw["ano"] = ano(sw.t0.values)
    # controle: as 4 epocas do desenho reproduzem os pontos
    for r in pts[pts.fonte.str.startswith("TESS")].itertuples():
        m = tess[tess.fonte == r.fonte]; assert len(m) == 1 and abs(float(m.t0.iloc[0]) - r.t0) * 1440 < 0.01 and int(m.E.iloc[0]) == int(r.E), r.fonte
    depois = tess[tess.ano > ANO_CORTE]
    print(f"TIC {TIC}: {len(tess)} setores TESS no cache, s{tess.setor.min()}-s{tess.setor.max()}, {tess.ano.min():.3f}-{tess.ano.max():.3f}; > {ANO_CORTE}: {len(depois)} -> item 3 {'feito' if len(depois) else 'NAO feito (sem setores > 2025,0 no cache)'}")
    print(tess[["setor", "t0", "ano", "sig_min", "E"]].round(4).to_string(index=False))
    hip = json.loads((BASE / "deriva_gp_26.json").read_text(encoding="utf-8"))["hiper_agrupados"]["matern32"]
    A, tc = hip["A_min"], hip["tau_c_d"]
    tudo = pd.concat([sw, tess], ignore_index=True).sort_values("t0").reset_index(drop=True)
    local = tess[tess.setor.isin(LOCAIS)].reset_index(drop=True)
    res = {"tic": TIC, "expectativa_commit": EXPECTATIVA, "n_tess": int(len(tess)), "setores": tess.setor.tolist(), "anos": [round(float(a), 4) for a in tess.ano],
           "n_depois_2025": int(len(depois)), "item3": "nao feito: sem setores > 2025,0 no cache", "gp": {"kernel": "matern32", "A_min": A, "tau_c_d": tc}, "P_escada_d": P}
    # 1. so TESS
    q_t = ajustar(tess, 2, "matern32", A, tc, P); l_t = ajustar(tess, 1, "matern32", A, tc, P)
    ref = j["parabola"]
    z = (q_t["dPdt"] - ref["dPdt_s_por_ano"]) / np.hypot(q_t["sdPdt"], ref["sdPdt_s_por_ano"])
    res["1_so_tess"] = {"quad": q_t, "lin": l_t, "delta_chi2": l_t["chi2"] - q_t["chi2"], "dPdt_7epocas": ref["dPdt_s_por_ano"], "s_7epocas": ref["sdPdt_s_por_ano"], "z_vs_7epocas": float(z)}
    print(f"\n1. so TESS ({len(tess)} epocas, GP A {A:.2f} min tau_c {tc:.0f} d): dP/dt {q_t['dPdt']:+.4f} +- {q_t['sdPdt']:.4f} s/ano (7 epocas: {ref['dPdt_s_por_ano']:+.4f} +- {ref['sdPdt_s_por_ano']:.4f}; diferenca {z:+.1f} sigma); "
          f"chi2 quad {q_t['chi2']:.1f}/{q_t['dof']} p_gof {q_t['p_gof']:.3f}; Delta chi2 lin - quad {l_t['chi2'] - q_t['chi2']:.1f}")
    # 2. tudo
    q_a = ajustar(tudo, 2, "matern32", A, tc, P); l_a = ajustar(tudo, 1, "matern32", A, tc, P)
    res["2_tudo"] = {"quad": q_a, "lin": l_a, "delta_chi2": l_a["chi2"] - q_a["chi2"]}
    print(f"2. tudo ({len(tudo)} epocas): dP/dt {q_a['dPdt']:+.4f} +- {q_a['sdPdt']:.4f} s/ano; chi2 {q_a['chi2']:.1f}/{q_a['dof']} = {q_a['chi2r']:.2f}, p_gof {q_a['p_gof']:.4f}; Delta chi2 lin - quad {l_a['chi2'] - q_a['chi2']:.1f}; p_curv {stats.chi2.sf(l_a['chi2'] - q_a['chi2'], 1):.2e}")
    # 4. previsao com tres efemerides + local 2024
    l_loc = ajustar(local, 1, "matern32", A, tc, P)
    P_loc = l_loc["coef"][0]; sP_loc = np.sqrt(l_loc["cov"][0][0])
    res["efemeride_local_2024"] = {"n": int(len(local)), "setores": local.setor.tolist(), "P_d": float(P_loc), "sP_d": float(sP_loc), "P_menos_escada_s": float((P_loc - P) * 86400.0),
                                   "chi2": l_loc["chi2"], "dof": l_loc["dof"], "p_gof": l_loc["p_gof"]}
    prevB = json.loads((BASE / "deriva_gp_26.json").read_text(encoding="utf-8"))["previsao_232634196"]["matern32"]
    res["4_previsao"] = {"E": E_PRED, "linear_tudo": l_a["previsao"], "quadratica_tudo": q_a["previsao"], "local_2024": l_loc["previsao"],
                         "passoB_7_quad": prevB["7"]["quad"], "passoB_18_quad": prevB["18"]["quad"], "passoB_18_lin": prevB["18"]["lin"]}
    print(f"\n4. previsao E = {E_PRED} (minimo primario mais proximo de 2027,5), BJD_TDB:")
    base = q_a["previsao"]["t_bjd"]
    for rot, pv in (("linear em tudo", l_a["previsao"]), ("quadratica em tudo", q_a["previsao"]), ("local 2024 (11 setores 76-86)", l_loc["previsao"]),
                    ("Passo B, 7 epocas quad", prevB["7"]["quad"]), ("Passo B, 18 epocas quad", prevB["18"]["quad"]), ("Passo B, 18 epocas lin", prevB["18"]["lin"])):
        tb = pv.get("t_bjd", pv.get("t_pred_bjd"))
        print(f"   {rot:32s}: {tb:.4f} +- {pv['sig_total_min']:.2f} min (modelo {pv['sig_modelo_min']:.2f}, GP {pv['sig_gp_min']:.2f}, mu_gp {pv['mu_gp_min']:+.2f}) | contra a quadratica-tudo: {(tb - base) * 1440:+.1f} min")
    # diagnostico: O-C das 43 epocas TESS por grupo de ano contra a quadratica das 7 (desenho), a quadratica-tudo e a linear-tudo
    cq7 = np.array(j["parabola"]["coef"]) if "coef" in j["parabola"] else None
    grupos = pd.cut(tess.ano, [2019, 2021, 2022, 2023.5, 2026], labels=["2019,6-2020,5", "2021,5", "2022,0-2023,0", "2024,0-2024,9"])
    oc_q = (tess.t0.values - np.polyval(q_a["coef"], tess.E.values)) * 1440; oc_l = (tess.t0.values - np.polyval(l_a["coef"], tess.E.values)) * 1440
    res["diagnostico_oc_tess_min"] = {}
    print("   O-C das 43 TESS por grupo (min): contra quadratica-tudo / linear-tudo")
    for g, idx in tess.groupby(grupos, observed=True).groups.items():
        res["diagnostico_oc_tess_min"][str(g)] = {"n": int(len(idx)), "quad_tudo_media": float(np.mean(oc_q[idx])), "lin_tudo_media": float(np.mean(oc_l[idx]))}
        print(f"     {g:14s} n {len(idx):2d}: quad {np.mean(oc_q[idx]):+.2f} (sd {np.std(oc_q[idx]):.2f}) | lin {np.mean(oc_l[idx]):+.2f}")
    sep = (l_a["previsao"]["t_bjd"] - q_a["previsao"]["t_bjd"]) * 1440.0
    res["4_previsao"]["separacao_lin_quad_min"] = sep; res["4_previsao"]["separacao_sobre_sigma"] = sep / np.hypot(l_a["previsao"]["sig_total_min"], q_a["previsao"]["sig_total_min"])
    sep_loc = (l_loc["previsao"]["t_bjd"] - q_a["previsao"]["t_bjd"]) * 1440.0
    res["4_previsao"]["separacao_local_quad_min"] = sep_loc; res["4_previsao"]["separacao_local_quad_sobre_sigma"] = sep_loc / np.hypot(l_loc["previsao"]["sig_total_min"], q_a["previsao"]["sig_total_min"])
    print(f"   separacao linear - quadratica {sep:+.1f} min = {res['4_previsao']['separacao_sobre_sigma']:.1f} sigma; local - quadratica {sep_loc:+.1f} min = {res['4_previsao']['separacao_local_quad_sobre_sigma']:.1f} sigma; "
          f"P local - P escada = {res['efemeride_local_2024']['P_menos_escada_s']:+.3f} s (sigma {sP_loc * 86400:.3f} s)")
    # POS-HOC (sugestao do revisor interno, sem expectativa previa): (i) cubica so-TESS - se o coeficiente quadratico
    # se move mais que o proprio sigma quando o termo cubico entra, o -0,70 e a parabola absorvendo a excursao, nao
    # uma taxa; (ii) janelas do item 4 tambem escaladas por sqrt(chi2/nu) do ajuste correspondente (p_gof ~ 0).
    c_t = ajustar(tess, 3, "matern32", A, tc, P)
    # curvatura local da cubica no CENTRO do bloco TESS (E = 0 da escada fica em 2004, fora dos dados): d2t/dE2 = 6 c3 E + 2 c2
    E_mid = float(np.mean(tess.E.values)); cc = np.asarray(c_t["coef"]); cv = np.asarray(c_t["cov"])
    g = np.array([6 * E_mid, 2.0, 0.0, 0.0])
    dPdt_cub = float(g @ cc) / P * S_POR_ANO; s_cub = float(np.sqrt(g @ cv @ g)) / P * S_POR_ANO
    res["pos_hoc"] = {"nota": "sugestao do revisor interno apos o Passo F; sem expectativa previa",
                      "cubica_so_tess": {"chi2": c_t["chi2"], "dof": c_t["dof"], "delta_chi2_quad_menos_cub": q_t["chi2"] - c_t["chi2"], "E_mid": E_mid, "dPdt_em_E_mid": dPdt_cub, "sdPdt_em_E_mid": s_cub,
                                         "dPdt_quad": q_t["dPdt"], "sdPdt_quad": q_t["sdPdt"], "moveu_sigmas": (dPdt_cub - q_t["dPdt"]) / q_t["sdPdt"]},
                      "janelas_escaladas_min": {rot: {"sig_total_min": pv["sig_total_min"], "fator": float(np.sqrt(max(f["chi2r"], 1.0))), "sig_escalada_min": pv["sig_total_min"] * float(np.sqrt(max(f["chi2r"], 1.0)))}
                                                for rot, pv, f in (("linear_tudo", l_a["previsao"], l_a), ("quadratica_tudo", q_a["previsao"], q_a), ("local_2024", l_loc["previsao"], l_loc))}}
    pj = res["pos_hoc"]["janelas_escaladas_min"]
    print(f"   POS-HOC cubica so-TESS: chi2 {c_t['chi2']:.1f}/{c_t['dof']} (quad {q_t['chi2']:.1f}/{q_t['dof']}; Delta chi2 {q_t['chi2'] - c_t['chi2']:.1f}); curvatura local no centro do bloco TESS: dP/dt {dPdt_cub:+.3f} +- {s_cub:.3f} "
          f"(quad {q_t['dPdt']:+.3f} +- {q_t['sdPdt']:.3f}; moveu {(dPdt_cub - q_t['dPdt']) / q_t['sdPdt']:+.1f} sigma)")
    print("   POS-HOC janelas escaladas por sqrt(chi2/nu): " + "; ".join(f"{k} +- {v['sig_total_min']:.2f} -> {v['sig_escalada_min']:.2f} min (x{v['fator']:.2f})" for k, v in pj.items())
          + f"; separacao lin - quad {sep:+.1f} min = {sep / np.hypot(pj['linear_tudo']['sig_escalada_min'], pj['quadratica_tudo']['sig_escalada_min']):.1f} sigma escalada")
    (BASE / "caso_232634196_completo.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    tudo.to_parquet(BASE / "caso_232634196_completo.parquet", index=False)
    print(f"  -> {BASE / 'caso_232634196_completo.json'}, .parquet")
