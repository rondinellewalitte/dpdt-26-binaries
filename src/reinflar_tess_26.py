# -*- coding: utf-8 -*-
"""Reinfla as barras das epocas TESS pela deriva entre setores e refaz o
ajuste dos 26 e o epsilon por alvo. Terceira rodada da revisao, item 1.

A hipotese, que a nota tratava como dois itens separados: a barra do TESS
e a dispersao meia-a-meia DENTRO do setor (~13 d) e nao ve deriva ENTRE
setores; os secundarios derivam a 1-2,6 min/ano em cinco alvos. Se as
barras do TESS sao pequenas demais, (i) epsilon esta superestimado
(1-3 orbitas vira teto) e (ii) o numero de deteccoes por ruido esta
subestimado. As duas empurram na mesma direcao.

Cenarios (declarados, nao escolhidos por resultado; nao ha valor central):
  d = 1,0 e d = 2,6 min/ano, os extremos da deriva medida nos secundarios.
  Para cada epoca TESS do alvo: sigma_novo = hypot(sigma, d x L), com L a
  alavanca TESS do proprio alvo em anos (primeiro ao ultimo setor). As
  epocas SuperWASP nao mudam. E uma extrapolacao (deriva de secundario
  medida em 5 alvos aplicada as epocas primarias dos 26) e e dita assim.

O que se refaz a partir dos `pontos` gravados de cada alvo (epocas, E e
barras - nada e re-extraido): (a) o ajuste linear so nas epocas TESS para
sP e o criterio do vao sP x N < 0,25 ciclo (se falhar, a contagem de
ciclos ficaria ambigua e o alvo sai); (b) o ajuste conjunto linear e
quadratico com `cadeia_oc.ajustar`, dP/dt +- sigma, p_curv, p_adv
(`cadeia_oc.adversarial`), p_gof; (c) o epsilon por alvo de
`epsilon_26.epsilon_alvo` com as barras reinfladas, mesma populacao e
mesma semente. Conta quantos de D11 e de D9 sobrevivem em cada cenario.

EXPECTATIVA (escrita e commitada antes de rodar):
  d = 1,0 (barras TESS ~4 min): a curvatura interna ao TESS desaparece e
    sobra a alavanca SuperWASP-TESS; sobrevivem os alvos cuja excursao
    quadratica em 20 anos e >= ~30 min - espero 6 a 9 dos 11 de D11;
    Sigma(epsilon) cai para ~2-4.
  d = 2,6 (barras TESS ~11 min): a incerteza do periodo extrapolada ao
    SuperWASP chega a ~30 min; sobrevivem so as excursoes > ~50 min
    (232634196, 359552377, 392536812, 424461577 e talvez 2-3 mais) -
    espero 3 a 6 dos 11; Sigma(epsilon) < 2.
  Em nenhum cenario espero que o vao fique ambiguo (sP x N < 0,25).
  Se a maioria dos 11 cair ja em d = 1,0, as deteccoes sao governadas
  pelas barras e isso vira a manchete da nota. Desvio em qualquer direcao
  vai para a nota como esta.

RODADA NO ESTADO B (2026-09-13, barras SuperWASP corrigidas do sqrt n; D11
virou 10): EXPECTATIVA antes de rodar: d = 1,0 -> 1 dos 10 sobrevive
(232634196: alavanca de 20 anos, temporadas a -116/-84/-94 min com barras
de 19,8 min); d = 2,6 -> 0. As barras SuperWASP maiores tiram peso das
temporadas, mas a excursao de 232634196 e de ~100 min contra ~20 de barra.
`--resumo`: mediana do chi2_red esperada com barras corretas continua
0,73 (so depende dos gl); a observada passa a 0,67 (estado B).
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cadeia_oc as C  # noqa: E402
import config  # noqa: E402
import epsilon_26 as EPS  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
CENARIOS = (1.0, 2.6)     # min/ano
S_POR_ANO = 365.25 * 86400.0
P_DET = 0.05


def reinflar(pontos, d):
    p = pontos.copy()
    tess = p.fonte.str.startswith("TESS")
    L = (p.t0[tess].max() - p.t0[tess].min()) / 365.25
    p.loc[tess, "sig_min"] = np.hypot(p.sig_min[tess], d * L)
    return p, float(L)


def refazer(pontos, P):
    E = pontos.E.values.astype(float)
    t = pontos.t0.values
    sig = pontos.sig_min.values / 1440.0
    tess = pontos.fonte.str.startswith("TESS").values
    # (a) sP so do TESS e o criterio do vao
    _, _, _, cov_t = C.ajustar(E[tess], t[tess], sig[tess], 1)
    sP = float(np.sqrt(cov_t[0, 0]))
    N_gap = float(np.max(np.abs(E[~tess] - E[tess].mean())))
    vao_ok = sP * N_gap < 0.25
    # (b) ajuste conjunto
    _, _, chi2l, _ = C.ajustar(E, t, sig, 1)
    coef, _, chi2p, cov = C.ajustar(E, t, sig, 2)
    dof = len(E) - 3
    p_curv = float(stats.chi2.sf(chi2l - chi2p, 1))
    ag = pd.DataFrame({"E": E, "t0": t, "sig_min": pontos.sig_min.values})
    p_adv = float(C.adversarial(ag, P, float(t[tess][-1]))["p_adversarial"])
    return {"dPdt": 2 * coef[0] / P * S_POR_ANO, "sdPdt": 2 * np.sqrt(cov[0, 0]) / P * S_POR_ANO,
            "chi2r": chi2p / dof, "p_gof": float(stats.chi2.sf(chi2p, dof)), "p_curv": p_curv, "p_adv": p_adv,
            "sP_vezes_N": sP * N_gap, "vao_ok": bool(vao_ok)}


def mediana_esperada_chi2r(dofs, n_sim=20000, semente=1):
    """Mediana de chi2/gl que barras CORRETAS dariam nesta amostra de gl -
    com 1-4 gl ela fica bem abaixo de 1 (mediana de chi2_1 e 0,45). E o
    ponto de comparacao que a nota nao tinha: 1,41 contra 0,73, nao contra 1."""
    rng = np.random.default_rng(semente)
    sims = np.array([np.median(stats.chi2.rvs(dofs, random_state=rng) / dofs) for _ in range(n_sim)])
    return float(np.median(sims)), sims


def resumo(R, t26):
    dofs = t26.dof.values
    med, sims = mediana_esperada_chi2r(dofs)
    lo, hi = np.percentile(sims, [2.5, 97.5])
    obs = float(t26.chi2r.median())
    print(f"\n== calibracao das barras pelo chi2_red: mediana esperada com barras corretas (gl {np.bincount(dofs)[1:].tolist()}) = {med:.2f} "
          f"[{lo:.2f}, {hi:.2f}]; observada com as barras originais {obs:.2f} (P(>= {obs:.2f}) = {(sims >= obs).mean():.4f}; "
          f"fator de barra implicado sqrt({obs:.2f}/{med:.2f}) = {np.sqrt(obs / med):.2f})")
    for d, g in R.groupby("cenario_min_por_ano"):
        m = float(g.chi2r.median())
        print(f"   d = {d}: mediana {m:.2f} (P(<= {m:.2f}) = {(sims <= m).mean():.3f})")


if __name__ == "__main__":
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    D11 = set(t26[t26.curvatura].index)
    D9 = set(t26[t26.curvatura & ~t26.inadequada].index)
    if "--resumo" in sys.argv:
        resumo(pd.read_parquet(BASE / "reinflar_tess_26.parquet"), t26)
        sys.exit()
    linhas = []
    for d in CENARIOS:
        rng = np.random.default_rng(EPS.SEMENTE)
        for tic in sorted(t26.index):
            j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
            pontos = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
            P = float(j["P_escada_d"])
            pr, L = reinflar(pontos, d)
            r = refazer(pr, P)
            e = EPS.epsilon_alvo(pr, P, rng)
            det = (r["p_curv"] < P_DET) and (r["p_adv"] < P_DET) and r["vao_ok"]
            linhas.append({"cenario_min_por_ano": d, "tic": int(tic), "alavanca_tess_anos": L,
                           "barra_tess_mediana_min": float(np.median(pr.sig_min[pr.fonte.str.startswith("TESS")])),
                           **r, "detecta": det, "detecta_adequado": det and r["p_gof"] >= P_DET,
                           "em_D11": tic in D11, "em_D9": tic in D9, "epsilon": e["epsilon"], "epsilon3": e["epsilon3"],
                           "dPdt_original": float(t26.loc[tic, "dPdt"]), "sdPdt_original": float(t26.loc[tic, "s"])})
            print(f"d={d}: TIC {tic} barra TESS {linhas[-1]['barra_tess_mediana_min']:.1f} min | dP/dt {r['dPdt']:+.4f}+-{r['sdPdt']:.4f} "
                  f"(era {t26.loc[tic, 'dPdt']:+.4f}+-{t26.loc[tic, 's']:.4f}) | p_curv {r['p_curv']:.3f} p_adv {r['p_adv']:.3f} p_gof {r['p_gof']:.3f} "
                  f"| vao sPxN {r['sP_vezes_N']:.3f} {'ok' if r['vao_ok'] else 'AMBIGUO'} | eps {e['epsilon']:.2f} | "
                  f"{'D11' if tic in D11 else '   '} {'sobrevive' if det else 'cai' if tic in D11 else ''}", flush=True)
    R = pd.DataFrame(linhas)
    R.to_parquet(BASE / "reinflar_tess_26.parquet", index=False)
    print("\n== resumo (expectativa: d=1,0 -> 6-9 de 11 sobrevivem, Sigma eps 2-4; d=2,6 -> 3-6 de 11, Sigma eps < 2; vao nunca ambiguo)")
    for d, g in R.groupby("cenario_min_por_ano"):
        d11 = g[g.em_D11]; d9 = g[g.em_D9]
        print(f"  d = {d} min/ano: barra TESS mediana {g.barra_tess_mediana_min.median():.1f} min | D11 sobrevivem {int(d11.detecta.sum())} de 11 "
              f"({', '.join(str(t) for t in d11[d11.detecta].tic)}) | D9 sobrevivem com ajuste adequado {int(d9.detecta_adequado.sum())} de 9 | "
              f"novas deteccoes fora de D11: {int(g[~g.em_D11].detecta.sum())} | vao ambiguo: {int((~g.vao_ok).sum())} | "
              f"Sigma epsilon {g.epsilon.sum():.2f} (eps3 {g.epsilon3.sum():.2f}) | chi2_red mediana {g.chi2r.median():.2f}")
    resumo(R, t26)
    print(f"  -> {BASE / 'reinflar_tess_26.parquet'}")
