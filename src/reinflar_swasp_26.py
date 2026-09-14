# -*- coding: utf-8 -*-
"""Infla as barras das epocas SuperWASP pelo fator MEDIDO em
`qual_barra_26.py` e refaz o ajuste dos 26 e o epsilon por alvo. Quarta
rodada, decisao 1 do revisor: e o cenario apoiado por medida, ao contrario
da reinflacao do TESS.

Fator: sigma_SW,novo = F x sigma_SW, com F = sqrt(2,73) = 1,65, a raiz do
fator de variancia Sigma z2 / Sigma(1 - h) das epocas SuperWASP (IC de
chi2 1,79-4,64 em variancia, i.e. 1,34-2,15 em barra; o valor central e o
que se roda, e os dois extremos do IC de chi2 tambem, como sensibilidade).
As epocas TESS nao mudam. Tudo o mais como em `reinflar_tess_26.refazer`:
ajuste linear e quadratico, p_curv, p_adv, p_gof, periodo so do TESS e o
criterio do vao (que nao muda, porque so envolve o TESS), e o epsilon por
alvo com as barras novas.

EXPECTATIVA (escrita e commitada antes de rodar):
  (a) a mediana de chi2_red vai para perto de 0,73 (o esperado com 1-4
      gl) e o chi2_red global de 1,92 para perto de 1,0 - por construcao
      aproximada, ja que o fator foi medido nesses mesmos residuos;
  (b) as dez deteccoes internas ao TESS FICAM (D11 continua com 9-11):
      a curvatura delas nao depende das temporadas SuperWASP, e tirar peso
      do SuperWASP pode ate aumentar a significancia de algumas;
  (c) TIC 232634196, cuja curvatura e carregada pelas tres temporadas
      SuperWASP a -116/-84/-94 min (barras 11,5 -> 19 min), continua com
      p_curv < 0,05 e passa no adversarial: cada temporada fica a ~5 sigma
      da efemeride linear e o deslocamento adversarial de 19 min nao as
      traz de volta; no extremo superior do IC (F = 2,15, barras 25 min)
      a margem encolhe e ele pode cair;
  (d) Sigma(epsilon) muda pouco (as orbitas injetadas sao detectadas pelo
      padrao TESS), talvez caia de 6,7 para 5-6.
  Desvio em qualquer direcao vai para a nota como esta.

RODADA NO ESTADO B (2026-09-13, depois do conserto do sqrt n em oc_lote):
`fator_medido()` le o fator novo de qual_barra_26.parquet: variancia 1,35
(chi2 0,89-2,29) -> barra 1,16 (0,94-1,51). O extremo inferior e < 1 e
DESINFLA: roda-se como esta, como sensibilidade, e diz-se que e isso.
EXPECTATIVA (escrita e commitada antes de rodar): D11 do estado A ja
virou 10 no estado B (359552377 caiu com a barra certa); com F = 1,16 os
10 ficam (a curvatura deles nao depende das temporadas); 232634196 fica
(barras 19,8 -> 23 min, temporadas a -116/-84/-94 min); em F = 1,51,
9-10; nenhuma nova; mediana do chi2_red 0,67 -> ~0,55; Sigma epsilon ~5,8.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402
import epsilon_26 as EPS  # noqa: E402
import reinflar_tess_26 as RT  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
P_DET = 0.05


def fator_medido():
    d = pd.read_parquet(BASE / "qual_barra_26.parquet")
    g = d[d.fonte == "SuperWASP"]
    from scipy import stats
    esperado = float((1 - g.h).sum())
    razao = float((g.z ** 2).sum()) / esperado
    lo, hi = esperado * razao / stats.chi2.ppf(0.975, esperado), esperado * razao / stats.chi2.ppf(0.025, esperado)
    return np.sqrt(razao), np.sqrt(lo), np.sqrt(hi)


def inflar_swasp(pontos, F):
    p = pontos.copy()
    sw = ~p.fonte.str.startswith("TESS")
    p.loc[sw, "sig_min"] = F * p.sig_min[sw]
    return p


if __name__ == "__main__":
    F, Flo, Fhi = fator_medido()
    print(f"fator de barra SuperWASP medido: {F:.2f} (IC chi2 {Flo:.2f}-{Fhi:.2f})")
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    D11 = set(t26[t26.curvatura].index)
    D9 = set(t26[t26.curvatura & ~t26.inadequada].index)
    linhas = []
    for rot, f in (("central", F), ("IC baixo", Flo), ("IC alto", Fhi)):
        rng = np.random.default_rng(EPS.SEMENTE)
        for tic in sorted(t26.index):
            j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
            pontos = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
            P = float(j["P_escada_d"])
            pr = inflar_swasp(pontos, f)
            r = RT.refazer(pr, P)
            e = EPS.epsilon_alvo(pr, P, rng) if rot == "central" else {"epsilon": np.nan, "epsilon3": np.nan}
            det = (r["p_curv"] < P_DET) and (r["p_adv"] < P_DET) and r["vao_ok"]
            linhas.append({"cenario": rot, "fator": f, "tic": int(tic), **r, "detecta": det,
                           "detecta_adequado": det and r["p_gof"] >= P_DET, "em_D11": tic in D11, "em_D9": tic in D9,
                           "epsilon": e["epsilon"], "epsilon3": e["epsilon3"],
                           "dPdt_original": float(t26.loc[tic, "dPdt"]), "sdPdt_original": float(t26.loc[tic, "s"])})
            if rot == "central":
                print(f"F={f:.2f}: TIC {tic} | dP/dt {r['dPdt']:+.4f}+-{r['sdPdt']:.4f} (era {t26.loc[tic, 'dPdt']:+.4f}+-{t26.loc[tic, 's']:.4f}) "
                      f"| p_curv {r['p_curv']:.3f} p_adv {r['p_adv']:.3f} p_gof {r['p_gof']:.3f} chi2r {r['chi2r']:.2f} | eps {e['epsilon']:.2f} | "
                      f"{'D11' if tic in D11 else '   '} {'sobrevive' if det else 'cai' if tic in D11 else ('NOVA' if det else '')}", flush=True)
    R = pd.DataFrame(linhas)
    R.to_parquet(BASE / "reinflar_swasp_26.parquet", index=False)
    print("\n== resumo (expectativa: mediana chi2r ~0,73; D11 fica com 9-11; 232634196 sobrevive no central; Sigma eps 5-6)")
    for rot, g in R.groupby("cenario", sort=False):
        d11 = g[g.em_D11]; d9 = g[g.em_D9]
        alvo = g[g.tic == 232634196].iloc[0]
        print(f"  {rot} (F = {g.fator.iloc[0]:.2f}): chi2_red mediana {g.chi2r.median():.2f} | global {(g.chi2r * (t26.loc[g.tic].dof.values)).sum() / t26.dof.sum():.2f} "
              f"| D11 sobrevivem {int(d11.detecta.sum())} de 11 | D9 com ajuste adequado {int(d9.detecta_adequado.sum())} de 9 | novas fora de D11: {int(g[~g.em_D11].detecta.sum())} "
              f"({', '.join(str(t) for t in g[~g.em_D11 & g.detecta].tic)}) | 232634196: dP/dt {alvo.dPdt:+.4f}+-{alvo.sdPdt:.4f}, p_curv {alvo.p_curv:.4f}, p_adv {alvo.p_adv:.4f} -> "
              f"{'sobrevive' if alvo.detecta else 'cai'}"
              + (f" | Sigma epsilon {g.epsilon.sum():.2f} (eps3 {g.epsilon3.sum():.2f})" if rot == "central" else ""))
    print(f"  -> {BASE / 'reinflar_swasp_26.parquet'}")
