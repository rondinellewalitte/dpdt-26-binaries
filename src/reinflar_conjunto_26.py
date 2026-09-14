# -*- coding: utf-8 -*-
"""As rodadas que faltavam (revisao 5, itens 2 e 3): a inflacao do TESS que
os PROPRIOS residuos do TESS autorizam, e as configuracoes conjuntas.

  TESS x1,35            limite superior a 95% do fator de barra do TESS pelo
                        teste de residuos (variancia 0,67-1,81 -> barra
                        0,82-1,35). E a unica inflacao do TESS que os dados
                        obrigam a considerar; 1,0 min/ano e teto extrapolado.
  SW x1,65 + TESS x1,35 a correcao conservadora que os dados sustentam.
  SW x1,65 + TESS +1,0  o teto conservador (SuperWASP corrigido E deriva
                        entre setores no TESS).

Sem rodada conjunta o leitor tem de supor que os efeitos sao independentes e
aditivos, o que num ajuste de 3 parametros sobre 4-7 pontos nao sao. Tudo o
mais como em reinflar_tess_26.refazer / reinflar_swasp_26.inflar_swasp;
epsilon por alvo re-medido nos tres cenarios.

EXPECTATIVA (escrita e commitada antes de rodar):
  TESS x1,35: as deteccoes internas ao TESS tem p_curv de 1e-3 a 1e-30 e
    barras de 0,3-1,5 min; multiplicar por 1,35 tira pouco - espero 10-11
    de 11 em D11 (359552377, que depende do SuperWASP, fica), mediana de
    chi2_red ~1,1-1,2, Sigma eps ~6.
  SW x1,65 + TESS x1,35: 9-10 de 11 (359552377 cai pelo SuperWASP;
    229914020 e 377253090 sao os marginais), mediana ~0,5-0,55 (abaixo do
    0,73: os dois fatores juntos sobrecorrigem um pouco, ja que o TESS nao
    exigia correcao), Sigma eps ~5.
  SW x1,65 + TESS +1,0: 1 de 11 (232634196: as temporadas SuperWASP a -100
    min ficam a ~5 sigma mesmo com barras de 19 min) ou 0 se o adversarial
    o derrubar; mediana < 0,5.
  Desvio em qualquer direcao vai para a nota como esta.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402
import epsilon_26 as EPS  # noqa: E402
import reinflar_swasp_26 as RS  # noqa: E402
import reinflar_tess_26 as RT  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
P_DET = 0.05


def fator_tess_residuos():
    """limite superior a 95% (chi2) do fator de barra do TESS em qual_barra_26."""
    from scipy import stats
    d = pd.read_parquet(BASE / "qual_barra_26.parquet")
    g = d[d.fonte == "TESS"]
    esperado = float((1 - g.h).sum())
    razao = float((g.z ** 2).sum()) / esperado
    hi = esperado * razao / stats.chi2.ppf(0.025, esperado)
    return np.sqrt(hi), np.sqrt(razao)


def aplicar(pontos, sw_mult=1.0, tess_mult=1.0, tess_add=0.0):
    p = pontos.copy()
    tess = p.fonte.str.startswith("TESS")
    p.loc[~tess, "sig_min"] = sw_mult * p.sig_min[~tess]
    p.loc[tess, "sig_min"] = tess_mult * p.sig_min[tess]
    if tess_add:
        L = (p.t0[tess].max() - p.t0[tess].min()) / 365.25
        p.loc[tess, "sig_min"] = np.hypot(p.sig_min[tess], tess_add * L)
    return p


if __name__ == "__main__":
    F_sw, _, _ = RS.fator_medido()
    F_tess_hi, F_tess = fator_tess_residuos()
    print(f"fator SuperWASP medido {F_sw:.2f}; TESS: central {F_tess:.2f}, limite superior 95% {F_tess_hi:.2f}")
    cenarios = [("TESS x%.2f" % F_tess_hi, dict(tess_mult=F_tess_hi)),
                ("SW x%.2f + TESS x%.2f" % (F_sw, F_tess_hi), dict(sw_mult=F_sw, tess_mult=F_tess_hi)),
                ("SW x%.2f + TESS +1.0" % F_sw, dict(sw_mult=F_sw, tess_add=1.0))]
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    D11 = set(t26[t26.curvatura].index)
    D9 = set(t26[t26.curvatura & ~t26.inadequada].index)
    linhas = []
    for rot, kw in cenarios:
        rng = np.random.default_rng(EPS.SEMENTE)
        for tic in sorted(t26.index):
            j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
            pontos = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
            P = float(j["P_escada_d"])
            pr = aplicar(pontos, **kw)
            r = RT.refazer(pr, P)
            e = EPS.epsilon_alvo(pr, P, rng)
            det = (r["p_curv"] < P_DET) and (r["p_adv"] < P_DET) and r["vao_ok"]
            linhas.append({"cenario": rot, "tic": int(tic), **r, "detecta": det, "detecta_adequado": det and r["p_gof"] >= P_DET,
                           "em_D11": tic in D11, "em_D9": tic in D9, "epsilon": e["epsilon"], "epsilon3": e["epsilon3"]})
            if tic in D11:
                print(f"{rot}: TIC {tic} dP/dt {r['dPdt']:+.4f}+-{r['sdPdt']:.4f} p_curv {r['p_curv']:.3f} p_adv {r['p_adv']:.3f} chi2r {r['chi2r']:.2f} -> {'sobrevive' if det else 'cai'}", flush=True)
    R = pd.DataFrame(linhas)
    R.to_parquet(BASE / "reinflar_conjunto_26.parquet", index=False)
    med_esp, sims = RT.mediana_esperada_chi2r(t26.dof.values)
    print("\n== resumo (expectativa: TESS x1,35 -> 10-11; SW+TESS x1,35 -> 9-10, mediana ~0,5; SW+TESS +1,0 -> 0-1)")
    for rot, g in R.groupby("cenario", sort=False):
        d11 = g[g.em_D11]; d9 = g[g.em_D9]; m = float(g.chi2r.median())
        print(f"  {rot}: chi2_red mediana {m:.2f} (P(<= {m:.2f} | barras corretas) = {(sims <= m).mean():.3f}), global {(g.chi2r * t26.loc[g.tic].dof.values).sum() / t26.dof.sum():.2f} "
              f"| D11 sobrevivem {int(d11.detecta.sum())} de 11 ({', '.join(str(t) for t in d11[~d11.detecta].tic)} caem) | D9 adequados {int(d9.detecta_adequado.sum())} de 9 "
              f"| novas: {int(g[~g.em_D11].detecta.sum())} | Sigma eps {g.epsilon.sum():.2f}")
    print(f"  -> {BASE / 'reinflar_conjunto_26.parquet'}")
