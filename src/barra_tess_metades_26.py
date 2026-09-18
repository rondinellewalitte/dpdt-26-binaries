# -*- coding: utf-8 -*-
"""Quarta rodada (critico externo, B4): a barra TESS e max(formal, |t_a - t_b| / sqrt 2),
e |t_a - t_b| / sqrt 2 e o desvio-padrao de UMA METADE (uma amostra de duas
metades com erro sigma_h cada: |a - b| ~ sigma_h sqrt 2). A epoca do ponto e o
ajuste do SETOR INTEIRO, cujo erro e sigma_h / sqrt 2 = |t_a - t_b| / 2. E o
espelho do erro de composicao do SuperWASP (std/sqrt n dado a cada temporada,
corrigido em 5.2.2), na direcao oposta: onde a dispersao entre metades manda
(36 das 84 epocas TESS), a barra esta sqrt 2 GRANDE demais.

Aqui: refaz os ajustes dos 26 a partir dos registros (como reinflar_tess_26),
com a barra TESS = max(formal, metades / sqrt 2) - i.e. |a - b| / 2 - e as
barras SuperWASP como estao; tabela dP/dt, sigma, p_curv, p_adv (regra de
cadeia_oc.adversarial), p_gof; e o fator de residuos TESS e SuperWASP com a
correcao de alavanca (qual_barra_26). NAO altera a cadeia nem as tabelas: e a
medida para decidir se a convencao muda (estado D) - decisao do autor.

EXPECTATIVA (escrita e commitada antes de rodar):
  - barras TESS menores em 36 das 84 epocas (as em que metades > formal),
    por um fator mediano ~0,8 nelas; mediana geral inalterada (formal manda);
  - D10: 10 de 10 sobrevivem (barras so encolhem -> Delta chi2 so cresce;
    o adversarial desloca por 1 sigma menor);
  - deteccoes NOVAS entre os 16: 0 ou 1 (candidatos: 126945917, 207496824,
    229461186, que falham so o adversarial);
  - sigma(dP/dt) dos D10: encolhe <= 10% (mediana ~3%);
  - fator TESS nos 26: 1,02 -> 1,2-1,4; nos 23 adequados: 0,55 -> 0,65-0,8
    (mais perto do nulo condicionado 0,70); SuperWASP inalterado (1,35);
  - mediana do chi2_red: 0,67 -> 0,70-0,75; p_gof < 0,05: 3 -> 3 ou 4.
Se um D10 CAIR com barras menores, ha erro de implementacao - parar.
Desvio em qualquer direcao vai para o humano antes de qualquer outra coisa.

RESULTADO (estado B, 2026-09-16): 36 de 84 barras encolhem (fator mediano
0,71); D10 10/10; UMA deteccao nova (390021728, p_adv 0,013); sigma dos D10
razao mediana 0,993, min 0,81; chi2_red 0,67 -> 0,70; marcados 3 -> 3; fator
TESS 1,02 -> 1,59 (esperado 1,2-1,4 - DESVIO), 23 adequados 0,55 -> 0,68;
SuperWASP inalterado. Decisao do autor: corrigir a cadeia (ESTADO D). Este
script fica como registro da medida e RECUSA registros do estado D.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cadeia_oc as co  # noqa: E402
import config  # noqa: E402
from qual_barra_26 import alavancas  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
S_POR_ANO = 86400.0 * 365.25
P_DET = 0.05

if __name__ == "__main__":
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    q = pd.read_parquet(BASE / "qual_barra_26.parquet")
    adequados = set(t26.index[~t26.inadequada]); D10 = set(t26.index[t26.curvatura])
    linhas, zs = [], []
    n_muda = n_tess = 0; fatores = []
    for tic in sorted(t26.index):
        j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
        if "/ 2" in str(j.get("barra_tess", "")):
            raise RuntimeError(f"TIC {tic}: registro do estado D (barra TESS ja e |a - b| / 2) - este script mede "
                               "o efeito da regra /2 sobre registros do estado B e nao se aplica aqui")
        p = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
        tess = p.fonte.str.startswith("TESS").values
        met = p.sig_metades_min.values.astype(float); formal = p.sig_formal_min.values.astype(float)
        nova = p.sig_min.values.astype(float).copy()
        alt = np.maximum(formal, np.where(np.isfinite(met), met / np.sqrt(2.0), 0.0))
        n_tess += int(tess.sum()); n_muda += int((alt[tess] < nova[tess] - 1e-9).sum())
        fatores.extend((alt[tess] / nova[tess]).tolist())
        nova[tess] = alt[tess]
        E, t, sig = p.E.values.astype(float), p.t0.values.astype(float), nova / 1440.0
        P = float(j["P_escada_d"])
        _, _, chi2l, _ = co.ajustar(E, t, sig, 1); coef, res, chi2p, cov = co.ajustar(E, t, sig, 2)
        dof = len(E) - 3
        p_curv = float(stats.chi2.sf(chi2l - chi2p, 1)); p_gof = float(stats.chi2.sf(chi2p, dof))
        p_adv = co.adversarial(pd.DataFrame({"E": E, "t0": t, "sig_min": nova}), P, float(t[0]))["p_adversarial"]
        det = p_curv < P_DET and p_adv < P_DET
        z = res / sig; h = alavancas(E, sig)
        for k in range(len(E)):
            zs.append({"tic": tic, "fonte": "TESS" if tess[k] else "SuperWASP", "z": z[k], "h": h[k]})
        linhas.append({"tic": int(tic), "em_D10": tic in D10, "adequado": tic in adequados, "dPdt_B": t26.loc[tic, "dPdt"], "s_B": t26.loc[tic, "s"],
                       "dPdt": 2 * coef[0] / P * S_POR_ANO, "s": 2 * np.sqrt(cov[0, 0]) / P * S_POR_ANO, "chi2r": chi2p / dof, "p_gof": p_gof,
                       "p_curv": p_curv, "p_adv": p_adv, "detecta_B": bool(t26.loc[tic, "curvatura"]), "detecta": det})
    R = pd.DataFrame(linhas); Z = pd.DataFrame(zs)
    R.to_parquet(BASE / "barra_tess_metades_26.parquet", index=False)
    f = np.array(fatores)
    print(f"barras TESS que encolhem: {n_muda} de {n_tess}; fator mediano entre elas {np.median(f[f < 1 - 1e-9]):.2f}; mediana geral {np.median(f):.2f}")
    d10 = R[R.em_D10]
    print(f"D10: sobrevivem {int(d10.detecta.sum())} de {len(d10)}; caem: {sorted(d10[~d10.detecta].tic.tolist())}; NOVAS: {sorted(R[~R.em_D10 & R.detecta].tic.tolist())}")
    print(f"sigma(dP/dt) D10: razao s/s_B mediana {(d10.s / d10.s_B).median():.3f}, min {(d10.s / d10.s_B).min():.3f}")
    print(f"mediana do chi2_red: {t26.chi2r.median():.2f} -> {R.chi2r.median():.2f}; p_gof < 0,05: {int(t26.inadequada.sum())} -> {int((R.p_gof < 0.05).sum())} {sorted(R[R.p_gof < 0.05].tic.tolist())}")
    for rot, sel in (("26", R.tic), ("23 adequados", R[R.adequado].tic)):
        zz = Z[Z.tic.isin(sel)]
        qq = q if rot == "26" else q[q.tic.isin(adequados)]
        for fonte in ("TESS", "SuperWASP"):
            g = zz[zz.fonte == fonte]; g0 = qq[qq.fonte == fonte]
            print(f"  fator {fonte:9s} {rot:12s}: estado B {float((g0.z ** 2).sum() / (1 - g0.h).sum()):.2f} -> com |a-b|/2 {float((g.z ** 2).sum() / (1 - g.h).sum()):.2f}")
    for r in R[R.em_D10 | R.detecta].itertuples():
        print(f"  {r.tic} {'D10' if r.em_D10 else 'NOVA'}: {r.dPdt_B:+.4f}+-{r.s_B:.4f} -> {r.dPdt:+.4f}+-{r.s:.4f} | p_curv {r.p_curv:.3g} p_adv {r.p_adv:.3g} p_gof {r.p_gof:.3f}")
    print(f"  -> {BASE / 'barra_tess_metades_26.parquet'}")
