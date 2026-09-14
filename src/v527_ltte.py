# -*- coding: utf-8 -*-
"""V527 Dra (TIC 424461577): o LTTE publicado, amostrado nas NOSSAS seis
epocas, vira curvatura secular na parabola?

A publicacao (2024MNRAS.532.3582C) descreve o ETV como orbita de terceiro
corpo com P3 = 2,733 a e semi-amplitude 6,6 min, SEM termo secular. A nossa
cadeia deu dP/dt = -0,0359 +- 0,0038 s/ano, p_curv ~ 0, p_adv ~ 0 e
chi2_red = 2,18. Se a leitura "a quadratica comeu o ciclo subamostrado"
esta certa, injetar o LTTE publicado numa efemeride LINEAR nos mesmos seis
instantes, com as mesmas barras, tem de reproduzir esse tipo de resultado.

Injecao: orbita circular (a excentricidade nao entra no que se testa),
tau(t) = A sin(2 pi (t - t_p) / P3), A = 6,6 min, P3 = 2,733 a, fase t_p
varrida em 360 passos. A cada fase: parabola ponderada com `cadeia_oc.ajustar`
(a mesma do lote), dP/dt = 2Q/P, p_curv por delta-chi2, chi2_red, e o
adversarial de `cadeia_oc.adversarial`. Nao ha ruido injetado alem do LTTE:
o que se mede e o que a parabola faz com o ciclo, nao com a barra.

EXPECTATIVA (antes de rodar): (a) |dP/dt| varre ate >= 0,036 s/ano ao longo
da fase; (b) p_curv < 0,05 na maioria das fases; (c) chi2_red < 3 na maioria
das fases - porque as seis epocas caem em tres grupos de fase de LTTE
(2008,4-2008,6; 2019,9-2020,1; 2024,1) e tres parametros passam por tres
grupos, sobrando so a estrutura intra-grupo. Se (a) falhar - o LTTE
publicado nao alcanca 0,036 nas nossas datas - a leitura cai e o desacordo
com a literatura fica sem explicacao. Desvio em qualquer direcao vai para a
nota como esta.

RODADA NO ESTADO B (2026-09-13): as duas temporadas SuperWASP do V527 Dra
passam de barra 1,58 para 2,09 min (x1,32). EXPECTATIVA antes de rodar:
(a) e (b) iguais - o LTTE de 6,6 min e detectado pelo padrao TESS; a fracao
de fases que sobrevive ao adversarial cai pouco, de 83% para 75-83%;
chi2_red observado 2,18 -> 2,03 (ja recomputado pela cadeia).
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

TIC = 424461577
A_MIN = 6.6
P3_ANOS = 2.733
S_POR_ANO = 365.25 * 86400.0


def injetar(pontos, P, fase):
    E = pontos.E.values.astype(float)
    sig = pontos.sig_min.values / 1440.0
    t_lin = pontos.t0.values[pontos.fonte.values == "TESS s19"][0] + P * E     # efemeride linear exata
    tau = (A_MIN / 1440.0) * np.sin(2 * np.pi * (t_lin - t_lin[2]) / (P3_ANOS * 365.25) + 2 * np.pi * fase)
    t = t_lin + tau
    _, _, chi2l, _ = C.ajustar(E, t, sig, 1)
    coef, _, chi2p, cov = C.ajustar(E, t, sig, 2)
    Q, sQ = coef[0], np.sqrt(cov[0, 0])
    ag = pd.DataFrame({"E": E, "t0": t, "sig_min": pontos.sig_min.values})
    adv = C.adversarial(ag, P, float(t[2]))
    return {"fase": fase, "dPdt": 2 * Q / P * S_POR_ANO, "sdPdt": 2 * sQ / P * S_POR_ANO,
            "p_curv": float(stats.chi2.sf(chi2l - chi2p, 1)), "chi2_red": chi2p / (len(E) - 3),
            "p_adv": adv["p_adversarial"]}


if __name__ == "__main__":
    j = json.loads((config.DATA / "orquestra" / "oc_lote" / f"oc__{TIC}.json").read_text(encoding="utf-8"))
    pontos = pd.DataFrame(j["pontos"])
    P = float(j["P_escada_d"])
    obs = j["parabola"]
    res = pd.DataFrame([injetar(pontos, P, f) for f in np.arange(0, 1, 1 / 360)])
    out = config.DATA / "orquestra" / "oc_lote" / "v527_ltte_injetado.parquet"
    res.to_parquet(out, index=False)
    print(f"V527 Dra: LTTE publicado (A {A_MIN} min, P3 {P3_ANOS} a) injetado em efemeride linear nas 6 epocas, 360 fases")
    print(f"  observado: dP/dt {obs['dPdt_s_por_ano']:+.4f} +- {obs['sdPdt_s_por_ano']:.4f}  chi2_red {obs['chi2'] / obs['dof']:.2f}  "
          f"p_curv {j['p_curvatura']:.1e}  p_adv {j['adversarial']['p_adversarial']:.1e}")
    print(f"  injetado:  dP/dt de {res.dPdt.min():+.4f} a {res.dPdt.max():+.4f} s/ano (|dP/dt| >= 0,036 em {(res.dPdt.abs() >= 0.036).mean():.0%} das fases)")
    print(f"             chi2_red mediana {res.chi2_red.median():.2f}, < 3 em {(res.chi2_red < 3).mean():.0%} das fases, maximo {res.chi2_red.max():.2f}")
    print(f"             p_curv < 0,05 em {(res.p_curv < 0.05).mean():.0%} das fases; p_adv < 0,05 em {(res.p_adv < 0.05).mean():.0%}")
    perto = res[(res.dPdt - obs["dPdt_s_por_ano"]).abs() < 0.005]
    print(f"             fases com dP/dt a < 0,005 do observado: {len(perto)} de 360; nelas chi2_red {perto.chi2_red.min():.2f}-{perto.chi2_red.max():.2f}"
          if len(perto) else "             nenhuma fase chega a < 0,005 do observado")
    print(f"  -> {out}")
