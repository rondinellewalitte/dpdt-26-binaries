# -*- coding: utf-8 -*-
"""TIC 232634196, a unica deteccao que sobrevive a reinflacao das barras do
TESS em d = 1,0 min/ano: o primario e o secundario derivam JUNTOS ou em
OPOSICAO? Juntos = variacao de periodo real (a deriva diferencial mede
outra coisa); oposicao = migracao de manchas / movimento apsidal, e a
curvatura do primario e suspeita pelo mesmo mecanismo.

Referencia: a efemeride LINEAR de 20 anos ajustada aos primarios (a mesma
contra a qual a curvatura foi medida; `cadeia_oc.ajustar` grau 1 sobre os
`pontos` do JSON). O-C do primario por setor = t_p - (T0 + P E); O-C do
secundario = (t_s - P/2) - (T0 + P E_s), com E_s arredondado do proprio t_s
e as epocas/barras de `secundarios_26_setores.parquet` (mesmo estimador,
barra das metades). A grandeza que decide e a MUDANCA entre 2019-20
(s19-21, media ponderada) e 2024 (s75) em cada conjunto: Delta_p e
Delta_s, com barras. Mesmo sinal dentro das barras = juntos; sinais
opostos = oposicao. Nota de desenho: uma efemeride ajustada so ao TESS e
aos dois conjuntos ao mesmo tempo repartiria a deriva diferencial em
+-metade por construcao e diria "oposicao" sempre - por isso a referencia
e a de 20 anos, dos primarios, e a comparacao e entre as mudancas.

EXPECTATIVA - declarada como NAO cega: os numeros gravados ja estao a
vista (residuos lineares do primario +2,0/+3,1/+2,9 min em 2020 e -4,5 em
2024; d = t_s - t_p - P/2 de +5,0/+5,9/+7,5 para -2,5). Deles decorre
Delta_p ~ -7 min e Delta_s ~ -15 min: mesmo sinal, secundario duas vezes o
primario. O script formaliza com barras e decide se o "mesmo sinal" e
significativo; se as barras do secundario (1,3-1,4 min) deixarem Delta_s
compativel com zero ou com sinal oposto, a leitura "juntos" cai.

RODADA NO ESTADO B (2026-09-13): as barras das tres temporadas passam de
11,5 para 19,8 min (sqrt 3). EXPECTATIVA antes de rodar: a media das
temporadas (-98 min) e o modo comum / diferencial primario-secundario NAO
mudam (sao epocas, nao barras); o "espalhamento de 2,8 sigma" das tres
temporadas vira ~1,6 sigma e deixa de merecer a frase 'circular' - o
espalhamento e a barra sao a mesma medida, agora sem o sqrt n. p_curv e
p_adv continuam < 1e-3.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cadeia_oc as C  # noqa: E402
import config  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
TIC = 232634196


def media_ponderada(x, s):
    w = 1 / s ** 2
    return float(np.sum(w * x) / np.sum(w)), float(1 / np.sqrt(np.sum(w)))


if __name__ == "__main__":
    j = json.loads((BASE / f"oc__{TIC}.json").read_text(encoding="utf-8"))
    p = pd.DataFrame(j["pontos"])
    P = float(j["P_escada_d"])
    E, t, sig = p.E.values.astype(float), p.t0.values, p.sig_min.values / 1440.0
    coef, res, chi2, cov = C.ajustar(E, t, sig, 1)          # efemeride linear de 20 anos, primarios
    P_lin, T0 = float(coef[0]), float(coef[1])
    S = pd.read_parquet(BASE / "secundarios_26_setores.parquet")
    s = S[S.tic == TIC].sort_values("setor").reset_index(drop=True)
    tp = p[p.fonte.str.startswith("TESS")].copy()
    tp["setor"] = tp.fonte.str.replace("TESS s", "").astype(int)
    m = tp.merge(s, on="setor")
    Ep = m.E.values.astype(float)
    oc_p = (m.t0_btjd_p.values - (T0 + P_lin * Ep)) * 1440
    # o E do secundario vem do proprio instante: o estimador pode ter
    # escolhido o secundario um ciclo antes ou depois do primario do setor
    Es = np.round((m.t0_btjd_s.values - 0.5 * P_lin - T0) / P_lin)
    oc_s = (m.t0_btjd_s.values - 0.5 * P_lin - (T0 + P_lin * Es)) * 1440
    sp, ss = m.sigma_min_p.values, m.sigma_min_s.values
    print(f"TIC {TIC}: P_lin {P_lin:.7f} d; efemeride linear de 20 anos sobre {len(p)} epocas (chi2 {chi2:.1f}, {len(p) - 2} gl)")
    print("  setor   O-C prim (min)      O-C sec (min)      d = sec - prim")
    for k in range(len(m)):
        print(f"  s{int(m.setor[k]):<4d} {oc_p[k]:+8.2f} +- {sp[k]:.2f}   {oc_s[k]:+8.2f} +- {ss[k]:.2f}   {oc_s[k] - oc_p[k]:+7.2f}")
    sw = p[~p.fonte.str.startswith("TESS")]
    print("  SuperWASP (primario) na mesma efemeride:", ", ".join(f"{r:+.1f} +- {b:.1f}" for r, b in zip(sw.res_linear_min, sw.sig_min)), "min")
    cedo = m.setor < 75
    mp0, ep0 = media_ponderada(oc_p[cedo], sp[cedo]); ms0, es0 = media_ponderada(oc_s[cedo], ss[cedo])
    k75 = int(np.where(m.setor == 75)[0][0])
    dp, sdp = oc_p[k75] - mp0, np.hypot(sp[k75], ep0)
    ds, sds = oc_s[k75] - ms0, np.hypot(ss[k75], es0)
    print(f"\n  mudanca 2019-20 -> 2024: primario {dp:+.2f} +- {sdp:.2f} min ({dp / sdp:+.1f} sigma); secundario {ds:+.2f} +- {sds:.2f} min ({ds / sds:+.1f} sigma)")
    print(f"  diferencial (sec - prim): {ds - dp:+.2f} +- {np.hypot(sdp, sds):.2f} min; razao sec/prim {ds / dp:.2f}")
    juntos = (np.sign(dp) == np.sign(ds)) and (abs(dp / sdp) > 2) and (abs(ds / sds) > 2)
    opostos = (np.sign(dp) != np.sign(ds)) and (abs(dp / sdp) > 2) and (abs(ds / sds) > 2)
    veredito = ("JUNTOS: mesmo sinal, ambos significativos" if juntos else
                "OPOSICAO: sinais opostos, ambos significativos" if opostos else
                "indecidivel: pelo menos um dos dois nao e significativo")
    print(f"  veredito: {veredito}")
    out = {"tic": TIC, "P_lin_d": P_lin, "oc_prim_min": oc_p.tolist(), "oc_sec_min": oc_s.tolist(), "setores": m.setor.tolist(),
           "delta_prim_min": dp, "s_delta_prim_min": sdp, "delta_sec_min": ds, "s_delta_sec_min": sds, "veredito": veredito}
    (BASE / f"caso_{TIC}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"  -> {BASE / f'caso_{TIC}.json'}")
