# -*- coding: utf-8 -*-
"""O que poe o nulo do estimador de residuos abaixo de 1: os PISOS da barra.

Quarta rodada (critico externo, B1): a nota explicava o nulo SuperWASP de
0,77 (Seccao 5.2.3) por "z^2 = 0,5 por construcao" - residuo +-Delta/2 e
barra |Delta|/sqrt 2 do mesmo par de temporadas. A razao corrigida por
alavanca de um par assim e exatamente 1, nao 0,77. O que deflaciona e o
piso: a barra e max(formal, dispersao) e leva o 0,78 min em quadratura, e
quando a dispersao de 2-3 amostras cai abaixo do piso a barra sobe e o
residuo nao. Este script liga e desliga cada piso sozinho no mesmo nulo
(`nulo_barras_reestimadas_26.realizacao` com o gerador reescrito aqui com
chaves), 500 realizacoes, e grava `nulo_pisos_26.parquet`.

Foi rodado primeiro como diagnostico fora do repositorio (estado B, barra
TESS |a - b| / sqrt 2): sem pisos SW 1,08 (mediana) e TESS 5,9 (diverge:
metade com dispersao ~0 faz z ilimitado); so o piso formal SW 0,89; so o
hypot 0,78: 0,86; todos os pisos 0,78 / 0,96. Promovido ao repositorio no
estado D (barra TESS |a - b| / 2, o gerador de `barras_reestimadas`).

EXPECTATIVA (estado D, antes de rodar): sem pisos SW ~1,05-1,10 (a parte
SuperWASP nao muda) e TESS ainda divergente (> 3); so piso formal SW
~0,89; so hypot ~0,86; todos os pisos = o nulo de nulo_barras_reestimadas_26
no estado D (SW ~0,77; TESS o valor que aquele script medir).

RESULTADO (estado D, 2026-09-16, 500 realizacoes, medianas): sem pisos SW
1.15 (esperado 1,05-1,10 - DESVIO: com a barra TESS /2 e sem piso o lado
TESS diverge mais e arrasta o ajuste), TESS 4.8 (diverge, esperado); so piso
formal SW 0.93 (esperado ~0,89); so hypot 0,78: 0.90 (esperado ~0,86); so piso
formal TESS: SW 1.04, TESS 1.01; todos os pisos 0.78 / 0.95
(nulo_barras_reestimadas_26 na mesma rodada, 2 000 realizacoes: 0,77 / 0,98 -
duas simulacoes independentes da mesma grandeza). Os dois desvios SW sem
piso ficam na mesma direcao (acima do previsto) e nao mudam a leitura: os
pisos deflacionam o nulo.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cadeia_oc as co  # noqa: E402
import config  # noqa: E402
import nulo_barras_reestimadas_26 as N  # noqa: E402
from qual_barra_26 import alavancas  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
N_SIM = 500


def barras_com_chaves(d, r, rng, piso_formal_sw=False, piso_vies=False, piso_formal_tess=False):
    """A regra de `N.barras_reestimadas`, com cada piso ligado ou desligado. Com as tres chaves ligadas
    tem de coincidir com o gerador do pipeline - afirmado abaixo em `_confere_gerador`."""
    sig = d["sig"].copy()
    tess, formal = d["tess"], d["formal"]
    isw = np.where(~tess)[0]
    if len(isw) >= 2:
        disp = float(np.std(r[isw], ddof=1))
        b = np.maximum(formal[isw], disp) if piso_formal_sw else np.full(len(isw), disp)
        sig[isw] = np.hypot(b, co.VIES_CADEIA_SIG) if piso_vies else b
    it = np.where(tess)[0]
    a = rng.normal(0, d["sig"][it] * np.sqrt(2))
    bb = rng.normal(0, d["sig"][it] * np.sqrt(2))
    est = np.abs(a - bb) / 2.0                      # a regra do estado D
    sig[it] = np.maximum(formal[it], est) if piso_formal_tess else est
    return sig


def _confere_gerador(des):
    rng1, rng2 = np.random.default_rng(11), np.random.default_rng(11)
    for d in des[:5]:
        r = np.random.default_rng(1).normal(0, d["sig"])
        s1 = N.barras_reestimadas(d, r, rng1)
        s2 = barras_com_chaves(d, r, rng2, True, True, True)
        assert np.allclose(s1, s2), (d["tic"], s1, s2)


def realizacao(des, rng, **kw):
    zT, hT, zS, hS = [], [], [], []
    for d in des:
        r = rng.normal(0, d["sig"])
        sig = np.maximum(barras_com_chaves(d, r, rng, **kw), 1e-3)
        w = 1 / sig
        X = np.vstack([np.ones_like(d["E"]), d["E"], d["E"] ** 2]).T
        beta, *_ = np.linalg.lstsq(X * w[:, None], r * w, rcond=None)
        z = (r - X @ beta) / sig
        h = alavancas(d["E"], sig)
        zT.extend(z[d["tess"]]); hT.extend(h[d["tess"]]); zS.extend(z[~d["tess"]]); hS.extend(h[~d["tess"]])
    zT, hT, zS, hS = map(np.array, (zT, hT, zS, hS))
    return (zT ** 2).sum() / (1 - hT).sum(), (zS ** 2).sum() / (1 - hS).sum()


CASOS = (("sem pisos", {}),
         ("so piso formal SW", dict(piso_formal_sw=True)),
         ("so hypot 0,78 SW", dict(piso_vies=True)),
         ("so piso formal TESS", dict(piso_formal_tess=True)),
         ("todos os pisos (= pipeline)", dict(piso_formal_sw=True, piso_vies=True, piso_formal_tess=True)))

if __name__ == "__main__":
    des = N.desenhos()
    _confere_gerador(des)
    linhas = []
    for rot, kw in CASOS:
        rng = np.random.default_rng(3)
        f = np.array([realizacao(des, rng, **kw) for _ in range(N_SIM)])
        linhas.append({"caso": rot, **kw, "SW_media": f[:, 1].mean(), "SW_mediana": float(np.median(f[:, 1])),
                       "SW_p2.5": float(np.percentile(f[:, 1], 2.5)), "SW_p97.5": float(np.percentile(f[:, 1], 97.5)),
                       "TESS_media": f[:, 0].mean(), "TESS_mediana": float(np.median(f[:, 0])),
                       "TESS_p2.5": float(np.percentile(f[:, 0], 2.5)), "TESS_p97.5": float(np.percentile(f[:, 0], 97.5)), "n_sim": N_SIM})
        print(f"{rot:30s}: SW media {f[:, 1].mean():.2f} mediana {np.median(f[:, 1]):.2f} | "
              f"TESS media {f[:, 0].mean():.2f} mediana {np.median(f[:, 0]):.2f}", flush=True)
    R = pd.DataFrame(linhas)
    R.to_parquet(BASE / "nulo_pisos_26.parquet", index=False)
    print("gravado:", BASE / "nulo_pisos_26.parquet")
