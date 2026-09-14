# -*- coding: utf-8 -*-
"""Tabela legivel por maquina das epocas dos 26 alvos (segunda rodada de
revisao, item 8: "tabela de epocas para o CDS"). Uma linha por epoca que entra
nos ajustes da Tabela 3: TIC, fonte (setor TESS ou temporada SuperWASP), ciclo
E contado do primeiro setor TESS, instante do minimo primario em BJD_TDB
(= BTJD + 2 457 000), barra atribuida (min), barra formal (min), e a barra das
metades do setor (TESS) ou a dispersao entre temporadas (SuperWASP) de onde a
barra atribuida sai (Secao 3.3). Le os registros por alvo do estado B' e nao
recalcula nada; a sanidade confere que o numero de epocas por alvo e o da
Tabela 3 (n = nT + nS) e que as 145 epocas (84 TESS + 61 SuperWASP) batem com
qual_barra_26.
Sai: reports/epocas_26.csv (UTF-8, cabecalho comentado com '#').
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
OUT = config.ROOT / "reports" / "epocas_26.csv"
BTJD0 = 2457000.0

if __name__ == "__main__":
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    linhas = []
    for tic in sorted(t26.index):
        j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
        pts = sorted(j["pontos"], key=lambda q: q["t0"])
        assert len(pts) == int(t26.loc[tic, "n"]), (tic, len(pts), t26.loc[tic, "n"])
        for q in pts:
            tess = q["fonte"].startswith("TESS")
            linhas.append({"TIC": int(tic), "source": q["fonte"], "E": int(q["E"]), "BJD_TDB": q["t0"] + BTJD0,
                           "sigma_min": q["sig_min"], "sigma_formal_min": q["sig_formal_min"],
                           "sigma_halves_min": q["sig_metades_min"] if tess else np.nan,
                           "season_dispersion_min": np.nan if tess else float(j["superwasp"]["dispersao_entre_temporadas_min"]),
                           "P_d": float(j["P_escada_d"])})
    d = pd.DataFrame(linhas)
    q26 = pd.read_parquet(BASE / "qual_barra_26.parquet")
    nT, nS = int(d.source.str.startswith("TESS").sum()), int((~d.source.str.startswith("TESS")).sum())
    assert (nT, nS) == (int((q26.fonte == "TESS").sum()), int((q26.fonte == "SuperWASP").sum())), (nT, nS, q26.fonte.value_counts().to_dict())
    cab = ("# Primary-minimum epochs of the 26 targets (this paper, Sect. 3.1-3.3). One row per epoch entering the fits of Table 3.\n"
           "# TIC: TESS Input Catalog number. source: TESS sector (2-min cadence, epoch of the sector; two halves give sigma_halves_min)\n"
           "#   or SuperWASP season (label = first JD of the season). E: cycle number counted from the first TESS sector of the target.\n"
           "# BJD_TDB: time of primary minimum, barycentric TDB. sigma_min: bar used in the fits (min) = max(formal, halves) for TESS,\n"
           "#   hypot(max(formal, season dispersion), 0.78) for SuperWASP (Sect. 3.3-3.4; SuperWASP epochs have the chain bias of -2.14 min removed, i.e. +2.14 min applied to the raw fit).\n"
           "# sigma_formal_min: formal (photon-noise) error; sigma_halves_min: |t_a - t_b|/sqrt 2 of the two sector halves (TESS only);\n"
           "#   season_dispersion_min: std (ddof = 1) of the target's season deviations (SuperWASP only). P_d: period of the ladder (d).\n")
    OUT.write_text(cab + d.to_csv(index=False, float_format="%.6f", lineterminator="\n"), encoding="utf-8")
    print(f"{len(d)} epocas ({nT} TESS + {nS} SuperWASP) de {d.TIC.nunique()} alvos -> {OUT}")
