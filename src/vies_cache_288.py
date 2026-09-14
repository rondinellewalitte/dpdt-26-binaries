# -*- coding: utf-8 -*-
"""Quantifica o vies do cache (revisao 5, item 3): para os 288 alvos que
passaram pelo passo 3 do funil (>= 2 setores 2-min separados por >= 1 ano
NOS 22 MANIFESTOS LOCAIS), pergunta ao MAST quantos setores 2-min SPOC cada
um tem de verdade, e reporta a razao real/cache e como ela varia com a
latitude ecliptica. Nao refaz a amostra: transforma "nao medimos o vies" em
"medimos o vies e ele e X". Se passar de ~1 h, para e reporta o que tiver.

Consulta: astroquery Observations.query_criteria(target_name=TIC,
obs_collection="TESS", dataproduct_type="timeseries"), filtrada em
t_exptime == 120 e provenance_name == "SPOC", setores unicos (a mesma
consulta de recuperar_14.setores_intermediarios). Falha de rede e contada
e repetida (3 tentativas), nao escondida.

EXPECTATIVA (escrita e commitada antes de rodar), a partir dos 8 ja
consultados na TAREFA A (4 na zona de visibilidade continua norte com
40-43 setores reais contra 4 em cache; 4 em latitude menor com 7-13 contra
3-8): razao real/cache mediana entre 2 e 3; correlacao positiva forte com
|latitude ecliptica| (Spearman > 0,5); razao >= 8 em quase todos os alvos
com |beta| > 78 graus (CVZ) e < 3 na maioria com |beta| < 40 graus; entre
15% e 25% dos 288 com razao >= 5. Desvio em qualquer direcao vai para a
nota como esta.
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
SAIDA = config.DATA / "results" / "vies_cache_288.parquet"
LIMITE_S = 3600.0


def setores_reais(tic):
    from astroquery.mast import Observations
    ult = None
    for tentativa in range(3):
        try:
            o = Observations.query_criteria(target_name=str(tic), obs_collection="TESS", dataproduct_type="timeseries").to_pandas()
            o = o[(o.t_exptime == 120) & (o.provenance_name == "SPOC")]
            return sorted(int(x) for x in o.sequence_number.unique())
        except Exception as e:  # noqa: BLE001 - contado e repetido
            ult = e
            time.sleep(3.0 * (tentativa + 1))
    raise RuntimeError(f"MAST falhou 3x para TIC {tic}: {str(ult)[:80]}")


if __name__ == "__main__":
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    d = pd.read_parquet(config.DATA / "results" / "dimensionamento_oc.parquet")
    assert len(d) == 288, len(d)
    c = SkyCoord(d.ra.values * u.deg, d.dec.values * u.deg).barycentrictrueecliptic
    d["beta_ecl"] = c.lat.deg
    if SAIDA.exists():                      # retomavel: nao repete o que ja tem
        feito = pd.read_parquet(SAIDA).set_index("tic")
    else:
        feito = pd.DataFrame(columns=["n_real", "setores_reais", "erro"]).set_index(pd.Index([], name="tic"))
    t0 = time.time()
    linhas = []
    for r in d.itertuples():
        if int(r.tic) in feito.index:
            continue
        if time.time() - t0 > LIMITE_S:
            print(f"limite de {LIMITE_S / 60:.0f} min atingido: parando com {len(feito) + len(linhas)} de 288", flush=True)
            break
        try:
            s = setores_reais(int(r.tic))
            linhas.append({"tic": int(r.tic), "n_real": len(s), "setores_reais": ",".join(map(str, s)), "erro": ""})
        except Exception as e:  # noqa: BLE001
            linhas.append({"tic": int(r.tic), "n_real": np.nan, "setores_reais": "", "erro": str(e)[:100]})
        if len(linhas) % 20 == 0:
            print(f"  {len(feito) + len(linhas)} de 288 ({(time.time() - t0) / 60:.1f} min)", flush=True)
            pd.concat([feito.reset_index(), pd.DataFrame(linhas)]).to_parquet(SAIDA, index=False)
    novo = pd.concat([feito.reset_index(), pd.DataFrame(linhas)]) if linhas else feito.reset_index()
    novo.to_parquet(SAIDA, index=False)
    m = d.merge(novo, on="tic", how="left")
    ok = m[m.n_real.notna()].copy()
    ok["razao"] = ok.n_real / ok.n_setores_2min
    ok["abs_beta"] = ok.beta_ecl.abs()
    print(f"\n== vies do cache: {len(ok)} de 288 consultados ({int(m.erro.fillna('').ne('').sum())} com erro de rede)")
    print(f"  setores em cache: mediana {ok.n_setores_2min.median():.0f} | reais: mediana {ok.n_real.median():.0f} (min {ok.n_real.min():.0f}, max {ok.n_real.max():.0f})")
    q = ok.razao.quantile([.1, .25, .5, .75, .9])
    print(f"  razao real/cache: mediana {q[.5]:.2f}, quartis {q[.25]:.2f}-{q[.75]:.2f}, p10-p90 {q[.1]:.2f}-{q[.9]:.2f}; razao >= 5 em {(ok.razao >= 5).mean():.0%}; == 1 em {(ok.razao == 1).mean():.0%}")
    from scipy import stats
    rho = stats.spearmanr(ok.abs_beta.astype(float).values, ok.razao.astype(float).values)
    print(f"  Spearman razao x |beta ecliptica|: rho {rho.statistic:+.2f} (p {rho.pvalue:.1e})")
    for lo, hi in ((0, 40), (40, 60), (60, 78), (78, 90)):
        g = ok[(ok.abs_beta >= lo) & (ok.abs_beta < hi)]
        if len(g):
            print(f"    |beta| {lo:2d}-{hi:2d}: n={len(g):3d}  cache mediana {g.n_setores_2min.median():.0f}  real mediana {g.n_real.median():.0f}  razao mediana {g.razao.median():.1f}  (razao >= 5: {(g.razao >= 5).mean():.0%})")
    # e os 26 medidos?
    t26 = pd.read_parquet(config.DATA / "orquestra" / "oc_lote" / "tabela_26.parquet")
    g = ok[ok.tic.isin(t26.TIC)]
    print(f"  os 26 medidos: cache mediana {g.n_setores_2min.median():.0f}, real mediana {g.n_real.median():.0f}, razao mediana {g.razao.median():.1f} (n={len(g)})")
    print(f"  -> {SAIDA}")
