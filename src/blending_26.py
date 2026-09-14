# -*- coding: utf-8 -*-
"""Segunda rodada de revisao, item 6: blending. Duas medidas por alvo, nos 26:

(a) CROWDSAP e FLFRCSAP do cabecalho SPOC (HDU 1) de CADA setor de 2 min que
    entrou nas epocas da Tabela 3 - os FITS ja estao em cache (oc_lote.CACHE);
    o setor usado e lido dos `pontos` do registro por alvo ("TESS sNN").
    CROWDSAP = fracao do fluxo na abertura que vem do alvo (1 = sem
    contaminacao); FLFRCSAP = fracao do fluxo do alvo capturada pela abertura.
(b) Vizinhos Gaia DR3 a menos de 48" (cerca de dois pixeis do TESS, 21"/px)
    com Delta G < 2 em relacao ao alvo (o alvo = a fonte mais proxima da
    posicao do TIC, a menos de 3"); cone em I/355/gaiadr3 via VizieR (o TAP
    da ESA respondeu 500/503 na primeira tentativa). Sem rede
    ou sem resposta: FALHA alta, nada de tabela parcial.

O que blending pode e nao pode fazer ao TEMPO de um minimo: um vizinho
constante dilui a profundidade e nao move o instante do minimo (o perfil
diluido e simetrico em torno do mesmo t0); so um vizinho VARIAVEL com
estrutura na escala do eclipse pode mover ou fabricar estrutura de timing.
Esta medida diz quem tem vizinho capaz; nao diz se o vizinho varia - isso e
uma limitacao declarada, nao um teste feito.

EXPECTATIVA (escrita e commitada antes de rodar; a do autor):
  - CROWDSAP: mediana sobre os 26 (mediana dos setores de cada um) > 0,90;
    minimo entre os 26 > 0,5 (Tmag 8,8-13,8 em campos de |b| moderada).
    Nenhum alvo com CROWDSAP < 0,5 - se houver, o alvo entra na nota nomeado.
  - Vizinhos com Delta G < 2 a menos de 48": zero na maioria; no maximo 4 dos
    26 com pelo menos um. Se algum membro de D10 tiver vizinho a menos de 48"
    com Delta G < 2, ele e nomeado na nota como nao separavel por esta via.
  - Sanidade: o numero de setores lidos por alvo = nT da Tabela 3 (84 no total);
    a fonte Gaia mais proxima esta a < 3" em 26 de 26 e tem G compativel com
    Tmag (|G - Tmag| < 1,5).
Desvio em qualquer direcao vai para o humano antes de qualquer outra coisa.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402
import oc_lote  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
RAIO_ARCSEC = 48.0
DELTA_G = 2.0
CASA_ARCSEC = 3.0


def crowdsap_por_setor(tic, setores):
    out = []
    for s in setores:
        cands = sorted(oc_lote.CACHE.glob(f"*-s{s:04d}-{tic:016d}-*_lc.fits"))
        if len(cands) != 1:
            raise RuntimeError(f"FALHA [blending]: TIC {tic} setor {s}: {len(cands)} FITS em cache (esperado 1)")
        with fits.open(cands[0]) as h:
            assert int(h[0].header["TICID"]) == tic and int(h[0].header["SECTOR"]) == s
            out.append({"tic": tic, "setor": s, "crowdsap": float(h[1].header["CROWDSAP"]), "flfrcsap": float(h[1].header["FLFRCSAP"])})
    return out


def vizinhos_gaia(tic, ra, dec):
    """Gaia DR3 (I/355/gaiadr3) via VizieR: o arquivo ESA (TAP) respondeu 500/503 (sobrecarga) na primeira tentativa;
    mesmo catalogo, mesmo criterio. Sem resposta: falha alta."""
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier
    import astropy.units as u
    v = Vizier(columns=["Source", "RA_ICRS", "DE_ICRS", "Gmag", "+_r"], row_limit=500)
    r = v.query_region(SkyCoord(ra * u.deg, dec * u.deg), radius=RAIO_ARCSEC * u.arcsec, catalog="I/355/gaiadr3")
    if len(r) == 0 or len(r[0]) == 0:
        raise RuntimeError(f"FALHA [blending]: TIC {tic}: cone Gaia DR3 (VizieR) vazio a {RAIO_ARCSEC}\" (nem o proprio alvo)")
    t = r[0].to_pandas().rename(columns={"_r": "sep_arcsec", "Gmag": "phot_g_mean_mag"})
    t = t.sort_values("sep_arcsec").reset_index(drop=True)
    alvo = t.iloc[0]
    if alvo.sep_arcsec > CASA_ARCSEC:
        raise RuntimeError(f"FALHA [blending]: TIC {tic}: fonte Gaia mais proxima a {alvo.sep_arcsec:.1f}\" (> {CASA_ARCSEC}\")")
    viz = t.iloc[1:]
    viz = viz[np.isfinite(viz.phot_g_mean_mag)]
    capazes = viz[viz.phot_g_mean_mag - alvo.phot_g_mean_mag < DELTA_G]
    return {"tic": tic, "G_alvo": float(alvo.phot_g_mean_mag), "sep_alvo_arcsec": float(alvo.sep_arcsec), "n_viz_48": int(len(viz)),
            "n_viz_dG2": int(len(capazes)), "viz_dG2": [(float(r.sep_arcsec), float(r.phot_g_mean_mag - alvo.phot_g_mean_mag)) for r in capazes.itertuples()],
            "G_viz_mais_brilhante": float(viz.phot_g_mean_mag.min()) if len(viz) else np.nan}


if __name__ == "__main__":
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    alvos = pd.read_parquet(config.DATA / "orquestra" / "alvos_oc_88.parquet").set_index("tic")
    D10 = set(t26.index[t26.curvatura])
    cs, vz = [], []
    for tic in sorted(t26.index):
        j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
        setores = sorted(int(q["fonte"].split("s")[-1]) for q in j["pontos"] if q["fonte"].startswith("TESS"))
        assert len(setores) == int(t26.loc[tic, "nT"]), (tic, setores, t26.loc[tic, "nT"])
        cs.extend(crowdsap_por_setor(int(tic), setores))
        v = vizinhos_gaia(int(tic), float(alvos.loc[tic, "ra"]), float(alvos.loc[tic, "dec"]))
        v["tmag"] = float(alvos.loc[tic, "tmag"]); v["em_D10"] = tic in D10
        vz.append(v)
        c = [x for x in cs if x["tic"] == tic]
        print(f"TIC {tic} {'D10' if tic in D10 else '   '} Tmag {v['tmag']:5.2f} G {v['G_alvo']:5.2f} @ {v['sep_alvo_arcsec']:.2f}\" | "
              f"CROWDSAP {' '.join(f'{x['crowdsap']:.3f}' for x in c)} | viz 48\": {v['n_viz_48']}, dG<2: {v['n_viz_dG2']} {v['viz_dG2'] if v['n_viz_dG2'] else ''}", flush=True)
    C = pd.DataFrame(cs); V = pd.DataFrame(vz)
    assert len(C) == 84, len(C)
    assert (V.sep_alvo_arcsec < CASA_ARCSEC).all() and (abs(V.G_alvo - V.tmag) < 1.5).all(), V[["tic", "sep_alvo_arcsec", "G_alvo", "tmag"]]
    med = C.groupby("tic").crowdsap.median()
    print(f"\n== resumo (expectativa: mediana > 0,90, minimo > 0,5; <= 4 alvos com vizinho dG<2 a 48\")")
    print(f"  CROWDSAP por alvo (mediana dos setores): mediana {med.median():.3f}, min {med.min():.3f} (TIC {med.idxmin()}), max {med.max():.3f}; setores < 0,9: {int((C.crowdsap < 0.9).sum())} de 84, < 0,5: {int((C.crowdsap < 0.5).sum())}")
    print(f"  FLFRCSAP: mediana {C.flfrcsap.median():.3f}, min {C.flfrcsap.min():.3f}")
    com = V[V.n_viz_dG2 > 0]
    print(f"  alvos com vizinho a < 48\" e dG < 2: {len(com)} de 26: {[(int(r.tic), r.viz_dG2, 'D10' if r.em_D10 else '') for r in com.itertuples()]}")
    print(f"  D10 com vizinho capaz: {sorted(int(t) for t in com[com.em_D10].tic)}")
    C.to_parquet(BASE / "blending_crowdsap_26.parquet", index=False)
    V.drop(columns=["viz_dG2"]).assign(viz_dG2=V.viz_dG2.map(json.dumps)).to_parquet(BASE / "blending_gaia_26.parquet", index=False)
    print(f"  -> {BASE / 'blending_crowdsap_26.parquet'}, {BASE / 'blending_gaia_26.parquet'}")
