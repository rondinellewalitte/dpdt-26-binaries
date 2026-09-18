# -*- coding: utf-8 -*-
"""Rodada cinco, Passo A: o sistema de tempo do SuperWASP, conferido de forma
independente de `superwasp.hjd_utc_para_bjd_tdb`.

O QUE O ARQUIVO DECLARA (registrado aqui porque a pagina nao diz mais nada):
  - CSV do CERIT-SC (https://wasp.cerit-sc.cz/csv?object=...): cabecalho
    `HJD,camera,magnitude,"magnitude error"`. A pagina inicial nao declara
    escala de tempo nem referencial; so o reconhecimento "DR1 of the WASP
    data (Butters et al. 2010)".
  - Butters et al. 2010 (A&A 520, L10; arXiv:1009.5306), Sec. 3: "TMID is the
    heliocentrically corrected mid-point of the exposure in seconds after
    2004-01-01T00:00:00, and can be converted to HJD using
    HJD = TMID/86400 + 2453005.5 (1)". Nenhuma escala (UTC/TT) e nomeada; a
    origem 2004-01-01T00:00:00 = JD 2453005.5 e uma data civil (UTC), e a
    convencao dos levantamentos da epoca era HJD em UTC.
  - `superwasp.py` assume HJD_UTC heliocentrico, desfaz a correcao
    heliocentrica (astropy, `light_travel_time` heliocentric avaliada no
    proprio HJD), converte UTC -> TDB e aplica a correcao baricentrica.
    Localizacao: SAAO para todos os alvos (os 21 de Draco sao do
    SuperWASP-Norte, La Palma; a diferenca e topocentrica, <= 21 ms).

O QUE ESTE SCRIPT FAZ: 20 timestamps brutos (3 por calibrador HJ - primeiro,
mediano, ultimo - e 3/3/2 de tres alvos dos 26), lidos como texto direto do
CSV; conversao INDEPENDENTE com astropy (Time + light_travel_time, helio e
bary, iterada uma vez para avaliar a correcao heliocentrica no instante
geocentrico, sitio correto por hemisferio: La Palma para dec > 0, SAAO para
dec < 0); comparacao linha a linha com o valor de `superwasp.py`: diferenca
(codigo - independente) em segundos, sinal, e a decomposicao em TDB - UTC,
(bary - helio), efeito do sitio e efeito da iteracao. Tambem as diferencas
que DUAS hipoteses alternativas produziriam contra o codigo: H1, o arquivo ja
em TT/TDB (o codigo somaria 65 s a mais: +65,2 s); H2, sinal trocado em
TDB - UTC no codigo (-130,4 s, a hipotese do autor). Grava
`data/orquestra/oc_lote/checa_tempo_swasp.json` e `.parquet`.

EXPECTATIVA (escrita e commitada ANTES de rodar):
  - Hipotese do autor: se houver erro de sinal em TDB - UTC, a diferenca e
    -130,4 +- 1 s em 2006-08.
  - Minha leitura do codigo antes de rodar: `t2.tdb.jd` SOMA TDB - UTC
    (+65,18 s em 2006-08, +64,18 s em 2004-05 - salto de segundo de
    2005-12-31), logo NAO ha erro de sinal: codigo - independente = 0,00 +-
    0,10 s nas 20 linhas (a unica diferenca esperada e a iteracao da
    correcao heliocentrica, |d| < 0,1 s, e o sitio, |d| < 0,03 s).
    bary - helio dentro de +-6 s. A hipotese -130,4 s NAO se confirma.
  - Se confirmar (|codigo - independente| > 60 s em qualquer linha): PARAR e
    reportar com patch proposto, sem aplicar.
Desvio em qualquer direcao vai ao relatorio.

RESULTADO (2026-09-17, contra a expectativa de 4a0f79c): codigo - independente
= +0,002 s em media, |max| 0,019 s (= a iteracao da correcao heliocentrica);
TDB - UTC 64,185 s (2004-05) e 65,18 s (2006-08); bary - helio -0,91..+1,09
s. A hipotese do sinal trocado (-130,4 s) NAO se confirma. O efeito do sitio
e 0 por construcao: o termo topocentrico entra igual em desfazer a helio e
em aplicar a bary e cancela; o do arquivo (La Palma/SAAO verdadeiros) fica
dentro do HJD e esta certo. H1 (arquivo ja em TT) nao e testavel pelo
arquivo e daria O-C de +1,1 min, sinal oposto ao -2,14.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402
import superwasp as sw  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
LA_PALMA = dict(lat=28.7606, lon=-17.8792, alt=2396.0)     # SuperWASP-Norte, ORM
SAAO = dict(lat=sw.SAAO_LAT, lon=sw.SAAO_LON, alt=sw.SAAO_ALT)
HJ = ("WASP-18 b", "WASP-4 b", "WASP-19 b", "WASP-6 b")
ALVOS = ((232634196, 3), (424461577, 3), (229914020, 2))
DECLARADO = {"coluna_csv": 'HJD,camera,magnitude,"magnitude error"', "escala_declarada_cerit": "nenhuma (pagina so cita Butters et al. 2010, DR1)",
             "butters2010_sec3": "TMID is the heliocentrically corrected mid-point of the exposure in seconds after 2004-01-01T00:00:00, "
                                 "and can be converted to HJD using HJD = TMID/86400 + 2453005.5 (1)", "escala_nomeada_em_butters": "nenhuma"}


def sitio(dec):
    s = LA_PALMA if dec > 0 else SAAO
    return EarthLocation(lat=s["lat"] * u.deg, lon=s["lon"] * u.deg, height=s["alt"] * u.m), ("La Palma" if dec > 0 else "SAAO")


def independente(hjd, ra, dec, loc, iterar=True):
    """HJD_UTC (heliocentrico) -> BJD_TDB, escrito do zero. Devolve pecas em segundos."""
    c = SkyCoord(ra * u.deg, dec * u.deg)
    t = Time(hjd, format="jd", scale="utc", location=loc)
    ltt_h = t.light_travel_time(c, kind="heliocentric").to_value(u.s)
    jd_geo = hjd - ltt_h / 86400.0
    if iterar:                                   # a correcao heliocentrica avaliada no instante geocentrico
        tg = Time(jd_geo, format="jd", scale="utc", location=loc)
        ltt_h = tg.light_travel_time(c, kind="heliocentric").to_value(u.s)
        jd_geo = hjd - ltt_h / 86400.0
    tg = Time(jd_geo, format="jd", scale="utc", location=loc)
    tdb_menos_utc = (tg.tdb.jd - tg.utc.jd) * 86400.0
    ltt_b = tg.light_travel_time(c, kind="barycentric").to_value(u.s)
    bjd_tdb = tg.tdb.jd + ltt_b / 86400.0
    return {"bjd_tdb": float(bjd_tdb), "tdb_menos_utc_s": float(tdb_menos_utc), "ltt_helio_s": float(ltt_h), "ltt_bary_s": float(ltt_b),
            "bary_menos_helio_s": float(ltt_b - ltt_h), "bjd_menos_hjd_s": float((bjd_tdb - hjd) * 86400.0)}


def linhas_brutas(path, n):
    """Le o CSV como TEXTO e devolve (string, indice) do primeiro, mediano(s) e ultimo timestamp."""
    txt = [l for l in path.read_text(encoding="utf-8").splitlines()[1:] if l.strip()]
    idx = [0, len(txt) // 2, len(txt) - 1] if n == 3 else [0, len(txt) - 1]
    return [(txt[i].split(",")[0], i) for i in idx]


if __name__ == "__main__":
    cal = pd.read_parquet(config.DATA / "cache" / "superwasp" / "hj_calibradores.parquet").set_index("pl_name")
    alv = pd.read_parquet(config.DATA / "orquestra" / "alvos_oc_88.parquet").set_index("tic")
    curvas = [(n, cal.loc[n, "swasp_id"], float(cal.loc[n, "ra"]), float(cal.loc[n, "dec"]), 3) for n in HJ]
    curvas += [(f"TIC {t}", alv.loc[t, "sourceid"], float(alv.loc[t, "ra"]), float(alv.loc[t, "dec"]), n) for t, n in ALVOS]
    linhas = []
    for nome, sid, ra, dec, n in curvas:
        p = sw.baixar_curva(sid)
        loc, nome_sitio = sitio(dec)
        loc_saao = EarthLocation(lat=SAAO["lat"] * u.deg, lon=SAAO["lon"] * u.deg, height=SAAO["alt"] * u.m)
        for s, i in linhas_brutas(p, n):
            hjd = float(s)
            cod = float(sw.hjd_utc_para_bjd_tdb(np.array([hjd]), ra, dec)[0])
            ind = independente(hjd, ra, dec, loc)
            ind_saao = independente(hjd, ra, dec, loc_saao)
            ind_sem_iter = independente(hjd, ra, dec, loc, iterar=False)
            d = (cod - ind["bjd_tdb"]) * 86400.0
            linhas.append({"curva": nome, "sourceid": sid, "linha_csv": i, "hjd_bruto": s, "ano": 2000.0 + (hjd - 2451544.5) / 365.25, "sitio": nome_sitio,
                           "bjd_tdb_codigo": cod, "bjd_tdb_independente": ind["bjd_tdb"], "codigo_menos_independente_s": d,
                           "sinal": "+" if d > 0 else ("-" if d < 0 else "0"),
                           "tdb_menos_utc_s": ind["tdb_menos_utc_s"], "bary_menos_helio_s": ind["bary_menos_helio_s"], "bjd_menos_hjd_independente_s": ind["bjd_menos_hjd_s"],
                           "efeito_sitio_s": (ind["bjd_tdb"] - ind_saao["bjd_tdb"]) * 86400.0, "efeito_iteracao_s": (ind["bjd_tdb"] - ind_sem_iter["bjd_tdb"]) * 86400.0,
                           "H1_arquivo_em_TT_codigo_menos_correto_s": ind["tdb_menos_utc_s"], "H2_sinal_trocado_codigo_menos_correto_s": -2 * ind["tdb_menos_utc_s"]})
    R = pd.DataFrame(linhas)
    assert len(R) == 20, len(R)
    pd.set_option("display.width", 250)
    print("== 20 timestamps brutos: codigo (superwasp.py) contra conversao independente (astropy), em segundos")
    print(R[["curva", "linha_csv", "hjd_bruto", "ano", "sitio", "codigo_menos_independente_s", "sinal", "tdb_menos_utc_s", "bary_menos_helio_s",
             "bjd_menos_hjd_independente_s", "efeito_sitio_s", "efeito_iteracao_s"]].round(4).to_string(index=False))
    dmax = float(R.codigo_menos_independente_s.abs().max())
    print(f"\n   |codigo - independente| max {dmax:.4f} s (esperado < 0,10); TDB - UTC {R.tdb_menos_utc_s.min():.3f}..{R.tdb_menos_utc_s.max():.3f} s; "
          f"bary - helio {R.bary_menos_helio_s.min():+.2f}..{R.bary_menos_helio_s.max():+.2f} s; sitio max {R.efeito_sitio_s.abs().max():.4f} s; iteracao max {R.efeito_iteracao_s.abs().max():.4f} s")
    print(f"   H2 (sinal trocado em TDB - UTC) daria codigo - correto = {R.H2_sinal_trocado_codigo_menos_correto_s.mean():+.1f} s; observado {R.codigo_menos_independente_s.mean():+.3f} s "
          f"-> hipotese {'CONFIRMADA' if dmax > 60 else 'NAO confirmada'}")
    print(f"   H1 (arquivo ja em TT/TDB) daria codigo - correto = {R.H1_arquivo_em_TT_codigo_menos_correto_s.mean():+.1f} s: nao testavel pelo arquivo (escala nao declarada); "
          f"produziria O-C de +1,1 min nos HJ, o oposto do -2,14 medido")
    out = {"declarado": DECLARADO, "n": int(len(R)), "max_abs_codigo_menos_independente_s": dmax, "media_s": float(R.codigo_menos_independente_s.mean()),
           "tdb_menos_utc_s": [float(R.tdb_menos_utc_s.min()), float(R.tdb_menos_utc_s.max())], "hipotese_sinal_confirmada": bool(dmax > 60),
           "H2_previsto_s": float(R.H2_sinal_trocado_codigo_menos_correto_s.mean()), "H1_previsto_s": float(R.H1_arquivo_em_TT_codigo_menos_correto_s.mean()),
           "linhas": R.to_dict("records")}
    (BASE / "checa_tempo_swasp.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    R.to_parquet(BASE / "checa_tempo_swasp.parquet", index=False)
    print(f"   -> {BASE / 'checa_tempo_swasp.json'}, .parquet")
