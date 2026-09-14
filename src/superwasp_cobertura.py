# -*- coding: utf-8 -*-
"""Cobertura do SuperWASP por cone search, SEM baixar curva.

A busca por coordenada do arquivo publico devolve, para cada fonte, o
identificador 1SWASP, o NUMERO DE PONTOS e a PRIMEIRA e ULTIMA data. Isso
converte a condicao (c) do dimensionamento de predicao em medida para a parte
que importa mais - se existe dado - e da de graca duas coisas que o
dimensionamento precisava:

  - o numero de temporadas (aglomerados de epoca dentro do proprio arquivo),
    que decide "medido" contra "testado";
  - o span do arquivo por alvo, que e a alavanca real e nao a nominal.

O que ela NAO da e o rms fotometrico: para isso a curva tem que ser baixada, e
e o que a calibracao do lado brilhante vai fazer numa amostra.

A pagina traz por fonte: nome, NPTS, inicio e fim da observacao, R.A., Dec,
MAGNITUDE e distancia. A magnitude e a do proprio SuperWASP - e o eixo certo
para a relacao rms(mag), sem passar por Tmag como proxy de V.

A posicao vem da tabela E do identificador (que a codifica), e as duas sao
comparadas: discordancia derruba a linha em vez de virar casamento por nome.
"""
import re
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

import config

BASE = "https://wasp.cerit-sc.cz/search"
CACHE = Path(config.CACHE) / "superwasp" / "cone"
RAIO_PADRAO_DEG = 0.02          # 72 arcsec
TIMEOUT_S = 45
TENTATIVAS = 4
ESPERA_BASE_S = 2.0
ESPERA_MAX_S = 30.0

_ID = re.compile(r"1SWASP\s+J(\d{2})(\d{2})(\d{2}\.\d{2})([+-])(\d{2})(\d{2})(\d{2}\.\d)")


def _coord_do_id(sid):
    """RA, Dec em graus a partir do identificador. Ele codifica a posicao."""
    m = _ID.search(sid)
    if not m:
        raise RuntimeError(f"FALHA [swasp cone]: identificador fora do formato: {sid!r}")
    hh, mm, ss, sinal, dd, dm, ds = m.groups()
    ra = 15.0 * (int(hh) + int(mm) / 60.0 + float(ss) / 3600.0)
    dec = int(dd) + int(dm) / 60.0 + float(ds) / 3600.0
    return ra, (-dec if sinal == "-" else dec)


def _sep_arcsec(ra1, dec1, ra2, dec2):
    d = np.radians
    x = np.cos(d(dec1)) * np.cos(d(dec2)) * np.cos(d(ra1 - ra2)) + np.sin(d(dec1)) * np.sin(d(dec2))
    return float(np.degrees(np.arccos(np.clip(x, -1, 1))) * 3600.0)


def _baixar_html(ra, dec, raio):
    CACHE.mkdir(parents=True, exist_ok=True)
    p = CACHE / f"{ra:.5f}_{dec:+.5f}_{raio:.4f}.html"
    if p.exists() and p.stat().st_size > 500:
        return p.read_text(encoding="utf-8", errors="replace")
    url = f"{BASE}?ra={ra:.6f}&dec={dec:.6f}&radius={raio:.5f}"
    req = urllib.request.Request(url, headers={"User-Agent": "nasa_discover/0.1"})
    # Falha de rede NAO pode virar "sem fonte": a primeira passada nos 288
    # devolveu 150 `getaddrinfo failed` seguidos, e sem repeticao isso teria
    # entrado no numero como ausencia de cobertura.
    ult = None
    for tentativa in range(TENTATIVAS):
        try:
            with urllib.request.urlopen(req, timeout=TIMEOUT_S) as r:
                html = r.read().decode("utf-8", errors="replace")
            break
        except Exception as e:  # noqa: BLE001
            ult = e
            time.sleep(min(ESPERA_BASE_S * (2 ** tentativa), ESPERA_MAX_S))
    else:
        raise RuntimeError(f"FALHA [swasp cone]: {ra},{dec} nao respondeu em "
                           f"{TENTATIVAS} tentativas: {ult}")
    if "SuperWASP" not in html:
        raise RuntimeError(f"FALHA [swasp cone]: resposta inesperada para {ra},{dec} "
                           f"({len(html)} bytes) - sem resposta NAO e ausencia de fonte")
    p.write_text(html, encoding="utf-8")
    return html


def _celulas(linha_html):
    return [re.sub(r"<[^>]*>", " ", c).strip()
            for c in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", linha_html, re.S)]


def consultar(ra, dec, raio=RAIO_PADRAO_DEG):
    """Fontes 1SWASP dentro do raio. Lista de dicts, vazia se nao ha nenhuma.

    Le a TABELA pelo cabecalho, celula a celula. A primeira versao usava uma
    regex que atravessava linhas atras de "id ... numero ... data ... data", e
    ela se desalinhava sempre que o NPTS tinha um digito so: pegava o numero da
    linha seguinte. O sintoma visivel foram 12 spans NEGATIVOS - data final
    antes da inicial, impossivel - e o invisivel seria todo alvo cujo par de
    datas veio da fonte errada. Cabecalho e a fonte da ordem das colunas.

    A pagina da mais do que eu supunha: alem de NPTS e das duas datas, ela traz
    R.A., Declination, Magnitude e "Distance in deg". A separacao NAO precisa
    ser deduzida do identificador - mas continua sendo, como CONFERENCIA: as
    duas concordam ou a linha e descartada.
    """
    html = _baixar_html(ra, dec, raio)
    m = re.search(r"<table[^>]*>(.*?)</table>", html, re.S)
    if not m:
        # A pagina de "nenhuma fonte" e a MESMA pagina, com titulo "found
        # objects" e sem tabela. Isso e ausencia medida; resposta que nem chega
        # a ser a pagina certa ja levantou em `_baixar_html`.
        if "found objects" in html:
            return []
        raise RuntimeError(f"FALHA [swasp cone]: {ra},{dec} respondeu algo que nao e "
                           "a pagina de resultados nem a de zero fontes")
    linhas = re.findall(r"<tr[^>]*>(.*?)</tr>", m.group(1), re.S)
    if not linhas:
        return []
    cab = _celulas(linhas[0])
    try:
        col = {n: cab.index(n) for n in ("Object name", "NPTS", "Observation start",
                                         "Observation stop", "R.A.", "Declination",
                                         "Magnitude")}
    except ValueError as e:
        raise RuntimeError(f"FALHA [swasp cone]: cabecalho mudou ({cab}): {e}")
    out = []
    for l in linhas[1:]:
        c = _celulas(l)
        if len(c) <= max(col.values()):
            continue
        sid = re.sub(r"\s+", " ", c[col["Object name"]]).strip()
        if not sid.startswith("1SWASP"):
            continue
        rs, ds = float(c[col["R.A."]]), float(c[col["Declination"]])
        # a posicao da pagina e a do identificador tem que bater; se nao
        # baterem, a linha esta desalinhada e nao entra
        ri, di = _coord_do_id(sid)
        if _sep_arcsec(rs, ds, ri, di) > 5.0:
            raise RuntimeError(f"FALHA [swasp cone]: {sid} com posicao da tabela "
                               f"({rs},{ds}) discordando do identificador ({ri},{di})")
        out.append({"sourceid": sid, "npts": int(c[col["NPTS"]]),
                    "data_ini": c[col["Observation start"]][:10],
                    "data_fim": c[col["Observation stop"]][:10],
                    "ra_swasp": rs, "dec_swasp": ds,
                    "mag_swasp": float(c[col["Magnitude"]]) if c[col["Magnitude"]] else float("nan"),
                    "sep_arcsec": _sep_arcsec(ra, dec, rs, ds)})
    return sorted(out, key=lambda x: x["sep_arcsec"])


def melhor(ra, dec, raio=RAIO_PADRAO_DEG, sep_max_arcsec=20.0):
    """A fonte mais proxima dentro de `sep_max_arcsec`, ou None."""
    c = consultar(ra, dec, raio)
    if not c:
        return None
    return c[0] if c[0]["sep_arcsec"] <= sep_max_arcsec else None


def anos(data):
    a, m, d = (int(x) for x in str(data).split("-"))
    import datetime as dt
    return a + (dt.date(a, m, d).timetuple().tm_yday - 1) / 365.25


def cobrir(alvos, raio=RAIO_PADRAO_DEG, pausa_s=0.3, verboso=True):
    """Cone search em lote. `alvos` precisa de tic, ra, dec."""
    linhas = []
    for i, r in enumerate(alvos.itertuples(), 1):
        try:
            m = melhor(float(r.ra), float(r.dec), raio)
            erro = None
        except Exception as e:  # noqa: BLE001 - contado, nunca lido como ausencia
            m, erro = None, str(e)[:120]
        d = {"tic": int(r.tic), "erro": erro}
        if m:
            d.update(m)
            d["span_arquivo_anos"] = anos(m["data_fim"]) - anos(m["data_ini"])
        linhas.append(d)
        if verboso and i % 25 == 0:
            print(f"  [swasp cone] {i}/{len(alvos)}", flush=True)
        if erro is None and not (CACHE / f"{float(r.ra):.5f}_{float(r.dec):+.5f}_{raio:.4f}.html").exists():
            time.sleep(pausa_s)
    return pd.DataFrame(linhas)
