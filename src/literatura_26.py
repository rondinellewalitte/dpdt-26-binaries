# -*- coding: utf-8 -*-
"""Cruzamento dos 26 alvos medidos de dP/dt com o que ja existe publicado.

Ponto 2 da revisao da nota: "a nota nao diz se algum dos 26 ja tem O-C
publicado. Se voce nao checar, o referee checa."

Duas fontes, as duas por posicao (o TIC nao e identificador em nenhuma):

  1. Bibliografia do objeto no SIMBAD (TAP: basic -> has_ref -> ref). E o
     mesmo indice que a busca por objeto do ADS consulta; sem chave do ADS,
     e o caminho equivalente. Sai a contagem de referencias e a lista de
     titulos, com os titulos que citam periodo/O-C/minimos marcados.
  2. Catalogo de estrelas do VarAstro (var.astro.cz), sucessor do O-C
     Gateway de Brno, que herdou a base de tempos de minimo ("including
     records from the original O-C gateway"). Sai o nome, o numero de
     minimos e, se houver, o intervalo de anos dos minimos.

EXPECTATIVA (escrita antes de rodar): os 26 sao binarias de catalogo TESS,
SuperWASP mag 9-13, sem nome GCVS na maioria. Espero (a) minoria com
minimos na base de Brno - ate 5 dos 26 - e (b) nenhum ou ate 2 com trabalho
dedicado de O-C. Qualquer alvo com dP/dt publicado e comparado por numero
com a Tabela 3. Desvio em qualquer direcao vai para a nota como esta.

Saida: data/orquestra/oc_lote/literatura_26.parquet e o relato impresso.
"""
import http.cookiejar
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")  # titulos com letras gregas
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402

UA = {"User-Agent": "Mozilla/5.0 (nasa_discover; literatura_26)"}
SIMBAD_TAP = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
VARASTRO = "https://var.astro.cz"
RAIO_SIMBAD_ARCSEC = 5.0
RAIO_VARASTRO_ARCMIN = 1  # a API so aceita inteiro
PALAVRAS = re.compile(r"period|O-C|O−C|minim|timing|eclips|ephemer|binar", re.I)


def _tap(q):
    url = SIMBAD_TAP + "?" + urllib.parse.urlencode(
        {"request": "doQuery", "lang": "adql", "format": "json", "query": q})
    ult = None
    for tentativa in range(4):
        try:
            j = json.loads(urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90).read().decode())
            cols = [m["name"] for m in j["metadata"]]
            return pd.DataFrame(j["data"], columns=cols)
        except Exception as e:  # noqa: BLE001 - contado e repetido, nao escondido
            ult = e
            time.sleep(2.0 * (2 ** tentativa))
    raise RuntimeError(f"FALHA [simbad]: {ult}")


def simbad(ra, dec):
    b = _tap(f"SELECT b.oid, b.main_id, b.otype, b.nbref, "
             f"DISTANCE(POINT('ICRS', b.ra, b.dec), POINT('ICRS', {ra}, {dec}))*3600 AS sep "
             f"FROM basic b WHERE CONTAINS(POINT('ICRS', b.ra, b.dec), "
             f"CIRCLE('ICRS', {ra}, {dec}, {RAIO_SIMBAD_ARCSEC}/3600.))=1")
    if not len(b):
        return {"simbad_id": None, "simbad_otype": None, "simbad_nref": 0, "simbad_sep_arcsec": None,
                "simbad_ids": None, "refs": [], "refs_relevantes": [], "refs_dedicadas": []}
    b = b.sort_values("sep").iloc[0]
    ids = _tap(f"SELECT id FROM ident WHERE oidref = {int(b.oid)}")
    # "year" e palavra reservada no TAP: entre aspas e com alias
    refs = _tap('SELECT r.bibcode, r."year" AS yr, r.title, r.nbobject FROM has_ref h '
                f'JOIN ref r ON h.oidbibref = r.oidbib WHERE h.oidref = {int(b.oid)} ORDER BY yr')
    # relevante = titulo cita periodo/O-C/minimos E o trabalho nao e catalogo
    # de massa (nbobject = quantos objetos o SIMBAD liga ao artigo)
    rel = refs[refs.title.fillna("").str.contains(PALAVRAS) & (refs.nbobject.fillna(10**9) <= 100)]
    dedic = refs[refs.nbobject.fillna(10**9) <= 10]
    return {"simbad_id": b.main_id, "simbad_otype": b.otype, "simbad_nref": int(b.nbref),
            "simbad_sep_arcsec": float(b.sep), "simbad_ids": "; ".join(ids.id.tolist()),
            "refs": refs.to_dict("records"), "refs_relevantes": rel.to_dict("records"),
            "refs_dedicadas": dedic.to_dict("records")}


class VarAstro:
    """Token anonimo vem no cookie da propria pagina; a API exige Bearer."""

    def __init__(self):
        cj = http.cookiejar.CookieJar()
        self.op = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cj))
        self.op.open(urllib.request.Request(VARASTRO + "/en/Stars", headers=UA), timeout=60).read()
        tok = {c.name: c.value for c in cj}
        if "Token" not in tok:
            raise RuntimeError("FALHA [varastro]: pagina nao entregou o token anonimo")
        self.h = dict(UA, Authorization="Bearer " + tok["Token"], Accept="application/json")

    def api(self, path):
        ult = None
        for tentativa in range(4):
            try:
                r = self.op.open(urllib.request.Request(VARASTRO + path, headers=self.h), timeout=60)
                return json.loads(r.read().decode())
            except urllib.error.HTTPError as e:
                if e.code == 404:  # cone sem resultado devolve 404 com corpo '""'
                    return None
                ult = e
            except Exception as e:  # noqa: BLE001
                ult = e
            time.sleep(2.0 * (2 ** tentativa))
        raise RuntimeError(f"FALHA [varastro]: {path}: {ult}")

    def cone(self, ra, dec):
        j = self.api(f"/api/Search/Stars?pageId=1&pageSize=20&ra={ra:.6f}&dec={dec:.6f}"
                     f"&radiusArcmin={RAIO_VARASTRO_ARCMIN}")
        if not j or not j.get("data"):
            return None
        return sorted(j["data"], key=lambda d: d["distArcMin"] if d["distArcMin"] is not None else 1e9)[0]

    def minimos(self, star_id):
        j = self.api(f"/api/charts/ocdiagram/{star_id}")
        return j


def minimos_do_diagrama(j, tic, nome):
    """O diagrama O-C do VarAstro devolve, por minimo: JD, O-C (dias) contra a
    efemeride (m0, period) do proprio VarAstro, e o metodo. Vira uma linha
    por minimo, com o tipo (P/S) e a efemeride, para o cruzamento numerico."""
    if not j:
        return pd.DataFrame()
    linhas = []
    for tipo, chave in (("P", "dataPrimary"), ("S", "dataSecondary")):
        for d in j.get(chave) or []:
            linhas.append({"tic": tic, "va_nome": nome, "tipo": tipo, "jd": float(d["jd"]),
                           "oc_d": float(d["ocValue"]), "metodo": d.get("method"),
                           "P_va": float(j["period"]), "m0_va": float(j["m0"])})
    return pd.DataFrame(linhas)


def _resumo_minimos(m):
    if not len(m):
        return 0, 0, None, None
    ano = lambda jd: 2000.0 + (jd - 2451545.0) / 365.25  # noqa: E731
    return int((m.tipo == "P").sum()), int((m.tipo == "S").sum()), ano(m.jd.min()), ano(m.jd.max())


def cruzar(alvos):
    va = VarAstro()
    linhas, minimos = [], []
    for r in alvos.itertuples():
        s = simbad(float(r.ra), float(r.dec))
        v = va.cone(float(r.ra), float(r.dec))
        lin = {"tic": int(r.tic), "period_d": float(r.period_d), **{k: s[k] for k in s if not k.startswith("refs")},
               "simbad_titulos_relevantes": " || ".join(f"{x['yr']} {x['bibcode']} [{x['nbobject']} obj] {x['title']}" for x in s["refs_relevantes"]),
               "simbad_n_relevantes": len(s["refs_relevantes"]),
               "simbad_titulos_dedicados": " || ".join(f"{x['yr']} {x['bibcode']} [{x['nbobject']} obj] {x['title']}" for x in s["refs_dedicadas"]),
               "simbad_n_dedicados": len(s["refs_dedicadas"])}
        if v is None:
            lin.update({"va_id": None, "va_nome": None, "va_sep_arcmin": None, "va_periodo": None,
                        "va_tipo": None, "va_n_prim": 0, "va_n_sec": 0, "va_ano_ini": None, "va_ano_fim": None})
        else:
            m = minimos_do_diagrama(va.minimos(v["id"]), int(r.tic), v["name"])
            minimos.append(m)
            npri, nsec, a0, a1 = _resumo_minimos(m)
            lin.update({"va_id": int(v["id"]), "va_nome": v["name"], "va_sep_arcmin": v["distArcMin"],
                        "va_periodo": v.get("period"), "va_tipo": v.get("variabilityType"),
                        "va_n_prim": npri, "va_n_sec": nsec, "va_ano_ini": a0, "va_ano_fim": a1,
                        "va_crossids": "; ".join(str(c.get("name", c)) for c in v.get("crossIds", []))})
        linhas.append(lin)
        print(f"TIC {lin['tic']}: SIMBAD {lin['simbad_id']} ({lin['simbad_otype']}, {lin['simbad_nref']} refs, "
              f"{lin['simbad_n_relevantes']} com periodo/O-C no titulo e <=100 obj, {lin['simbad_n_dedicados']} com <=10 obj) | VarAstro "
              f"{lin['va_nome']} ({lin['va_n_prim']} min. primarios + {lin['va_n_sec']} secundarios"
              + (f", {lin['va_ano_ini']:.1f}-{lin['va_ano_fim']:.1f}" if lin.get("va_ano_ini") else "") + ")", flush=True)
        vistos = set()
        for x in s["refs_relevantes"] + s["refs_dedicadas"]:
            if x["bibcode"] in vistos:
                continue
            vistos.add(x["bibcode"])
            print(f"      {x['yr']} {x['bibcode']} [{x['nbobject']} obj]  {x['title']}", flush=True)
    return pd.DataFrame(linhas), (pd.concat(minimos, ignore_index=True) if minimos else pd.DataFrame())


if __name__ == "__main__":
    t26 = pd.read_parquet(config.DATA / "orquestra" / "oc_lote" / "tabela_26.parquet")
    alvos = pd.read_parquet(config.DATA / "orquestra" / "alvos_oc_88.parquet")
    alvos = alvos[alvos.tic.isin(t26.TIC)].reset_index(drop=True)
    assert len(alvos) == 26, len(alvos)
    df, minimos = cruzar(alvos)
    out = config.DATA / "orquestra" / "oc_lote" / "literatura_26.parquet"
    df.to_parquet(out, index=False)
    minimos.to_parquet(out.with_name("minimos_brno_26.parquet"), index=False)
    print("\n== resumo (expectativa: ate 5 com minimos em Brno; ate 2 com trabalho dedicado)")
    print(f"  no SIMBAD: {df.simbad_id.notna().sum()}/26; refs mediana {df.simbad_nref.median():.0f} "
          f"(min {df.simbad_nref.min()}, max {df.simbad_nref.max()})")
    print(f"  com titulo citando periodo/O-C/minimos (<=100 obj): {(df.simbad_n_relevantes > 0).sum()}/26; "
          f"com trabalho de <=10 objetos: {(df.simbad_n_dedicados > 0).sum()}/26")
    print(f"  no VarAstro: {df.va_id.notna().sum()}/26; com minimos primarios: {(df.va_n_prim > 0).sum()}/26; "
          f"com >= 4 primarios: {(df.va_n_prim >= 4).sum()}/26 ({len(minimos)} minimos guardados)")
    print(f"  -> {out}")
