# -*- coding: utf-8 -*-
"""Calibracao da cadeia do SuperWASP com JUPITERES QUENTES.

POR QUE TROCAR A CLASSE. A primeira tentativa usou binarias eclipsantes e
falhou por fisica, nao por medida: os tres calibradores que sobreviveram deram
O-C de +22,3, +71,0 e +81,8 min, que convertidos em dP/dt sao +0,051, +0,171 e
+0,202 s/ano - taxas tipicas de Algol com transferencia de massa. Binaria
eclipsante de periodo curto e JUSTAMENTE a classe onde o relogio anda. Extrair
um vies de ~1 min de baixo de 31,7 min de dispersao astrofisica exigiria mais
de mil calibradores.

Jupiter quente nao tem essa fonte de instabilidade. E o caso mais famoso de
decaimento orbital nao muda nada aqui: 30 ms/ano dao O-C de 0,005 min em 19
anos, contra a precisao de ~1 min que se busca. TTV de sistemas com companheira
tambem sao de segundos a poucos minutos e apareceriam como dispersao entre
setores do TESS - por isso o teste de estabilidade e INTERNO, pelo residuo do
ajuste linear multi-setor, e nao por curadoria de literatura.

O criterio de propagacao continua o mesmo: sigma_P vindo de N ciclos com epocas
de sigma_t custa sigma_t*M/N ao propagar M ciclos. Setores separados por anos
sao obrigatorios.
"""
import pathlib
import sys
import warnings

warnings.filterwarnings("ignore")
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import config
import ids as ids_mod
import superwasp as sw

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from astroquery.mast import Observations

pd.set_option("display.width", 240)

DEST = pathlib.Path(config.CACHE) / "hj_tess"
DEST.mkdir(parents=True, exist_ok=True)
OUT = config.RESULTS / "hj_calibracao.parquet"

ERR_PROP_MAX = 2.0
MIN_EPOCAS = 3
N_ALVOS = 12


def curvas_spoc(nome_alvo):
    o = Observations.query_criteria(objectname=nome_alvo, obs_collection="TESS",
                                    dataproduct_type="timeseries", radius="0.001 deg")
    if not len(o):
        return []
    pr = Observations.get_product_list(o).to_pandas()
    pr = pr[pr.productSubGroupDescription == "LC"]
    if not len(pr):
        return []
    m = Observations.download_products(Table.from_pandas(pr),
                                       download_dir=str(DEST), cache=True)
    out = []
    for p in m["Local Path"]:
        try:
            with fits.open(p) as h:
                s = int(h[0].header["SECTOR"])
                d = h[1].data
                ok = ((d["QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"])
                      & np.isfinite(d["TIME"]))
                if ok.sum() < 500:
                    continue
                t = np.asarray(d["TIME"][ok], dtype="<f8")
                f = np.asarray(d["PDCSAP_FLUX"][ok], dtype="<f8")
                out.append((s, t, f / np.percentile(f, 95)))
        except Exception as e:
            print("    erro em %s: %s" % (p, e), flush=True)
    return out


def main():
    c = pd.read_parquet("data/cache/superwasp/hj_calibradores.parquet").head(N_ALVOS)
    print("Jupiteres quentes a calibrar: %d" % len(c), flush=True)
    linhas = []
    for r in c.itertuples():
        P0, dur = float(r.P), float(r.dur_h)
        print("\n== %s   P %.6f d  profundidade %.2f%%  V %.2f"
              % (r.pl_name, P0, r.dep_pct, r.vmag), flush=True)
        rec = {"pl_name": r.pl_name, "P0": P0, "dep_pct": r.dep_pct,
               "dur_h": dur, "vmag": r.vmag, "swasp_id": r.swasp_id,
               "npts": int(r.npts), "snr_esperado": r.snr}
        try:
            cur = curvas_spoc(r.pl_name)
        except Exception as e:
            rec["status"] = "download falhou"
            linhas.append(rec)
            print("   download falhou: %s" % str(e)[:80], flush=True)
            continue
        if len(cur) < MIN_EPOCAS:
            rec["status"] = "menos de %d setores" % MIN_EPOCAS
            linhas.append(rec)
            print("   so %d setor(es)" % len(cur), flush=True)
            continue
        # DEDUPLICA POR SETOR ANTES DE ORDENAR. `cur.sort()` compara tuplas e,
        # quando dois produtos trazem o MESMO setor, cai na comparacao dos
        # arrays de tempo e levanta "truth value of an array is ambiguous".
        # Dois produtos do mesmo setor acontecem (reprocessamento, dois
        # segmentos), e sem a deduplicacao a mesma epoca entraria duas vezes no
        # ajuste, com peso dobrado.
        vistos, unicos = set(), []
        for s_, t_, f_ in cur:
            if s_ in vistos:
                continue
            vistos.add(s_)
            unicos.append((s_, t_, f_))
        cur = sorted(unicos, key=lambda x: x[0])
        rows = []
        for s, t, f in cur:
            t0, sg, pf = sw.medir_epoca(t, f, P0, sw.semente_local(t, f, P0),
                                        t14_h=dur)
            if np.isfinite(t0) and np.isfinite(sg) and 0 < sg < 8:
                rows.append({"setor": s, "t0": t0, "sig": sg, "prof": pf})
        if len(rows) < MIN_EPOCAS:
            rec["status"] = "menos de %d epocas boas" % MIN_EPOCAS
            linhas.append(rec)
            continue
        e = pd.DataFrame(rows)
        pr_ = e.prof.values
        if np.std(pr_) / max(abs(np.mean(pr_)), 1e-9) > 0.30:
            rec["status"] = "profundidade inconsistente"
            rec["profs"] = str(np.round(pr_ * 1e6).astype(int))
            linhas.append(rec)
            print("   profundidade inconsistente %s" % rec["profs"], flush=True)
            continue
        esc = sw.escada_de_periodo(e[["t0", "sig"]], P0)
        if esc is None:
            rec["status"] = "escada de ciclos nao fechou"
            linhas.append(rec)
            print("   escada nao fechou", flush=True)
            continue
        P, T0, sP, E, res = esc
        span = (e.t0.max() - e.t0.min()) / 365.25
        # ESTABILIDADE INTERNA: o residuo do ajuste linear multi-setor e o que
        # denuncia TTV ou deriva, sem depender de curadoria de literatura.
        chi2red = float(np.sum((res / e.sort_values("t0").sig.values) ** 2)
                        / max(len(e) - 2, 1))
        print("   %d epocas, span %.2f anos -> P %.8f +- %.4f s   res max %.2f min"
              "   chi2/dof %.1f"
              % (len(e), span, P, sP * 86400, np.max(np.abs(res)), chi2red),
              flush=True)
        rec.update({"n_epocas": len(e), "span_anos": span, "P": P,
                    "sP_s": sP * 86400, "res_max_min": float(np.max(np.abs(res))),
                    "chi2red_tess": chi2red,
                    "prof_tess_ppm": float(np.mean(pr_) * 1e6)})
        try:
            tw, fw, sig = sw.carregar(r.swasp_id, r.ra, r.dec)
        except Exception as ex:
            rec["status"] = "swasp falhou"
            linhas.append(rec)
            print("   swasp falhou: %s" % str(ex)[:80], flush=True)
            continue
        tmed = float(np.median(tw))
        M = abs(np.round((tmed - (T0 + 2457000)) / P))
        err_prop = sP * M * 1440
        rec.update({"ciclos_propagados": int(M), "err_prop_min": err_prop,
                    "t_medio_swasp": tmed, "n_swasp": len(tw),
                    "sig_fot_mmag": sig * 1000})
        if err_prop > ERR_PROP_MAX:
            rec["status"] = "propagacao imprecisa"
            linhas.append(rec)
            print("   propagacao %.2f min > %.1f - descartado"
                  % (err_prop, ERR_PROP_MAX), flush=True)
            continue
        t0w, sgw, ssub, subs, st_rv = sw.epoca_por_subconjunto(
            tw, fw, P, T0 + 2457000, t14_h=dur)
        if not np.isfinite(t0w):
            rec["status"] = "epoca swasp nao medivel"
            linhas.append(rec)
            print("   epoca do swasp nao medivel", flush=True)
            continue
        # POLITICA: ruido vermelho NAO TESTAVEL nao entra na calibracao.
        # Cair na barra formal aqui seria usar uma barra cuja premissa - a
        # ausencia de correlacao - nao foi verificada. E a correcao do gate de
        # cobertura AUMENTA a frequencia deste caso, entao sem esta politica
        # ela poderia MELHORAR barras por nao conseguir testa-las.
        rec["status_ruido_vermelho"] = st_rv
        if not np.isfinite(ssub):
            rec["status"] = "ruido vermelho nao testavel"
            linhas.append(rec)
            print("   %s - descartado" % st_rv, flush=True)
            continue
        kw = np.round((t0w - (T0 + 2457000)) / P)
        oc = (t0w - (T0 + 2457000 + kw * P)) * 1440
        sig_int = max(sgw, ssub)
        rec.update({"status": "ok", "oc_min": oc, "sig_formal_min": sgw,
                    "sig_entre_temporadas_min": ssub, "sig_interna_min": sig_int,
                    "sig_total_min": float(np.hypot(sig_int, err_prop)),
                    "n_temporadas": len(subs),
                    "desvios_temporadas": str([round(s["desvio_min"], 2) for s in subs]),
                    "datas_temporadas": str([round(s["t_medio"], 1) for s in subs]),
                    "prof_swasp_ppm": float(np.mean([s["prof"] for s in subs]) * 1e6)})
        print("   O-C = %+.2f min   formal %.2f  temporadas %s  propag %.2f"
              "  -> total %.2f"
              % (oc, sgw, ("%.2f" % ssub) if np.isfinite(ssub) else "n/a",
                 err_prop, rec["sig_total_min"]), flush=True)
        linhas.append(rec)

    d = pd.DataFrame(linhas)
    d.to_parquet(OUT, index=False)
    ok = d[d.status == "ok"] if "status" in d.columns else d.iloc[:0]
    print("\ngravado: %d   ok %d" % (len(d), len(ok)))
    if len(ok):
        print(ok[["pl_name", "oc_min", "sig_total_min", "chi2red_tess",
                  "t_medio_swasp"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()


# -----------------------------------------------------------------------------
# RESULTADO — a troca de classe resolveu o Passo 1
#
# Dos 12: 4 rejeitados (profundidade inconsistente, propagacao imprecisa, dois
# com menos de 3 setores) e 8 mediram. Destes, 4 passam a guarda de dispersao
# entre temporadas (barra < 5 min):
#
#   WASP-18 b   O-C  -2,85  +- 1,07     9 epocas, 7,81 anos, chi2/dof 1,03
#   WASP-4 b    O-C  -2,55  +- 1,52     7 epocas, 7,81 anos, chi2/dof 0,79
#   WASP-6 b    O-C  +0,85  +- 2,21     3 epocas, 5,01 anos, chi2/dof 0,19
#   WASP-19 b   O-C  -0,64  +- 2,85     7 epocas, 6,86 anos, chi2/dof 0,09
#
#   VIES DA CADEIA = -2,14 +- 0,78 min
#   chi2 2,62 / 3 dof = 0,87  -> consistente com offset CONSTANTE
#   dispersao bruta 1,73 min, contra 31,7 min dos calibradores de binaria
#
# A diferenca com o Passo 1 anterior e de uma ordem e meia de grandeza, e ela e
# inteiramente da classe: relogio de Jupiter quente nao anda, relogio de Algol
# anda. As datas medias cobrem 2454038 a 2454467 - 429 dias, abrangendo as duas
# temporadas do alvo (mediana 2454293) -, e o modelo constante e adequado, entao
# nao ha estrutura temporal exigida pelos dados.
#
# Os que a guarda condenou nao foram descartados em silencio: WASP-42 b saiu com
# O-C de -2.956 min e barra de 951, WASP-178 b com barra de 353, WASP-123 b com
# 16,7 e WASP-23 b com 23,5. A dispersao entre temporadas os marca como sem
# valor em vez de os deixar entrar com barra formal bonita.
#
# -----------------------------------------------------------------------------
# O PONTO DE 2006-2008 COM BARRA DEFENSAVEL
#
# CORRIGIDO EM 2026-09-12, depois de a cadeia virar codigo (`cadeia_oc.py`) e
# ser rodada contra este bloco como expectativa por componentes. O que estava
# escrito aqui e o que o codigo sustenta nao eram a mesma coisa, e num caso a
# diferenca inverte uma conclusao. O bloco antigo fica no git (2d46b3d^);
# este e o que vale.
#
#   t0 = 2454293,06732 BJD_TDB       E = -4699, 19,0 anos antes do s105
#   barra interna = 2,10 min  (dispersao entre as 2 temporadas; a formal da
#                              1,51, e a regra do desenho e "a maior vale")
#   barra = hypot(2,10 ; calibracao 0,78) = 2,24 min
#
#   residuo no LINEAR   dos 4 pontos TESS:  +22,72 +- 2,24 min  = 10,2 sigma
#   residuo na PARABOLA dos 4 pontos TESS:   -5,48 +- 2,24 min  =  2,4 sigma
#
# O QUE ESTAVA AQUI ANTES, E POR QUE CAIU. O bloco anterior dizia "residuo na
# parabola -16,04 = 7,9 sigma" e "a parabola previa +38,9 min para 2007". A
# parabola que o codigo ajusta aos 4 pontos TESS e A MESMA (chi2 2,04, Q
# 1,368e-9 - reproduzidos ao digito), mas avaliada em E = -4699 ela preve
# +28,2 min, nao +38,9. Testadas as convencoes possiveis: ponderada +28,2,
# nao ponderada +28,7, Q*E^2 puro +43,5 - NENHUMA da 38,9. O numero so
# existia neste comentario, sem codigo que o produzisse: era anotacao, nao
# resultado. E o mesmo padrao do docstring que descrevia a regra hibrida nao
# implementada, na direcao inversa - o comentario afirmando MAIS do que o
# codigo sustenta.
#
# O QUE CAI: "parabola pura rejeitada", e com ela a leitura de que a FORMA do
# O-C e inconsistente com dP/dt constante. A inferencia sobre senoide e
# terceiro corpo era construida inteiramente sobre a discrepancia 38,9 contra
# 22,9, e perde a base.
#
# O QUE FICA: o ponto de 2007 esta a 10 sigma da efemeride linear - variacao
# de periodo medida, robusta porque a alavanca e de 19 anos e a barra tem
# calibracao de cadeia propria. E ha curvatura dentro do TESS a p = 0,003
# (delta chi2 8,68, 1 dof). Isso NAO sustenta afirmacao sobre a forma da
# variacao.
#
# E O ADVERSARIAL ENFRAQUECE O QUE FICA, e por isso acompanha o numero em vez
# de ficar em nota separada: deslocando cada epoca TESS 1 sigma na direcao que
# mais enfraquece (s4 -0,51, s31 +0,30, s97 +0,30, s105 -0,30 min), delta chi2
# cai de 8,68 para 1,19 e p sobe de 0,003 para 0,28. A curvatura interna ao
# TESS depende de nenhuma das quatro epocas estar deslocada 1 sigma - o s4
# sozinho ja a derrubava. O elo forte e o ponto de 2007 contra o linear; o elo
# fraco e a curvatura dentro do TESS.
