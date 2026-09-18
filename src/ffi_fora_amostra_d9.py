# -*- coding: utf-8 -*-
"""Rodada oito, Passo J: as FFI posteriores ao ultimo setor da amostra (s86) como
teste FORA DA AMOSTRA para os 8 testaveis do D9 - e o controle do estimador de FFI.

1. INVENTARIO no MAST: para cada um dos 8, todos os produtos de serie temporal
   (HLSP TESS-SPOC e QLP, mais SPOC 2-min para referencia), com setor, datas e
   nome do produto. O corte da amostra e s86 (o ultimo setor do cache dos 26,
   2024,93). Se NAO houver nada depois de s86, o Passo J para aqui e diz.
2. MEDIDA (so se houver): epocas primarias com O MESMO estimador (`epocas_142874476.medir`,
   trapezio de profundidade livre, semente pela efemeride), barra =
   max(formal, |t0_a - t0_b| / 2 das metades), cadencia registrada.
3. COMPARACAO (so se houver): O-C dessas epocas contra tres modelos ajustados
   SO a 2004-2024,9 - quadratica de 7 epocas (Tabelas 3-5), quadratica-tudo
   (TESS completo + SuperWASP) e linear local de 2024 (setores 76-86) - com o
   chi2 de cada sob a covariancia do GP agrupado (Matern-3/2, A 0,835 min,
   tau_c 67,7 d), incluindo o termo do GP entre as epocas novas e as antigas.
4. PREVISAO POR CLASSE (declarada ANTES de rodar, Passo I):
   - "estavel" (329246824, 392536812, 230386284, 377253090): vence a
     QUADRATICA-TUDO; a de 7 epocas empata com ela dentro das barras (as duas
     concordam em dP/dt), e a linear local de 2024 perde por extrapolar sem
     curvatura. Previsao falsificavel: chi2/n da quadratica-tudo < 2 e menor
     que o da linear local.
   - "contraditado" (390021728, 232634196, 424461577): vence a LINEAR LOCAL DE
     2024, porque nesses o bloco TESS recente tem periodo proprio incompativel
     com a curvatura de 20 anos (em 232634196 a linear local ja coincide com a
     linear global; em 424461577 P_local - P_escada = +1,7 s contra +0,2 s
     exigidos; em 390021728 o dP/dt so-TESS tem sinal oposto). As duas
     quadraticas perdem, e a de 7 epocas perde por mais.
   - "ambiguo" (198408416): sem previsao - e o que o teste decide.
   Um alvo "estavel" em que a linear local vencer, ou um "contraditado" em que
   a quadratica de 7 epocas vencer, refuta a leitura do Passo I e vai ao
   relatorio como desvio.

EXPECTATIVA DO INVENTARIO (escrita e commitada ANTES de rodar; ja incorpora uma
SONDAGEM de orientacao feita em TIC 230386284, que mostrou o arquivo parando em
s86 para aquele alvo):
  Espero ZERO produtos depois de s86 nos 8 - os 8 estao em Draco (dec +55 a +75)
  e o TESS em 2025-2026 (setores 97+) aponta para o sul: o alvo de controle
  142874476 (dec -34) tem s97 e s105 no cache, e V Gru (dec -40) tem s105, mas
  nenhum alvo de Draco tem setor > 86. Ultimo setor esperado: s86 para os 8
  (2024,93). Se algum tiver s87+, o Passo J roda inteiro e isso e o resultado.
  Produtos por alvo esperados: QLP 200 s em s56-s86 (~16), TESS-SPOC 200 s no
  mesmo intervalo (~15, o TESS-SPOC as vezes atrasa um setor), QLP/TESS-SPOC
  600 s em s40-s55, 1800 s em s14-s26, SPOC 2-min em quase todos.

CONTROLE DO ESTIMADOR DE FFI (--controle; expectativa tambem commitada antes):
  Como o teste fora da amostra vai usar FFI quando os setores do norte voltarem,
  a maquinaria de 2 e 3 e exercitada agora nos setores que JA existem: as FFI de
  200 s de s84, s85 e s86 dos 8 alvos, medidas pelo mesmo estimador, contra as
  epocas de 2 min dos MESMOS setores (jitter_primarios_26.parquet, estimador do
  estado D). Esperado: |t0(FFI) - t0(2 min)| mediana 0,1-0,4 min e < 1 min em
  todas; barra formal da FFI 1,3-2,0x a de 2 min (200 s contra 120 s, mais o
  ruido de FFI); barra por metades da mesma ordem; QLP e TESS-SPOC concordando
  dentro de 0,5 min; cadencia lida 3,33 min. Um deslocamento SISTEMATICO (mesma
  direcao em todos os alvos) acima de 0,5 min significa que FFI e 2 min nao sao
  intercambiaveis no O-C e o teste fora da amostra precisa de um termo de
  offset - iria ao relatorio.

RESULTADO (2026-09-17, contra a expectativa acima; logs 34_passo_J_inventario.txt
e 35_passo_J_controle.txt):
  INVENTARIO: ZERO produtos depois de s86 nos 8, como esperado. Todos param em
  s86 (2024-12-18); 136-188 produtos por alvo; QLP e TESS-SPOC de 200 s cobrem
  s56-s86 (o TESS-SPOC entrega ate s85 em 7 dos 8: atrasa um setor), 600 s em
  s40-s55, 1800 s em s14-s26, SPOC 2-min em s14-s86. Dec +56 a +64: o TESS de
  2025-2026 aponta para o sul (o controle 142874476, dec -34, tem s97 e s105).
  Logo os itens 2, 3 e 4 NAO podem ser feitos hoje - o Passo J para aqui, com a
  previsao por classe ja registrada acima para quando os setores do norte
  voltarem. Desvio menor: apareceu uma reducao que eu nao tinha listado, TARS
  (200/600/1800 s), em 4 dos 8 alvos.
  CONTROLE do estimador (38 medidas, s84-s86, QLP e TESS-SPOC contra as epocas
  de 2 min dos MESMOS setores; cadencia lida 3,33 min, como esperado):
  - TESS-SPOC 200 s: |dt| mediana 0,30 min, max 1,25, media -0,06 +- 0,16 -
    dentro do esperado e SEM deslocamento sistematico; barra formal x1,33 (o
    esperado era 1,3-2,0). E utilizavel no teste fora da amostra.
  - QLP 200 s: |dt| mediana 0,48, max 7,20 (232634196 s85), media -0,64 - FORA
    do esperado ("< 1 min em todas"). E estrutura, nao ruido: em s85 o QLP fica
    NEGATIVO nos 7 alvos que o tem (media -2,02, mediana -1,20); tirando s85, a
    media do QLP volta a -0,05 +- 0,17 com |dt| mediana 0,30. E a barra formal
    do QLP e x1,00 a de 2 min - ou seja, nao cobre esse deslocamento.
  - 232634196 e o pior caso nos dois (QLP s85 -7,20 e s86 -3,76; TESS-SPOC
    +1,25 / -1,25): e o alvo de eclipse mais largo e P mais longo do conjunto.
  Leitura: quando os setores do norte voltarem, o teste fora da amostra deve
  usar TESS-SPOC 200 s; usar QLP exige conferir o setor contra 2 min antes, e
  um alvo sem 2 min no setor novo precisaria de um termo de offset declarado.
"""
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402
import deriva_gp_26 as GP  # noqa: E402
import download as DL  # noqa: E402
import epocas_142874476 as ep  # noqa: E402
import gls_vies_26 as G  # noqa: E402
import jitter_primarios_26 as J  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
CACHE_FFI = config.CACHE / "ffi_d9"
S_CORTE = 86          # ultimo setor da amostra (2024,93)
S_CONTROLE = (84, 85, 86)
D9_TESTAVEIS = [329246824, 392536812, 198408416, 230386284, 377253090, 390021728, 232634196, 424461577]
CLASSE_I = {329246824: "estavel", 392536812: "estavel", 230386284: "estavel", 377253090: "estavel",
            198408416: "ambiguo", 390021728: "contraditado", 232634196: "contraditado", 424461577: "contraditado"}
URL = "https://mast.stsci.edu/api/v0.1/Download/file?uri="


def setor_do_nome(s):
    m = re.search(r"[-_]s(\d{4})", str(s))
    return int(m.group(1)) if m else None


def mjd_para_data(mjd):
    from astropy.time import Time
    return Time(float(mjd), format="mjd").iso[:10] if np.isfinite(mjd) else ""


def inventario(tic, tentativas=3):
    """Todos os produtos de serie temporal do alvo no MAST (HLSP + TESS), com setor e datas."""
    from astroquery.mast import Observations
    ult = None
    for k in range(tentativas):
        try:
            o = Observations.query_criteria(target_name=str(tic), obs_collection=["HLSP", "TESS"],
                                            dataproduct_type="timeseries").to_pandas()
            if not len(o):
                return pd.DataFrame()
            o = o[["obs_collection", "provenance_name", "t_exptime", "obs_id", "obsid", "t_min", "t_max", "calib_level"]].copy()
            o["tic"] = int(tic)
            o["setor"] = o.obs_id.map(setor_do_nome)
            o["data_ini"] = o.t_min.map(mjd_para_data); o["data_fim"] = o.t_max.map(mjd_para_data)
            return o
        except Exception as e:  # noqa: BLE001 - contado e repetido, nunca lido como ausencia
            ult = e
            time.sleep(3.0 * (k + 1))
    raise RuntimeError(f"MAST falhou {tentativas}x para TIC {tic}: {str(ult)[:120]}")


def produtos_ffi(tics, setores, tentativas=3):
    """URLs das curvas de luz FFI (QLP e TESS-SPOC) dos setores pedidos."""
    from astroquery.mast import Observations
    linhas = []
    for tic in tics:
        inv = inventario(tic)
        alvo = inv[(inv.obs_collection == "HLSP") & (inv.provenance_name.isin(["QLP", "TESS-SPOC"]))
                   & (inv.t_exptime <= 200) & (inv.setor.isin(setores))]
        if not len(alvo):
            continue
        ult = None
        for k in range(tentativas):
            try:
                pr = Observations.get_product_list(Observations.query_criteria(obsid=[str(x) for x in alvo.obsid])).to_pandas()
                break
            except Exception as e:  # noqa: BLE001
                ult = e; time.sleep(3.0 * (k + 1))
        else:
            raise RuntimeError(f"MAST (produtos) falhou para TIC {tic}: {str(ult)[:120]}")
        for r in pr.itertuples():
            nome = Path(str(r.dataURI)).name
            if not nome.endswith(".fits") or ("_llc" not in nome and "_lc" not in nome):
                continue
            s = setor_do_nome(nome)
            if s not in setores:
                continue
            linhas.append({"tic": int(tic), "provenance": "QLP" if "_qlp_" in nome else "TESS-SPOC",
                           "setor": s, "filename": nome, "url": URL + str(r.dataURI)})
    return pd.DataFrame(linhas).drop_duplicates(["tic", "filename"]).reset_index(drop=True)


def ler_ffi(caminho):
    """Curva FFI (QLP ou TESS-SPOC) com a MESMA normalizacao de epocas_142874476._norm."""
    from astropy.io import fits
    with fits.open(caminho) as h:
        d = h[1].data
        cols = [c.upper() for c in d.columns.names]
        if "KSPSAP_FLUX" in cols:
            red, fl = "QLP", d["KSPSAP_FLUX"]
        elif "DET_FLUX" in cols:
            red, fl = "QLP", d["DET_FLUX"]
        elif "PDCSAP_FLUX" in cols:
            red, fl = "TESS-SPOC", d["PDCSAP_FLUX"]
        else:
            raise RuntimeError(f"FALHA [ffi]: {Path(caminho).name} sem coluna de fluxo reconhecida: {cols}")
        setor = int(h[0].header.get("SECTOR"))
        cad = float(h[1].header.get("TIMEDEL", np.nan)) * 1440.0
        q = d["QUALITY"] if "QUALITY" in d.columns.names else None
        t, f = ep._norm(d["TIME"], fl, q)
    return setor, red, cad, t, f


def epoca_ffi(t, f, P, t14_h, semente):
    """A epoca do setor com o estimador da cadeia e a barra max(formal, |a - b| / 2)."""
    sp = max(120.0, 1.5 * t14_h * 60)
    t0, sig, prof = ep.medir(t, f, P=P, t14_h=t14_h, semente=float(semente), span_min=sp)
    meio = int(np.searchsorted(t, np.median(t)))
    dif = disp = np.nan
    if meio > 100 and len(t) - meio > 100:
        ta, _, _ = ep.medir(t[:meio], f[:meio], P=P, t14_h=t14_h, semente=t0, span_min=sp)
        tb, _, _ = ep.medir(t[meio:], f[meio:], P=P, t14_h=t14_h, semente=t0, span_min=sp)
        ka, kb = np.round((ta - t0) / P), np.round((tb - t0) / P)
        dif = abs((ta - ka * P) - (tb - kb * P)) * 1440.0
        disp = dif / 2.0
    sig_uso = max(float(sig), float(disp)) if np.isfinite(disp) else float(sig)
    return {"t0_btjd": float(t0), "sigma_min": float(sig_uso), "sigma_formal_min": float(sig),
            "sigma_metades_min": float(disp), "dif_metades_min": float(dif), "prof_ppm": float(prof) * 1e6, "n": int(len(t))}


def modelos_do_alvo(tic, A, tc):
    """Os tres modelos ajustados SO a 2004-2024,9, com a covariancia do GP agrupado."""
    import caso_completo_d9 as H
    j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
    P, sP, T0, E_ult = J.escada_do_registro(j)
    S = pd.read_parquet(BASE / "jitter_primarios_26.parquet")
    d = J.classificar_setores(S[(S.tic == tic) & (S.minimo == "primario")], P, sP, T0, E_ult, 0.0)
    d = d[d.ok]
    tess = pd.DataFrame({"fonte": ["TESS s%d" % s for s in d.setor], "t0": d.t0_btjd.values, "sig_min": d.sig_min.values,
                         "E": d.E.values, "setor": d.setor.values, "ano": H.ano(d.t0_btjd.values)}).sort_values("t0").reset_index(drop=True)
    pts = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
    sw = pts[pts.fonte.str.startswith("SuperWASP")][["fonte", "t0", "sig_min", "E"]].copy(); sw["setor"] = -1; sw["ano"] = H.ano(sw.t0.values)
    tudo = pd.concat([sw, tess], ignore_index=True).sort_values("t0").reset_index(drop=True)
    E0 = float(np.mean(tess.E.values))
    local = tess[tess.setor.isin(H.LOCAIS)].reset_index(drop=True)
    mods = {"quad_tudo": H.ajustar(tudo, 2, "matern32", A, tc, P, E0),
            "linear_local_2024": H.ajustar(local, 1, "matern32", A, tc, P, E0)}
    # quadratica de 7 epocas: os `pontos` do registro (Tabelas 3-5), no mesmo E0
    des = pts[["fonte", "t0", "sig_min", "E"]].copy(); des["setor"] = -1; des["ano"] = H.ano(des.t0.values)
    mods["quad_7_epocas"] = H.ajustar(des, 2, "matern32", A, tc, P, E0)
    return P, E0, T0, tess, tudo, mods


def oc_contra_modelos(novas, P, E0, mods, A, tc, ref_pontos):
    """O-C das epocas novas contra cada modelo e o chi2 sob a covariancia do GP agrupado.

    A covariancia entre as epocas NOVAS inclui o termo do GP (todas sao TESS);
    o kernel tambem as liga as antigas, mas o modelo ja foi ajustado sem elas -
    o chi2 aqui e o da predicao, com C_nova = diag(sigma^2) + K(nova, nova).
    """
    E = novas.E.values.astype(float) - E0
    t = novas.t0.values.astype(float)
    Knn = GP.kernel("matern32", np.abs(np.subtract.outer(novas.t0.values, novas.t0.values)), A, tc) / 1440.0 ** 2
    out = {}
    for nome, m in mods.items():
        coef = np.asarray(m["coef"]); cov = np.asarray(m["cov"])
        X = np.vander(E, len(coef))
        pred = X @ coef
        r = t - pred
        C = np.diag((novas.sig_min.values / 1440.0) ** 2) + Knn + X @ cov @ X.T
        chi2 = float(r @ np.linalg.solve(C, r))
        out[nome] = {"oc_min": (r * 1440.0).tolist(), "oc_medio_min": float(np.mean(r) * 1440.0),
                     "chi2": chi2, "n": int(len(r)), "chi2_por_n": chi2 / len(r)}
    return out


if __name__ == "__main__":
    CACHE_FFI.mkdir(parents=True, exist_ok=True)
    hip = json.loads((BASE / "deriva_gp_26.json").read_text(encoding="utf-8"))["hiper_agrupados"]["matern32"]
    A, tc = hip["A_min"], hip["tau_c_d"]
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")

    print(f"== PASSO J.1: inventario MAST dos {len(D9_TESTAVEIS)} testaveis do D9; corte da amostra: s{S_CORTE}\n")
    invs, resumo = [], []
    for tic in D9_TESTAVEIS:
        inv = inventario(tic)
        invs.append(inv)
        g = inv.dropna(subset=["setor"])
        novos = g[(g.setor > S_CORTE) & (g.obs_collection == "HLSP") & (g.provenance_name.isin(["QLP", "TESS-SPOC"]))]
        por = g.groupby(["obs_collection", "provenance_name", "t_exptime"]).setor.agg(["min", "max", "count"])
        resumo.append({"tic": int(tic), "nome": str(t26.loc[tic, "nome"]), "dec": float(t26.loc[tic, "dec"]), "classe_I": CLASSE_I[tic],
                       "n_produtos": int(len(g)), "setor_max": int(g.setor.max()), "data_fim": g.loc[g.setor.idxmax(), "data_fim"],
                       "n_ffi_depois_corte": int(len(novos)), "setores_depois_corte": sorted(set(int(s) for s in novos.setor))})
        R = resumo[-1]
        print(f"  TIC {tic} {R['nome']:<18s} dec {R['dec']:+.1f}  {R['n_produtos']:3d} produtos, ultimo setor s{R['setor_max']} ({R['data_fim']}) | FFI depois de s{S_CORTE}: {R['n_ffi_depois_corte']}")
        print("     " + "; ".join(f"{i[1]} {int(i[2])} s: s{int(r['min'])}-s{int(r['max'])} ({int(r['count'])})" for i, r in por.iterrows()))
    INV = pd.concat(invs, ignore_index=True); INV.to_parquet(BASE / "ffi_fora_amostra_d9_inventario.parquet", index=False)
    RES = pd.DataFrame(resumo); RES.to_parquet(BASE / "ffi_fora_amostra_d9.parquet", index=False)
    n_novos = int(RES.n_ffi_depois_corte.sum())
    out = {"corte_setor": S_CORTE, "gp": {"A_min": A, "tau_c_d": tc}, "classe_I": {str(k): v for k, v in CLASSE_I.items()},
           "inventario": resumo, "n_ffi_depois_corte": n_novos}
    print(f"\n  TOTAL de produtos FFI depois de s{S_CORTE} nos 8 alvos: {n_novos}")
    if n_novos == 0:
        out["veredito"] = (f"NAO HA nada depois de s{S_CORTE} para nenhum dos 8: o teste fora da amostra do Passo J "
                           "nao pode ser feito hoje. Os 8 estao em Draco e o TESS voltou a apontar para o sul.")
        print(f"  -> {out['veredito']}")
        print("     (itens 2, 3 e 4 ficam declarados no docstring, para quando os setores do norte voltarem)")
    else:
        out["veredito"] = f"HA {n_novos} produtos FFI depois de s{S_CORTE}: rodar --medir"
        print(f"  -> {out['veredito']}")

    if "--controle" in sys.argv:
        print(f"\n== CONTROLE do estimador de FFI: setores {S_CONTROLE} (dentro da amostra) medidos na FFI de 200 s "
              "contra as epocas de 2 min dos MESMOS setores")
        prods = produtos_ffi(D9_TESTAVEIS, set(S_CONTROLE))
        print(f"   {len(prods)} curvas FFI a baixar ({prods.groupby('provenance').size().to_dict() if len(prods) else {}})")
        S2 = pd.read_parquet(BASE / "jitter_primarios_26.parquet")
        linhas = []
        for r in prods.itertuples():
            dest = CACHE_FFI / r.filename
            if not dest.exists():
                ok, nb, err = DL.download_one(r.url, dest)
                if not ok:
                    raise RuntimeError(f"FALHA [ffi]: download de {r.filename}: {err}")
            setor, red, cad, t, f = ler_ffi(dest)
            j = json.loads((BASE / f"oc__{r.tic}.json").read_text(encoding="utf-8"))
            # P da escada (o do registro) e a duracao do catalogo que a cadeia usou para o trapezio
            alvo = pd.read_parquet(config.DATA / "orquestra" / "alvos_oc_88.parquet").set_index("tic").loc[r.tic]
            P = float(j["P_escada_d"]); t14_h = float(alvo.t14_h)
            ref = S2[(S2.tic == r.tic) & (S2.minimo == "primario") & (S2.setor == setor)]
            if not len(ref):
                print(f"   TIC {r.tic} s{setor} {red}: sem epoca de 2 min no cache - pulado")
                continue
            e = epoca_ffi(t, f, P, t14_h, float(ref.t0_btjd.iloc[0]))
            # DOBRAR NO CICLO: `medir` recentra no minimo mais proximo da MEDIA da curva, e a curva FFI
            # nao comeca/termina no mesmo ponto que a de 2 min - a diferenca crua sai em ciclos inteiros.
            d_d = e["t0_btjd"] - float(ref.t0_btjd.iloc[0])
            n_ciclos = int(np.round(d_d / P))
            dt = (d_d - n_ciclos * P) * 1440.0
            linhas.append({"tic": int(r.tic), "setor": int(setor), "provenance": red, "cadencia_min": cad, **e,
                           "t0_2min_btjd": float(ref.t0_btjd.iloc[0]), "sig_2min_min": float(ref.sig_min.iloc[0]),
                           "sig_formal_2min_min": float(ref.sig_formal_min.iloc[0]), "dt_ffi_menos_2min_min": dt, "n_ciclos_dobrados": n_ciclos,
                           "razao_sigma_formal": e["sigma_formal_min"] / float(ref.sig_formal_min.iloc[0])})
            print(f"   TIC {r.tic} s{setor} {red:<10s} cad {cad:.2f} min: t0 - t0(2 min) = {dt:+.3f} min | barra FFI {e['sigma_min']:.3f} "
                  f"(formal {e['sigma_formal_min']:.3f}, metades {e['sigma_metades_min']:.3f}) contra 2 min {float(ref.sig_min.iloc[0]):.3f}")
        C = pd.DataFrame(linhas); C.to_parquet(BASE / "ffi_fora_amostra_d9_controle.parquet", index=False)
        out["controle"] = {"n": int(len(C)), "dt_mediana_min": float(C.dt_ffi_menos_2min_min.median()), "dt_abs_mediana_min": float(C.dt_ffi_menos_2min_min.abs().median()),
                           "dt_abs_max_min": float(C.dt_ffi_menos_2min_min.abs().max()), "dt_media_min": float(C.dt_ffi_menos_2min_min.mean()),
                           "dt_sd_min": float(C.dt_ffi_menos_2min_min.std(ddof=1)), "razao_sigma_formal_mediana": float(C.razao_sigma_formal.median()),
                           "cadencia_min": float(C.cadencia_min.median()), "por_provenance": {k: float(v) for k, v in C.groupby("provenance").dt_ffi_menos_2min_min.median().items()}}
        c = out["controle"]
        print(f"\n   {c['n']} medidas: dt mediana {c['dt_mediana_min']:+.3f} min (|dt| mediana {c['dt_abs_mediana_min']:.3f}, max {c['dt_abs_max_min']:.3f}); "
              f"media {c['dt_media_min']:+.3f} +- {c['dt_sd_min'] / max(np.sqrt(c['n']), 1):.3f} (sd {c['dt_sd_min']:.3f}); barra formal FFI/2 min x{c['razao_sigma_formal_mediana']:.2f}; "
              f"cadencia {c['cadencia_min']:.2f} min; por reducao {c['por_provenance']}")
    (BASE / "ffi_fora_amostra_d9.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(f"\n  -> {BASE / 'ffi_fora_amostra_d9.json'}, _inventario.parquet, .parquet"
          + (", _controle.parquet" if "--controle" in sys.argv else ""))
