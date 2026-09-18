# -*- coding: utf-8 -*-
"""Passo 1 da rodada 'quatro medidas': jitter entre setores dos PRIMARIOS
(e secundarios) dos 26, medido em todo setor SPOC 2-min do MAST fora do
desenho, com o MESMO estimador (`oc_lote.epocas_2min`, estado D). Os setores
extras sao validacao e medida de ruido; NUNCA entram no ajuste das tabelas.

Expectativa em `expectativas/jitter_primarios_26.md`, commit cf0ef4e,
escrita ANTES de rodar. Nada em cadeia_oc / oc_lote / superwasp muda.

Etapas (ver o .md): (1) manifesto completo dos 26 (799 curvas), desenho =
`oc_lote.manifestos()`; (2) por alvo, em paralelo, `epocas_2min` no manifesto
desenho + extras, primario (semente do catalogo) e secundario (+ P/2), um
JSON por alvo em jitter_26/ (`--medir`); controle embutido: os setores do
desenho reproduzem os `pontos` do registro; (3) contagem de ciclos com P e T0
da escada do desenho, criterio sP x |Delta E| < 0,25, fase |phi| < 0,1,
profundidade > 0; (4) sigma_j por alvo, primario e secundario, MV contra a
linear local ao bloco TESS estendido, IC por perfil, n extras, cobertura de
lag; informativo = >= 4 extras; (5) deriva do secundario sobre todos os
setores contra 5.3.1; (6) variograma agregado com bootstrap por alvo e ajuste
c + q tau (q_RW = 2 q); (7) inflacao branca das barras TESS do DESENHO por
hypot(sigma, sigma_j) e `reinflar_tess_26.refazer` - quantos de D11/D9
sobrevivem (`--resumo`).

Saidas em data/orquestra/oc_lote/: jitter_primarios_26.parquet (alvo x setor
x minimo), jitter_primarios_26_resumo.parquet (alvo), variograma_26.parquet,
jitter_inflacao_26.parquet, jitter_primarios_26.json.
"""
import json
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402
import oc_lote as L  # noqa: E402
import reinflar_tess_26 as RT  # noqa: E402
from jitter_232634196 import menos2lnL, perfil_sigma_j  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
SAIDA = BASE / "jitter_26"
N_PROC = 6
MIN_EXTRAS_INFORMATIVO = 4
TOL_CONTROLE_MIN = 0.05
LIMITE_FASE = 0.1
LIMITE_SPN = 0.25
BINS_DIAS = np.array([0, 30, 100, 182.625, 365.25, 730.5, 1095.75, 1461.0, 1826.25, 2557.0])
ROTULOS_BINS = ["0-30 d", "30-100 d", "100 d-0,5 a", "0,5-1 a", "1-2 a", "2-3 a", "3-4 a", "4-5 a", "5-7 a"]
N_BOOT = 2000
SEMENTE_BOOT = 20260916
SJ_NAO_FECHA = 10.0
EXPECTATIVA = "cf0ef4e"


def manifesto_alvo(tic, man_desenho, man_completo):
    d = man_desenho[man_desenho.tic == tic][["tic", "filename", "url", "sector"]].copy()
    d["extra"] = False
    c = man_completo[(man_completo.tic == tic) & ~man_completo.filename.isin(d.filename)][["tic", "filename", "url", "sector"]].copy()
    c["extra"] = True
    m = pd.concat([d, c], ignore_index=True).drop_duplicates("filename").reset_index(drop=True)
    m["sector"] = m.sector.astype(int)
    return m


def medir_alvo(tic):
    """Trabalhador: primario e secundario em todos os setores do alvo; grava jitter__{tic}.json."""
    t_ini = time.time()
    try:
        j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
        alvo = pd.read_parquet(config.DATA / "orquestra" / "alvos_oc_88.parquet").set_index("tic").loc[tic]
        man = manifesto_alvo(tic, L.manifestos(), pd.read_parquet(config.CATALOGS / "pendente" / "completo26_2min.parquet"))
        P0, t14, t0cat = float(alvo.period_d), float(alvo.t14_h), float(alvo.t0_btjd)
        extras = set(man.sector[man.extra].tolist())
        prim = L.epocas_2min(tic, P0, t14, man, t0cat)
        desc_p = prim.attrs.get("descartados", [])
        sec = L.epocas_2min(tic, P0, t14, man, t0cat + 0.5 * P0)
        desc_s = sec.attrs.get("descartados", [])
        # controle embutido: os setores do desenho reproduzem os pontos do registro
        ref = {int(f.split("s")[1]): t for f, t in zip((p["fonte"] for p in j["pontos"]), (p["t0"] for p in j["pontos"])) if f.startswith("TESS")}
        dif = {int(r.setor): abs(r.t0_btjd - ref[int(r.setor)]) * 1440.0 for r in prim.itertuples() if int(r.setor) in ref}
        linhas = []
        for minimo, df in (("primario", prim), ("secundario", sec)):
            for r in df.itertuples():
                linhas.append({"tic": int(tic), "minimo": minimo, "setor": int(r.setor), "extra": int(r.setor) in extras,
                               "t0_btjd": float(r.t0_btjd), "sig_min": float(r.sigma_min), "sig_formal_min": float(r.sigma_formal_min),
                               "sig_metades_min": float(r.sigma_metades_min), "prof_ppm": float(r.prof_ajustada_ppm), "n_pontos": int(r.n)})
        out = {"tic": int(tic), "status": "ok", "P0_catalogo": P0, "t14_h": t14, "t0_catalogo": t0cat,
               "n_setores_manifesto": int(len(man)), "n_extras_manifesto": int(len(extras)),
               "setores_desenho": sorted(int(s) for s in man.sector[~man.extra]), "setores_extras": sorted(extras),
               "controle_max_dif_min": float(max(dif.values())) if dif else None, "controle_n": len(dif),
               "descartados_primario": desc_p, "descartados_secundario": desc_s, "linhas": linhas, "segundos": time.time() - t_ini}
    except Exception as e:  # falha alta no relatorio, nao no pool
        out = {"tic": int(tic), "status": f"erro: {e}", "traceback": traceback.format_exc(), "segundos": time.time() - t_ini}
    SAIDA.mkdir(parents=True, exist_ok=True)
    (SAIDA / f"jitter__{tic}.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    return out["tic"], out["status"], out.get("controle_max_dif_min"), out.get("n_extras_manifesto"), out["segundos"]


def medir_todos(tics, n_proc=N_PROC):
    pend = [t for t in tics if not (SAIDA / f"jitter__{t}.json").exists()
            or json.loads((SAIDA / f"jitter__{t}.json").read_text(encoding="utf-8")).get("status") != "ok"]
    print(f"medir: {len(pend)} de {len(tics)} alvos pendentes, {n_proc} processos", flush=True)
    with ProcessPoolExecutor(max_workers=n_proc) as ex:
        futs = {ex.submit(medir_alvo, t): t for t in pend}
        for f in as_completed(futs):
            tic, st, ctrl, nex, seg = f.result()
            print(f"  TIC {tic}: {st} | extras {nex} | controle max |dif| {ctrl if ctrl is None else round(ctrl, 4)} min | {seg / 60:.1f} min", flush=True)


def escada_do_registro(j):
    """P, sP da escada e T0 = t0 - E P do ultimo TESS do registro."""
    P, sP = float(j["P_escada_d"]), float(j["sP_escada_d"])
    tess = [p for p in j["pontos"] if p["fonte"].startswith("TESS")]
    ult = max(tess, key=lambda p: p["t0"])
    return P, sP, float(ult["t0"] - ult["E"] * P), int(ult["E"])


def classificar_setores(df, P, sP, T0, E_ult, deslocamento):
    """E pelo P da escada; fase; criterios. `deslocamento` = 0 (primario) ou P/2 (secundario)."""
    d = df.copy()
    x = (d.t0_btjd - T0 - deslocamento) / P
    d["E"] = np.round(x).astype(int)
    d["fase"] = x - d.E
    d["sPxN"] = sP * np.abs(d.E - E_ult)
    d["ok_ciclo"] = d.sPxN < LIMITE_SPN
    d["ok_fase"] = d.fase.abs() < LIMITE_FASE
    d["ok_prof"] = d.prof_ppm > 0
    d["ok_finito"] = np.isfinite(d.t0_btjd) & np.isfinite(d.sig_min) & (d.sig_min > 0)
    d["ok"] = d.ok_ciclo & d.ok_fase & d.ok_prof & d.ok_finito
    d["motivo"] = np.where(d.ok, "", np.where(~d.ok_finito, "nao finito", np.where(~d.ok_prof, "profundidade <= 0",
                           np.where(~d.ok_ciclo, "ciclo ambiguo", "fase inconsistente"))))
    return d


def sigma_j_alvo(d):
    """MV contra a linear local ao bloco TESS estendido; residuos no MV."""
    ok = d[d.ok].sort_values("t0_btjd")
    E = ok.E.values.astype(float)
    t_min = (ok.t0_btjd.values - ok.t0_btjd.values[0]) * 1440.0
    s = ok.sig_min.values
    if len(ok) < 3:
        return None, None
    perfil = perfil_sigma_j(E, t_min, s)
    _, r = menos2lnL(E, t_min, s, perfil["sigma_j_min"])
    return perfil, pd.DataFrame({"setor": ok.setor.values, "t0_btjd": ok.t0_btjd.values, "sig_min": s, "res_min": r, "extra": ok.extra.values})


def menos2lnL_grau(E, t_min, sig_min, sj, grau):
    v = sig_min ** 2 + sj ** 2
    w = 1.0 / v
    X = np.vstack([E ** k for k in range(grau + 1)]).T
    beta = np.linalg.solve(X.T @ (X * w[:, None]), X.T @ (w * t_min))
    r = t_min - X @ beta
    return float(np.sum(r ** 2 / v + np.log(v))), r


def sigma_j_alvo_quad(d, grade=np.arange(0.0, 10.0 + 1e-9, 0.005)):
    """POS-HOC (sugestao do revisor interno, sem expectativa previa): sigma_j contra QUADRATICA local ao bloco
    TESS estendido, para separar o jitter da curvatura real dentro do bloco (232634196, 424461577, 198388252)."""
    ok = d[d.ok].sort_values("t0_btjd")
    if len(ok) < 4:
        return None
    E = ok.E.values.astype(float); t_min = (ok.t0_btjd.values - ok.t0_btjd.values[0]) * 1440.0; s = ok.sig_min.values
    L = np.array([menos2lnL_grau(E, t_min, s, sj, 2)[0] for sj in grade])
    i = int(np.argmin(L)); dL = L - L[i]; ok95 = np.where(dL <= 3.84)[0]
    return {"sigma_j_quad": float(grade[i]), "ic95_lo": float(grade[ok95[0]]), "ic95_hi": float(grade[ok95[-1]]) if ok95[-1] < len(grade) - 1 else np.nan}


def diagnostico_quad(tics):
    """Cenario C (pos-hoc): MV contra quadratica local (informativos), IC95 sup (nao informativos)."""
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    A = pd.read_parquet(BASE / "jitter_primarios_26_resumo.parquet").set_index("tic")
    S = pd.read_parquet(BASE / "jitter_primarios_26.parquet")
    sjC, linhas = {}, []
    for tic in tics:
        j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
        P, sP, T0, E_ult = escada_do_registro(j)
        dp = classificar_setores(S[(S.tic == tic) & (S.minimo == "primario")], P, sP, T0, E_ult, 0.0)
        q = sigma_j_alvo_quad(dp)
        inf = bool(A.loc[tic, "informativo"])
        hi = A.loc[tic, "sj_prim_ic95_hi"]
        hi = hi if np.isfinite(hi) else SJ_NAO_FECHA
        sjC[tic] = (q["sigma_j_quad"] if q else hi) if inf else hi
        linhas.append({"tic": tic, "informativo": inf, "sj_lin": float(A.loc[tic, "sj_prim"]), "sj_quad": q["sigma_j_quad"] if q else np.nan,
                       "sj_quad_ic95_lo": q["ic95_lo"] if q else np.nan, "sj_quad_ic95_hi": q["ic95_hi"] if q else np.nan, "sj_usado_C": sjC[tic]})
    Q = pd.DataFrame(linhas)
    print("== POS-HOC: sigma_j contra quadratica local ao bloco TESS estendido (sem expectativa previa)")
    print(Q.round(3).to_string(index=False))
    inf = Q[Q.informativo]
    print(f"   mediana sj_quad nos informativos {inf.sj_quad.median():.2f} min (linear {inf.sj_lin.median():.2f}); razao quad/lin mediana {(inf.sj_quad / inf.sj_lin.replace(0, np.nan)).median():.2f}")
    IC = inflar_e_refazer(t26, sjC, "C (pos-hoc): MV quadratica (informativos), IC95 sup linear (nao informativos)")
    IC.to_parquet(BASE / "jitter_inflacao_26_quad.parquet", index=False)
    Q.to_parquet(BASE / "jitter_sigma_j_quad.parquet", index=False)
    d11 = IC[IC.em_D11]; d9 = IC[IC.em_D9]
    print(f"   cenario C: D11 sobrevivem {int(d11.detecta.sum())} de 11 ({', '.join(str(t) for t in d11[d11.detecta].tic)}) | caem: {', '.join(str(t) for t in d11[~d11.detecta].tic)} | "
          f"D9 adequados {int(d9.detecta_adequado.sum())} de 9 | novas {int(IC[~IC.em_D11].detecta.sum())}")
    print(d11[["tic", "sigma_j_usado_min", "barra_tess_mediana_min", "dPdt", "sdPdt", "p_curv", "p_adv", "p_gof", "detecta"]].round(4).to_string(index=False))
    j = json.loads((BASE / "jitter_primarios_26.json").read_text(encoding="utf-8"))
    j["pos_hoc_quadratica"] = {"nota": "sugestao do revisor interno apos o Passo 1; sem expectativa previa", "sj_quad_mediana_informativos": float(inf.sj_quad.median()),
                               "D11_sobrevivem": int(d11.detecta.sum()), "quais": sorted(int(t) for t in d11[d11.detecta].tic), "D9_adequados": int(d9.detecta_adequado.sum()),
                               "sj_quad": {int(r.tic): (None if not np.isfinite(r.sj_quad) else round(float(r.sj_quad), 3)) for r in Q.itertuples()}}
    (BASE / "jitter_primarios_26.json").write_text(json.dumps(j, indent=1, default=float), encoding="utf-8")


def deriva_secundario(dp, ds, P):
    m = dp[dp.ok][["setor", "t0_btjd", "sig_min"]].merge(ds[ds.ok][["setor", "t0_btjd", "sig_min"]], on="setor", suffixes=("_p", "_s"))
    if len(m) < 3:
        return {"n_setores_d": int(len(m))}
    d = (m.t0_btjd_s - m.t0_btjd_p - 0.5 * P) % P
    d = np.where(d > 0.5 * P, d - P, d) * 1440.0
    sd = np.hypot(m.sig_min_p, m.sig_min_s).values
    ano = (m.t0_btjd_p.values - m.t0_btjd_p.values.mean()) / 365.25
    A = np.vstack([ano, np.ones(len(ano))]).T
    w = 1 / sd ** 2
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    beta = cov @ (A.T @ (w * d))
    return {"n_setores_d": int(len(m)), "d_mediano_min": float(np.median(d)), "deriva_todos_min_por_ano": float(beta[0]),
            "s_deriva_todos_min_por_ano": float(np.sqrt(cov[0, 0])), "alavanca_todos_anos": float(ano.max() - ano.min())}


def variograma(res_por_alvo, rng):
    """Pares do mesmo alvo; bins de lag; bootstrap por alvo; ajuste c + q tau (anos) no nivel dos pares."""
    pares = []
    for tic, R in res_por_alvo.items():
        t, r, s = R.t0_btjd.values, R.res_min.values, R.sig_min.values
        i, k = np.triu_indices(len(t), 1)
        pares.append(pd.DataFrame({"tic": tic, "tau_d": np.abs(t[k] - t[i]), "d2": (r[k] - r[i]) ** 2, "s2": s[i] ** 2 + s[k] ** 2}))
    Pp = pd.concat(pares, ignore_index=True)
    Pp["bin"] = np.digitize(Pp.tau_d, BINS_DIAS) - 1
    Pp["gnorm"] = Pp.d2 / Pp.s2
    Pp["gmin2"] = 0.5 * (Pp.d2 - Pp.s2)
    Pp["tau_a"] = Pp.tau_d / 365.25
    tics = np.array(sorted(res_por_alvo))
    grupos = {t: Pp[Pp.tic == t] for t in tics}

    def estat(sel):
        X = pd.concat([grupos[t] for t in sel], ignore_index=True)
        g = X.groupby("bin").agg(gnorm=("gnorm", "mean"), gmin2=("gmin2", "mean"), n=("gnorm", "size"), tau_med_d=("tau_d", "median"))
        g = g.reindex(range(len(ROTULOS_BINS)))
        A = np.vstack([np.ones(len(X)), X.tau_a.values]).T
        c, q = np.linalg.lstsq(A, X.gmin2.values, rcond=None)[0]
        return g.gnorm.values, g.gmin2.values, g.n.values, g.tau_med_d.values, float(c), float(q)
    gn, gm, n, tau_med, c0, q0 = estat(tics)
    B = [estat(rng.choice(tics, len(tics), replace=True)) for _ in range(N_BOOT)]
    gn_b = np.array([b[0] for b in B]); gm_b = np.array([b[1] for b in B]); cq_b = np.array([[b[4], b[5]] for b in B])
    bins = pd.DataFrame({"bin": ROTULOS_BINS, "tau_min_d": BINS_DIAS[:-1], "tau_max_d": BINS_DIAS[1:], "tau_mediano_d": tau_med, "n_pares": n,
                         "gama_norm": gn, "gama_norm_p16": np.nanpercentile(gn_b, 16, axis=0), "gama_norm_p84": np.nanpercentile(gn_b, 84, axis=0),
                         "gama_min2": gm, "gama_min2_p16": np.nanpercentile(gm_b, 16, axis=0), "gama_min2_p84": np.nanpercentile(gm_b, 84, axis=0)})
    ajuste = {"c_min2": c0, "q_min2_por_ano": q0, "c_p2.5": float(np.percentile(cq_b[:, 0], 2.5)), "c_p97.5": float(np.percentile(cq_b[:, 0], 97.5)),
              "q_p2.5": float(np.percentile(cq_b[:, 1], 2.5)), "q_p16": float(np.percentile(cq_b[:, 1], 16)), "q_p84": float(np.percentile(cq_b[:, 1], 84)),
              "q_p97.5": float(np.percentile(cq_b[:, 1], 97.5)), "q_sigma_boot": float(np.std(cq_b[:, 1], ddof=1)),
              "q_positivo_2sigma": bool(np.percentile(cq_b[:, 1], 2.5) > 0), "q_RW_min2_por_ano": 2 * q0, "n_alvos": int(len(tics)), "n_pares": int(len(Pp)), "n_boot": N_BOOT}
    return bins, ajuste, Pp


def inflar_e_refazer(t26, sj_por_alvo, rotulo):
    linhas = []
    for tic in sorted(t26.index):
        j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
        pontos = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
        P = float(j["P_escada_d"])
        sj = float(sj_por_alvo[tic])
        p = pontos.copy()
        tess = p.fonte.str.startswith("TESS")
        p.loc[tess, "sig_min"] = np.hypot(p.sig_min[tess], sj)
        r = RT.refazer(p, P)
        det = (r["p_curv"] < RT.P_DET) and (r["p_adv"] < RT.P_DET) and r["vao_ok"]
        linhas.append({"cenario": rotulo, "tic": int(tic), "sigma_j_usado_min": sj, "barra_tess_mediana_min": float(np.median(p.sig_min[tess])),
                       "barra_tess_mediana_original_min": float(np.median(pontos.sig_min[tess])), **r, "detecta": det,
                       "detecta_adequado": det and r["p_gof"] >= RT.P_DET, "em_D11": bool(t26.loc[tic, "curvatura"]),
                       "em_D9": bool(t26.loc[tic, "curvatura"] and not t26.loc[tic, "inadequada"]),
                       "dPdt_original": float(t26.loc[tic, "dPdt"]), "sdPdt_original": float(t26.loc[tic, "s"]), "p_curv_original": float(t26.loc[tic, "p"])})
    return pd.DataFrame(linhas)


def resumo(tics):
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    sec531 = pd.read_parquet(BASE / "secundarios_26.parquet").set_index("tic")
    setores, alvos, res_prim = [], [], {}
    for tic in tics:
        jj = json.loads((SAIDA / f"jitter__{tic}.json").read_text(encoding="utf-8"))
        if jj["status"] != "ok":
            alvos.append({"tic": tic, "status": jj["status"]}); print(f"TIC {tic}: {jj['status']}"); continue
        j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
        P, sP, T0, E_ult = escada_do_registro(j)
        df = pd.DataFrame(jj["linhas"])
        dp = classificar_setores(df[df.minimo == "primario"], P, sP, T0, E_ult, 0.0)
        ds = classificar_setores(df[df.minimo == "secundario"], P, sP, T0, E_ult, 0.5 * P)
        setores.append(dp); setores.append(ds)
        n_extras_ok = int((dp.ok & dp.extra).sum())
        informativo = n_extras_ok >= MIN_EXTRAS_INFORMATIVO
        pp, Rp = sigma_j_alvo(dp)
        ps, Rs = sigma_j_alvo(ds)
        lag = np.abs(np.subtract.outer(dp.t0_btjd[dp.ok].values, dp.t0_btjd[dp.ok].values))
        lag = lag[np.triu_indices_from(lag, 1)]
        der = deriva_secundario(dp, ds, P)
        a = {"tic": tic, "status": "ok", "P_escada_d": P, "sP_escada_d": sP, "n_setores": int(len(dp)), "n_extras": int(dp.extra.sum()),
             "n_extras_ok": n_extras_ok, "n_ok_primario": int(dp.ok.sum()), "n_ok_secundario": int(ds.ok.sum()),
             "n_falhas_primario": int((~dp.ok).sum()), "n_falhas_secundario": int((~ds.ok).sum()),
             "falhas_primario": "; ".join(f"s{r.setor}: {r.motivo}" for r in dp[~dp.ok].itertuples()),
             "falhas_secundario": "; ".join(f"s{r.setor}: {r.motivo}" for r in ds[~ds.ok].itertuples()),
             "sPxN_max": float(dp.sPxN.max()), "controle_max_dif_min": jj["controle_max_dif_min"], "controle_n": jj["controle_n"],
             "informativo": informativo, "lag_min_d": float(lag.min()) if len(lag) else np.nan, "lag_max_d": float(lag.max()) if len(lag) else np.nan,
             "barra_tess_mediana_min": float(dp.sig_min[dp.ok].median()), "segundos": jj["segundos"]}
        for rot, pf in (("prim", pp), ("sec", ps)):
            if pf is None:
                a.update({f"sj_{rot}": np.nan, f"sj_{rot}_ic68_lo": np.nan, f"sj_{rot}_ic68_hi": np.nan, f"sj_{rot}_ic95_lo": np.nan, f"sj_{rot}_ic95_hi": np.nan, f"sj_{rot}_n": 0})
            else:
                a.update({f"sj_{rot}": pf["sigma_j_min"], f"sj_{rot}_ic68_lo": pf["ic68"][0], f"sj_{rot}_ic68_hi": pf["ic68"][1],
                          f"sj_{rot}_ic95_lo": pf["ic95"][0], f"sj_{rot}_ic95_hi": pf["ic95"][1] if pf["ic95_fecha"] else np.nan,
                          f"sj_{rot}_ic95_fecha": pf["ic95_fecha"], f"sj_{rot}_n": pf["n"], f"sj_{rot}_dof": pf["dof_linear"],
                          f"sj_{rot}_exigido95": pf["ic95"][0] > 0, f"sj_{rot}_delta2lnL_0": pf["menos2lnL_em_0"] - pf["menos2lnL_min"]})
        a.update(der)
        if tic in sec531.index and "deriva_min_por_ano" in sec531.columns and np.isfinite(sec531.loc[tic].get("deriva_min_por_ano", np.nan)):
            a["deriva_531_min_por_ano"] = float(sec531.loc[tic, "deriva_min_por_ano"]); a["s_deriva_531"] = float(sec531.loc[tic, "s_deriva_min_por_ano"])
        alvos.append(a)
        if informativo and Rp is not None:
            res_prim[tic] = Rp
        print(f"TIC {tic}: {a['n_setores']} setores ({a['n_extras']} extras, {n_extras_ok} ok) | falhas prim {a['n_falhas_primario']} sec {a['n_falhas_secundario']} | "
              f"sPxN max {a['sPxN_max']:.3f} | controle {a['controle_max_dif_min']:.4f} min | "
              f"sigma_j prim {a['sj_prim']:.2f} [{a['sj_prim_ic95_lo']:.2f}, {a['sj_prim_ic95_hi']:.2f}] sec {a['sj_sec']:.2f} [{a['sj_sec_ic95_lo']:.2f}, {a['sj_sec_ic95_hi']:.2f}] "
              f"| lag {a['lag_min_d']:.0f}-{a['lag_max_d']:.0f} d | {'informativo' if informativo else 'NAO informativo'}"
              + (f" | deriva sec todos {a['deriva_todos_min_por_ano']:+.2f}+-{a['s_deriva_todos_min_por_ano']:.2f} (5.3.1 {a.get('deriva_531_min_por_ano', np.nan):+.2f}+-{a.get('s_deriva_531', np.nan):.2f})" if "deriva_todos_min_por_ano" in a else ""), flush=True)
    S = pd.concat(setores, ignore_index=True)
    A = pd.DataFrame(alvos)
    S.to_parquet(BASE / "jitter_primarios_26.parquet", index=False)
    A.to_parquet(BASE / "jitter_primarios_26_resumo.parquet", index=False)
    inf = A[A.informativo.fillna(False)]
    med_prim = float(inf.sj_prim.median()); med_sec = float(inf.sj_sec.median())
    print(f"\n== sigma_j PRIMARIO nos {len(inf)} informativos: mediana {med_prim:.2f} min (IQR {inf.sj_prim.quantile(0.25):.2f}-{inf.sj_prim.quantile(0.75):.2f}; "
          f"min {inf.sj_prim.min():.2f}, max {inf.sj_prim.max():.2f}); exigido a 95% em {int(inf.sj_prim_exigido95.sum())} de {len(inf)}")
    print(f"== sigma_j SECUNDARIO: mediana {med_sec:.2f} min; razao sec/prim mediana {float((inf.sj_sec / inf.sj_prim).median()):.2f}")
    print(f"== REGRA DE PARADA (mediana do primario > 3 min): {'DISPARA' if med_prim > 3 else 'nao dispara'}")
    print(f"== falhas: primario {int(A.n_falhas_primario.sum())} de {int(A.n_setores.sum())} setores; secundario {int(A.n_falhas_secundario.sum())}; "
          f"sPxN max {A.sPxN_max.max():.3f}; controle max |dif| {A.controle_max_dif_min.max():.4f} min em {int(A.controle_n.sum())} epocas do desenho")
    # variograma
    rng = np.random.default_rng(SEMENTE_BOOT)
    bins, ajuste, Pp = variograma(res_prim, rng)
    bins.to_parquet(BASE / "variograma_26.parquet", index=False)
    print("\n== variograma agregado dos residuos primarios (informativos):")
    print(bins.round(3).to_string(index=False))
    print(f"   ajuste gama_min2 = c + q tau: c {ajuste['c_min2']:.3f} min^2 [{ajuste['c_p2.5']:.3f}, {ajuste['c_p97.5']:.3f}], "
          f"q {ajuste['q_min2_por_ano']:.3f} min^2/ano [{ajuste['q_p2.5']:.3f}, {ajuste['q_p97.5']:.3f}] (boot por alvo, {ajuste['n_alvos']} alvos, {ajuste['n_pares']} pares); "
          f"q > 0 a 2 sigma: {ajuste['q_positivo_2sigma']} -> Passo 2 {'RODA' if ajuste['q_positivo_2sigma'] else 'nao roda'}; q_RW = 2q = {ajuste['q_RW_min2_por_ano']:.3f}")
    # inflacao
    sjA, sjB = {}, {}
    for r in A.itertuples():
        hi = r.sj_prim_ic95_hi if np.isfinite(r.sj_prim_ic95_hi) else SJ_NAO_FECHA
        sjA[r.tic] = r.sj_prim if r.informativo else hi
        sjB[r.tic] = hi
    IA = inflar_e_refazer(t26, sjA, "A: MV (informativos), IC95 sup (nao informativos)")
    IB = inflar_e_refazer(t26, sjB, "B: IC95 sup em todos")
    I = pd.concat([IA, IB], ignore_index=True)
    I.to_parquet(BASE / "jitter_inflacao_26.parquet", index=False)
    print("\n== inflacao branca das barras TESS do desenho por hypot(sigma, sigma_j):")
    for cen, g in I.groupby("cenario", sort=False):
        d11 = g[g.em_D11]; d9 = g[g.em_D9]
        print(f"  {cen}: sigma_j mediano usado {g.sigma_j_usado_min.median():.2f} min, barra TESS mediana {g.barra_tess_mediana_original_min.median():.2f} -> {g.barra_tess_mediana_min.median():.2f} min | "
              f"D11 sobrevivem {int(d11.detecta.sum())} de 11 ({', '.join(str(t) for t in d11[d11.detecta].tic)}) | caem: {', '.join(str(t) for t in d11[~d11.detecta].tic)} | "
              f"D9 com ajuste adequado {int(d9.detecta_adequado.sum())} de 9 | novas fora de D11 {int(g[~g.em_D11].detecta.sum())} | vao ambiguo {int((~g.vao_ok).sum())} | sPxN max {g.sP_vezes_N.max():.3f}")
        print(g[g.em_D11][["tic", "sigma_j_usado_min", "barra_tess_mediana_original_min", "barra_tess_mediana_min", "dPdt", "sdPdt", "p_curv", "p_adv", "p_gof", "detecta"]].round(4).to_string(index=False))
    out = {"expectativa_commit": EXPECTATIVA, "n_alvos": len(A), "n_informativos": int(len(inf)), "informativos": sorted(int(t) for t in inf.tic),
           "sigma_j_prim_mediana_min": med_prim, "sigma_j_sec_mediana_min": med_sec, "regra_parada_3min_dispara": med_prim > 3,
           "falhas_primario": int(A.n_falhas_primario.sum()), "falhas_secundario": int(A.n_falhas_secundario.sum()), "n_setores_total": int(A.n_setores.sum()),
           "controle_max_dif_min": float(A.controle_max_dif_min.max()), "variograma_ajuste": ajuste,
           "inflacao": {cen: {"D11_sobrevivem": int(g[g.em_D11].detecta.sum()), "quais": sorted(int(t) for t in g[g.em_D11 & g.detecta].tic),
                              "D9_adequados": int(g[g.em_D9].detecta_adequado.sum()), "novas": int(g[~g.em_D11].detecta.sum()),
                              "sigma_j_mediano": float(g.sigma_j_usado_min.median())} for cen, g in I.groupby("cenario", sort=False)}}
    (BASE / "jitter_primarios_26.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(f"\n  -> jitter_primarios_26.parquet, _resumo.parquet, variograma_26.parquet, jitter_inflacao_26.parquet, jitter_primarios_26.json em {BASE}")


if __name__ == "__main__":
    t26 = pd.read_parquet(BASE / "tabela_26.parquet")
    tics = sorted(int(t) for t in t26.TIC)
    assert len(tics) == 26
    if "--medir" in sys.argv:
        n = int(sys.argv[sys.argv.index("--medir") + 1]) if len(sys.argv) > sys.argv.index("--medir") + 1 and sys.argv[sys.argv.index("--medir") + 1].isdigit() else N_PROC
        medir_todos(tics, n)
    if "--resumo" in sys.argv:
        resumo(tics)
    if "--quad" in sys.argv:
        diagnostico_quad(tics)
