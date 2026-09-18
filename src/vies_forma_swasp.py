# -*- coding: utf-8 -*-
"""Passo 4 da rodada 'quatro medidas': o vies do estimador SuperWASP por FORMA
do minimo, medido com formas sinteticas sobre instantes, cadencia, temporadas
e ruido REAIS, pelo `superwasp.py` real (nada nele muda).

Expectativa em `expectativas/vies_forma_swasp.md`, commit a9088a8, escrita
ANTES de rodar.

Curvas: os 4 calibradores (hj_calibradores + hj_calibracao) e os 26 alvos.
`superwasp.carregar` faz o clip assimetrico, HJD_UTC -> BJD_TDB e a
normalizacao uma vez por curva; a conversao de tempo e deterministica nos
instantes e igual para verdade e medida, logo aplicada uma vez e exatamente
o mesmo que passar cada realizacao pelo CSV. Residuos reais = fluxo - curva
media em fase (mediana em 200 bins, P da verdade e a epoca REAL de cada temporada
para o ponto - senao o eclipse deslocado da temporada fica no residuo); bootstrap em blocos
de NOITE (lacunas > 0,5 d): cada noite real recebe o bloco de residuos de
uma noite sorteada (encadeada se mais curta). Formas, todas simetricas:
(i) transito Mandel-Agol (batman, LD quadratico 0,45/0,20 declarado) com
rp/rs, a/R*, i publicados (hj_transito_params.parquet) nas 4 curvas HJ;
(ii) ocultacao com escurecimento de bordo (batman, i = 90, rp/rs =
sqrt(depth), a/R* = (1 + k)/sin(pi T14/P)) nas 26; (iii) contato: elipsoidal
1 - A_ell (1 + cos 4 pi phi)/2 (maximos em quadratura; o .md escreveu
1 - cos, um deslize sem efeito na simetria), A_ell = 0,5 x depth, mais
eclipses de perfil cosseno (depth e 0,8 depth) de duracao T14, nas 26.
Verdade: efemeride linear (P da escada / P do calibrador, T0 = epoca global
real medida na curva real); t0_previsto = a verdade. 500 realizacoes por
curva e forma. Vies = epoca ajustada - verdadeira (global e por temporada).

Controle (`--controle`): a curva REAL pelo mesmo caminho reproduz as epocas
por temporada dos `pontos` (26) e os desvios/barras de hj_calibracao (4), e
- com a escada TESS refeita como em run_hj_calib (cache MAST) - o O-C dos 4.
"""
import json
import sys
import time
import traceback
import zlib
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cadeia_oc as co  # noqa: E402
import config  # noqa: E402
import superwasp as sw  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
SAIDA = BASE / "vies_forma_26"
N_REAL = 500
SEMENTE = 20260917
N_PROC = 10
U_LD = (0.45, 0.20)
N_BINS_FASE = 200
GAP_NOITE_D = 0.5
EXPECTATIVA = "a9088a8"
TOL_CONTROLE_MIN = 0.05
HJ = ("WASP-18 b", "WASP-4 b", "WASP-19 b", "WASP-6 b")
SD_INSTAVEL_MIN = 10.0


# ----------------------------------------------------------------- curvas
def curva_alvo(tic):
    alvo = pd.read_parquet(config.DATA / "orquestra" / "alvos_oc_88.parquet").set_index("tic").loc[tic]
    j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
    pts = pd.DataFrame(j["pontos"]).sort_values("t0")
    P = float(j["P_escada_d"])
    ult = pts[pts.fonte.str.startswith("TESS")].iloc[-1]
    T0_t = float(ult.t0 - ult.E * P)                       # ~ intercepto da escada (residuo < 1 min; so centra a busca)
    t, f, _ = sw.carregar(alvo.sourceid, float(alvo.ra), float(alvo.dec))
    t0_prev = T0_t + 2457000.0 + np.round((np.median(t) - T0_t - 2457000.0) / P) * P
    sw_pts = pts[pts.fonte.str.startswith("SuperWASP")]
    return {"nome": f"TIC {tic}", "tipo": "alvo", "t": t, "f": f, "P": P, "t14_h": float(alvo.t14_h), "depth": float(alvo.depth_frac),
            "morph": float(alvo.morph), "t0_prev": float(t0_prev),
            "ref_temporadas_bjd": (sw_pts.t0.values + 2457000.0 + co.VIES_CADEIA_MIN / 1440.0).tolist(),
            "ref_formal_min": sw_pts.sig_formal_min.values.tolist()}


def curva_hj(nome):
    cal = pd.read_parquet(config.DATA / "cache" / "superwasp" / "hj_calibradores.parquet").set_index("pl_name").loc[nome]
    res = pd.read_parquet(config.RESULTS / "hj_calibracao.parquet").set_index("pl_name").loc[nome]
    geo = pd.read_parquet(config.CATALOGS / "hj_transito_params.parquet").set_index("pl_name").loc[nome]
    assert res.status == "ok", (nome, res.status)
    P = float(res.P)                                        # o P da escada TESS do calibrador
    t, f, _ = sw.carregar(cal.swasp_id, float(cal.ra), float(cal.dec))
    tm = float(cal.tranmid)
    t0_prev = tm + np.round((np.median(t) - tm) / P) * P
    return {"nome": nome, "tipo": "hj", "t": t, "f": f, "P": P, "t14_h": float(cal.dur_h), "depth": float(cal.dep_pct) / 100.0,
            "k": float(geo.pl_ratror), "aRs": float(geo.pl_ratdor), "inc": float(geo.pl_orbincl), "t0_prev": float(t0_prev),
            "ref_oc_min": float(res.oc_min), "ref_formal_min": float(res.sig_formal_min), "ref_disp_min": float(res.sig_entre_temporadas_min),
            "ref_desvios": json.loads(res.desvios_temporadas), "ref_sig_total_min": float(res.sig_total_min), "ra": float(cal.ra), "dec": float(cal.dec)}


def carregar_curva(chave):
    return curva_hj(chave) if chave in HJ else curva_alvo(int(chave))


# ----------------------------------------------------------------- formas
def forma_transito(t, t0, P, k, aRs, inc):
    import batman
    p = batman.TransitParams()
    p.t0, p.per, p.rp, p.a, p.inc, p.ecc, p.w = t0, P, k, aRs, inc, 0.0, 90.0
    p.u, p.limb_dark = list(U_LD), "quadratic"
    return batman.TransitModel(p, np.asarray(t, float)).light_curve(p)


def forma_destacado(t, t0, P, depth, t14_h):
    k = float(np.sqrt(depth))
    aRs = (1.0 + k) / np.sin(np.pi * (t14_h / 24.0) / P)
    return forma_transito(t, t0, P, k, aRs, 90.0)


def forma_contato(t, t0, P, depth, t14_h):
    phi = ((np.asarray(t, float) - t0) / P) % 1.0
    A_ell = 0.5 * depth
    f = 1.0 - A_ell * (1.0 + np.cos(4 * np.pi * phi)) / 2.0
    w = (t14_h / 24.0) / P                                   # duracao em fase
    for centro, prof in ((0.0, depth), (0.5, 0.8 * depth)):
        d = (phi - centro + 0.5) % 1.0 - 0.5
        dentro = np.abs(d) < w / 2
        f[dentro] *= 1.0 - prof * (1.0 + np.cos(2 * np.pi * d[dentro] / w)) / 2.0
    return f


def forma(c, rot, t0_verdade):
    if rot == "transito":
        return forma_transito(c["t"], t0_verdade, c["P"], c["k"], c["aRs"], c["inc"])
    if rot == "destacado":
        return forma_destacado(c["t"], t0_verdade, c["P"], c["depth"], c["t14_h"])
    if rot == "contato":
        return forma_contato(c["t"], t0_verdade, c["P"], c["depth"], c["t14_h"])
    raise ValueError(rot)


# ----------------------------------------------------------------- residuos e bootstrap
def t0_por_ponto(t, t0g, subs, tol_d=200.0):
    """Cada ponto recebe a epoca da SUA temporada real (a mais proxima em t_medio); fora de qualquer
    temporada medida, a epoca global. Sem isso o residuo de uma temporada deslocada (232634196: -26 e
    +12 min) carrega o proprio eclipse e o bootstrap injeta esses 'dips' na sintetica, com vies."""
    t0 = np.full(len(t), float(t0g))
    if subs:
        tm = np.array([s["t_medio"] for s in subs]); ts = np.array([s["t0"] for s in subs])
        i = np.argmin(np.abs(t[:, None] - tm[None, :]), axis=1)
        perto = np.abs(t - tm[i]) < tol_d
        t0[perto] = ts[i[perto]]
    return t0


def residuos_reais(t, f, P, t0):
    """t0: escalar ou vetor (uma epoca por ponto)."""
    phi = ((t - t0) / P) % 1.0
    b = np.minimum((phi * N_BINS_FASE).astype(int), N_BINS_FASE - 1)
    med = np.array([np.median(f[b == i]) if np.any(b == i) else np.nan for i in range(N_BINS_FASE)])
    med = pd.Series(med).interpolate(limit_direction="both").values
    return f - med[b]


def noites(t):
    o = np.argsort(t)
    corte = np.where(np.diff(t[o]) > GAP_NOITE_D)[0]
    return [o[g] for g in np.split(np.arange(len(t)), corte + 1)]


def bootstrap_blocos(r, blocos, rng):
    out = np.empty_like(r)
    nb = len(blocos)
    for g in blocos:
        n = len(g)
        acc = []
        while sum(len(a) for a in acc) < n:
            acc.append(r[blocos[rng.integers(nb)]])
        out[g] = np.concatenate(acc)[:n]
    return out


# ----------------------------------------------------------------- pipeline real
def medir(c, f_syn, t0_prev):
    t0g, sigf, ssub, subs, status = sw.epoca_por_subconjunto(c["t"], f_syn, c["P"], t0_prev, t14_h=c["t14_h"], n_sub=2, modo="temporada", exigir_cobertura=True)
    return t0g, sigf, ssub, subs, status


def tarefa(chave, rot):
    t_ini = time.time()
    try:
        c = carregar_curva(chave)
        # a verdade: epoca global da curva REAL (o controle) e o P da escada
        t0g_real, sigf_real, ssub_real, subs_real, st_real = medir(c, c["f"], c["t0_prev"])
        if not np.isfinite(t0g_real):
            # 237116051: a epoca GLOBAL da curva inteira encosta na borda (NaN) e as temporadas medem
            # normalmente - a cadeia so usa as temporadas. A verdade passa a ser a media das epocas
            # de temporada levadas ao ciclo mais proximo de t0_prev.
            assert len(subs_real) >= 1, (chave, st_real)
            ts = np.array([s["t0"] for s in subs_real])
            T0 = float(np.mean(ts - np.round((ts - c["t0_prev"]) / c["P"]) * c["P"]))
        else:
            T0 = float(t0g_real)
        modelo = forma(c, rot, T0)
        r = residuos_reais(c["t"], c["f"], c["P"], t0_por_ponto(c["t"], T0, subs_real))
        blocos = noites(c["t"])
        rng = np.random.default_rng(SEMENTE + zlib.crc32(f"{chave}|{rot}".encode()) % 100000)
        linhas = []
        for i in range(N_REAL):
            f_syn = modelo + bootstrap_blocos(r, blocos, rng)
            t0g, sigf, ssub, subs, status = medir(c, f_syn, T0)
            k = np.round((t0g - T0) / c["P"]) if np.isfinite(t0g) else np.nan
            vt = [float((s["t0"] - (T0 + np.round((s["t0"] - T0) / c["P"]) * c["P"])) * 1440.0) for s in subs]
            linhas.append({"curva": str(chave), "tipo": c["tipo"], "forma": rot, "real": i,
                           "vies_global_min": float((t0g - (T0 + k * c["P"])) * 1440.0) if np.isfinite(t0g) else np.nan,
                           "sig_formal_min": float(sigf), "n_temporadas": len(subs),
                           "disp_temporadas_min": float(ssub * np.sqrt(len(subs))) if np.isfinite(ssub) and len(subs) >= 2 else np.nan,
                           "status": status, "prof_ppm": float(np.mean([s["prof"] for s in subs]) * 1e6) if subs else np.nan,
                           "vies_temporadas_min": vt, "t_temporadas": [float(s["t_medio"]) for s in subs]})
        R = pd.DataFrame(linhas)
        R.to_parquet(SAIDA / f"vies__{str(chave).replace(' ', '_')}__{rot}.parquet", index=False)
        info = {"chave": str(chave), "forma": rot, "status": "ok", "n_pontos": int(len(c["t"])), "n_noites": len(blocos), "T0_verdade": T0,
                "real": {"t0g": T0, "sig_formal_min": float(sigf_real), "disp_min": float(ssub_real * np.sqrt(len(subs_real))) if np.isfinite(ssub_real) else np.nan,
                         "n_temporadas": len(subs_real), "status": st_real, "desvios_min": [float(s["desvio_min"]) for s in subs_real],
                         "t0_temporadas": [float(s["t0"]) for s in subs_real]},
                "segundos": time.time() - t_ini}
    except Exception as e:
        info = {"chave": str(chave), "forma": rot, "status": f"erro: {e}", "traceback": traceback.format_exc(), "segundos": time.time() - t_ini}
    (SAIDA / f"info__{str(chave).replace(' ', '_')}__{rot}.json").write_text(json.dumps(info, indent=1, default=float), encoding="utf-8")
    return info


def tarefas(tics):
    return [(h, "transito") for h in HJ] + [(t, rot) for t in tics for rot in ("destacado", "contato")]


def medir_todas(tics, n_proc=N_PROC):
    SAIDA.mkdir(parents=True, exist_ok=True)
    pend = [(k, r) for k, r in tarefas(tics) if not (SAIDA / f"info__{str(k).replace(' ', '_')}__{r}.json").exists()
            or json.loads((SAIDA / f"info__{str(k).replace(' ', '_')}__{r}.json").read_text(encoding="utf-8")).get("status") != "ok"]
    print(f"medir: {len(pend)} tarefas pendentes de {len(tarefas(tics))}, {n_proc} processos, {N_REAL} realizacoes cada", flush=True)
    with ProcessPoolExecutor(max_workers=n_proc) as ex:
        futs = [ex.submit(tarefa, k, r) for k, r in pend]
        for fu in as_completed(futs):
            i = fu.result()
            print(f"  {i['chave']} {i['forma']}: {i['status']} | {i.get('n_pontos')} pts, {i.get('n_noites')} noites | {i['segundos'] / 60:.1f} min", flush=True)


# ----------------------------------------------------------------- controle
def controle(tics):
    """A curva real pelo mesmo caminho: 26 alvos contra os `pontos`; 4 HJ contra hj_calibracao (desvios, barras e O-C)."""
    linhas = []
    for tic in tics:
        c = curva_alvo(tic)
        t0g, sigf, ssub, subs, st = medir(c, c["f"], c["t0_prev"])
        meus = np.array(sorted(s["t0"] for s in subs)); ref = np.array(sorted(c["ref_temporadas_bjd"]))
        dif = np.abs(meus - ref).max() * 1440.0 if len(meus) == len(ref) else np.nan
        linhas.append({"curva": f"TIC {tic}", "n_temporadas": len(subs), "n_ref": len(ref), "max_dif_temporadas_min": float(dif), "status": st,
                       "disp_min": float(ssub * np.sqrt(len(subs))) if np.isfinite(ssub) else np.nan})
        print(f"  TIC {tic}: {len(subs)} temporadas (ref {len(ref)}), max |dif| {dif:.4f} min, disp {linhas[-1]['disp_min']:.2f} min | {st}", flush=True)
    for nome in HJ:
        c = curva_hj(nome)
        t0g, sigf, ssub, subs, st = medir(c, c["f"], c["t0_prev"])
        meus = sorted(round(s["desvio_min"], 2) for s in subs); ref = sorted(c["ref_desvios"])
        dif = max(abs(a - b) for a, b in zip(meus, ref)) if len(meus) == len(ref) else np.nan
        disp = float(ssub * np.sqrt(len(subs))) if np.isfinite(ssub) else np.nan
        rec = {"curva": nome, "n_temporadas": len(subs), "n_ref": len(ref), "max_dif_desvios_min": float(dif), "dif_formal_min": float(abs(sigf - c["ref_formal_min"])),
               "dif_disp_min": float(abs(disp - c["ref_disp_min"])), "status": st, "oc_ref_min": c["ref_oc_min"]}
        # O-C contra a escada TESS refeita como em run_hj_calib (cache MAST; consulta de produtos na rede)
        try:
            import run_hj_calib as RHC
            cur = RHC.curvas_spoc(nome)
            vistos, unicos = set(), []
            for s_, t_, f_ in cur:
                if s_ not in vistos:
                    vistos.add(s_); unicos.append((s_, t_, f_))
            rows = []
            for s, t, f in sorted(unicos, key=lambda x: x[0]):
                t0, sg, pf = sw.medir_epoca(t, f, c["P"], sw.semente_local(t, f, c["P"]), t14_h=c["t14_h"])
                if np.isfinite(t0) and np.isfinite(sg) and 0 < sg < 8:
                    rows.append({"setor": s, "t0": t0, "sig": sg, "prof": pf})
            esc = sw.escada_de_periodo(pd.DataFrame(rows)[["t0", "sig"]], c["P"])
            P_e, T0_e, sP_e, _, _ = esc
            kw = np.round((t0g - (T0_e + 2457000)) / P_e)
            rec["oc_min"] = float((t0g - (T0_e + 2457000 + kw * P_e)) * 1440)
            rec["dif_oc_min"] = abs(rec["oc_min"] - c["ref_oc_min"])
            rec["P_escada_refeita"] = float(P_e); rec["P_ref"] = c["P"]
        except Exception as e:
            rec["oc_min"] = np.nan; rec["dif_oc_min"] = np.nan; rec["erro_escada"] = str(e)[:200]
        linhas.append(rec)
        print(f"  {nome}: desvios {meus} (ref {ref}) max |dif| {dif:.3f}; formal {sigf:.3f} (ref {c['ref_formal_min']:.3f}); disp {disp:.3f} (ref {c['ref_disp_min']:.3f}); "
              f"O-C {rec['oc_min']:+.3f} (ref {c['ref_oc_min']:+.3f}) | {st}", flush=True)
    C = pd.DataFrame(linhas)
    C.to_parquet(BASE / "vies_forma_swasp_controle.parquet", index=False)
    print(f"  -> {BASE / 'vies_forma_swasp_controle.parquet'}")
    return C


# ----------------------------------------------------------------- resumo
def resumo(tics):
    infos = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(SAIDA.glob("info__*.json"))]
    erros = [i for i in infos if i["status"] != "ok"]
    for i in erros:
        print(f"ERRO {i['chave']} {i['forma']}: {i['status']}")
    R = pd.concat([pd.read_parquet(p) for p in sorted(SAIDA.glob("vies__*.parquet"))], ignore_index=True)
    R.to_parquet(BASE / "vies_forma_swasp.parquet", index=False)
    real = {(i["chave"], i["forma"]): i["real"] for i in infos if i["status"] == "ok"}
    linhas = []
    for (cv, rot), g in R.groupby(["curva", "forma"]):
        v = g.vies_global_min.dropna().values
        vt = np.concatenate([np.asarray(x) for x in g.vies_temporadas_min]) if len(g) else np.array([])
        rr = real.get((cv, rot), {})
        linhas.append({"curva": cv, "tipo": g.tipo.iloc[0], "forma": rot, "n_real": int(len(g)), "n_ok": int(len(v)),
                       "vies_medio_min": float(v.mean()), "vies_mediano_min": float(np.median(v)), "vies_sd_min": float(v.std(ddof=1)), "vies_p2.5": float(np.percentile(v, 2.5)), "vies_p97.5": float(np.percentile(v, 97.5)),
                       "ic95_media_min": float(1.96 * v.std(ddof=1) / np.sqrt(len(v))), "vies_temporadas_medio_min": float(vt.mean()) if len(vt) else np.nan,
                       "vies_temporadas_sd_min": float(vt.std(ddof=1)) if len(vt) > 1 else np.nan,
                       "sig_formal_sint_med": float(g.sig_formal_min.median()), "disp_sint_med": float(g.disp_temporadas_min.median()),
                       "disp_real": rr.get("disp_min", np.nan), "sig_formal_real": rr.get("sig_formal_min", np.nan),
                       "n_temporadas_real": rr.get("n_temporadas", np.nan), "n_temporadas_sint_med": float(g.n_temporadas.median()),
                       "frac_nao_testavel": float((~g.status.str.startswith("testado")).mean())})
    S = pd.DataFrame(linhas)
    S["razao_disp_sint_real"] = S.disp_sint_med / S.disp_real
    # estimador instavel (sd entre realizacoes > 10 min: a temporada sintetica cai em minimo errado; 139256217,
    # 359552377, 229476285, 126763885 - alvos de barra real de 9-360 min) fica FORA dos agregados robustos, e dito.
    S["instavel"] = S.vies_sd_min > SD_INSTAVEL_MIN
    S.to_parquet(BASE / "vies_forma_swasp_resumo.parquet", index=False)
    print(S.round(3).to_string(index=False))
    out = {"expectativa_commit": EXPECTATIVA, "n_tarefas_ok": int(len(infos) - len(erros)), "erros": [(i["chave"], i["forma"], i["status"]) for i in erros], "formas": {}}
    for rot, g0 in S.groupby("forma"):
        g = g0[~g0.instavel]
        out["formas"][rot] = {"n_instaveis": int(g0.instavel.sum()), "instaveis": sorted(g0.curva[g0.instavel].tolist()),
                              "vies_mediano_das_curvas_estaveis": float(g.vies_mediano_min.median()), "media_dos_medianos_estaveis": float(g.vies_mediano_min.mean()),
                              "sd_dos_medianos_estaveis": float(g.vies_mediano_min.std(ddof=1)) if len(g) > 1 else np.nan,
                              "n_com_vies_acima_0.5": int((g.vies_mediano_min.abs() > 0.5).sum()), "n_com_vies_acima_1": int((g.vies_mediano_min.abs() > 1.0).sum()),
                              "curvas_vies_acima_1": {r.curva: round(float(r.vies_mediano_min), 2) for r in g[g.vies_mediano_min.abs() > 1.0].itertuples()},"n_curvas": int(len(g)), "vies_medio_das_curvas": float(g.vies_medio_min.mean()), "sd_entre_curvas": float(g.vies_medio_min.std(ddof=1)) if len(g) > 1 else np.nan,
                              "ic95_media": float(1.96 * g.vies_medio_min.std(ddof=1) / np.sqrt(len(g))) if len(g) > 1 else np.nan,
                              "mediana": float(g.vies_medio_min.median()), "max_abs": float(g.vies_medio_min.abs().max()),
                              "sd_realizacoes_mediana": float(g.vies_sd_min.median()), "razao_disp_sint_real_mediana": float(g.razao_disp_sint_real.median())}
        print(f"== {rot}: {len(g0)} curvas, {int(g0.instavel.sum())} instaveis (sd > {SD_INSTAVEL_MIN} min: {sorted(g0.curva[g0.instavel].tolist())}) fora dos agregados; "
              f"medianos por curva (estaveis): media {g.vies_mediano_min.mean():+.3f}, mediana {g.vies_mediano_min.median():+.3f}, sd {g.vies_mediano_min.std(ddof=1):.3f}; "
              f"|vies| > 0,5 min em {int((g.vies_mediano_min.abs() > 0.5).sum())}, > 1 min em {int((g.vies_mediano_min.abs() > 1.0).sum())} de {len(g)}: "
              f"{ {r.curva: round(float(r.vies_mediano_min), 2) for r in g[g.vies_mediano_min.abs() > 1.0].itertuples()} }")
        print(f"   (medias por curva, estaveis): vies medio {g.vies_medio_min.mean():+.3f} min (sd entre curvas {g.vies_medio_min.std(ddof=1):.3f}; mediana {g.vies_medio_min.median():+.3f}; max |.| {g.vies_medio_min.abs().max():.3f}); "
              f"sd entre realizacoes mediana {g.vies_sd_min.median():.2f} min; disp sint/real mediana {g.razao_disp_sint_real.median():.2f}")
    # a media ponderada dos 4 HJ como na cadeia (pesos 1/sig_total^2 de hj_calibracao)
    cal = pd.read_parquet(config.RESULTS / "hj_calibracao.parquet").set_index("pl_name")
    h = S[S.forma == "transito"].set_index("curva")
    w = np.array([1 / cal.loc[n, "sig_total_min"] ** 2 for n in h.index])
    vm = float(np.sum(w * h.vies_medio_min.values) / w.sum())
    # IC pela distribuicao das medias ponderadas por realizacao
    Rh = R[R.forma == "transito"].pivot(index="real", columns="curva", values="vies_global_min")[list(h.index)]
    vr = (Rh.values * w).sum(axis=1) / w.sum()
    vr = vr[np.isfinite(vr)]
    oc_real = float(np.sum(w * np.array([cal.loc[n, "oc_min"] for n in h.index])) / w.sum())
    out["hj_ponderado"] = {"vies_sintetico_min": vm, "p2.5": float(np.percentile(vr, 2.5)), "p97.5": float(np.percentile(vr, 97.5)), "sd_realizacoes": float(vr.std(ddof=1)),
                           "oc_real_ponderado_min": oc_real, "vies_cadeia_min": co.VIES_CADEIA_MIN, "sig_cadeia_min": co.VIES_CADEIA_SIG,
                           "reproduz_2_14": bool(abs(vm - co.VIES_CADEIA_MIN) < co.VIES_CADEIA_SIG)}
    print(f"== HJ (i) ponderado como na cadeia: vies sintetico {vm:+.3f} min (realizacoes 2,5-97,5%: {np.percentile(vr, 2.5):+.2f}..{np.percentile(vr, 97.5):+.2f}); "
          f"O-C real ponderado {oc_real:+.3f}; cadeia {co.VIES_CADEIA_MIN:+.2f} +- {co.VIES_CADEIA_SIG:.2f} -> {'REPRODUZ' if out['hj_ponderado']['reproduz_2_14'] else 'NAO reproduz'} o -2,14")
    # (ii) e (iii) contra (i): diferenca das medias com IC por bootstrap sobre curvas
    rng = np.random.default_rng(1)
    vi = h.vies_mediano_min.values
    for rot in ("destacado", "contato"):
        vj = S[(S.forma == rot) & ~S.instavel].vies_mediano_min.values
        d = vj.mean() - vi.mean()
        boots = [rng.choice(vj, len(vj)).mean() - rng.choice(vi, len(vi)).mean() for _ in range(5000)]
        out["formas"][rot]["dif_vs_transito_min"] = float(d); out["formas"][rot]["dif_ic95"] = [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]
        print(f"== ({rot}) - (transito): {d:+.3f} min, IC95 bootstrap [{np.percentile(boots, 2.5):+.3f}, {np.percentile(boots, 97.5):+.3f}]")
    (BASE / "vies_forma_swasp.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(f"  -> vies_forma_swasp.parquet, _resumo.parquet, .json em {BASE}")


if __name__ == "__main__":
    t26 = pd.read_parquet(BASE / "tabela_26.parquet")
    tics = sorted(int(t) for t in t26.TIC)
    assert len(tics) == 26
    if "--controle" in sys.argv:
        controle(tics)
    if "--medir" in sys.argv:
        i = sys.argv.index("--medir")
        n = int(sys.argv[i + 1]) if len(sys.argv) > i + 1 and sys.argv[i + 1].isdigit() else N_PROC
        medir_todas(tics, n)
    if "--resumo" in sys.argv:
        resumo(tics)
