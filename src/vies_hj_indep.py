# -*- coding: utf-8 -*-
"""Rodada seis, Passo E: a calibracao dos 4 Jupiteres quentes contra efemerides
PUBLICADAS que incluem dados de solo (independentes da escada TESS da cadeia).

A cadeia (run_hj_calib) mede a epoca SuperWASP global de cada HJ e a compara
com a efemeride LINEAR ajustada so aos setores TESS (2018-2026) propagada
12-13 anos para tras: O-C -2,85 / -2,55 / -0,64 / +0,85 min (WASP-18/4/19/6),
media ponderada -2,14 +- 0,78 = VIES_CADEIA. Um decaimento orbital real
(WASP-4 b) ou qualquer erro de periodo entra nessa propagacao. Aqui as MESMAS
epocas SuperWASP (curvas do cache, `superwasp.carregar` +
`epoca_por_subconjunto`, como o controle do Passo 4) sao comparadas com:
  (a) Ivshina & Winn 2022 (ApJS 259, 62; VizieR J/ApJS/259/62 table3):
      efemeride linear de todos os transitos publicados 1991-2021 (Ndata
      122 / 25 / 102 / 141), T0 no epoch de referencia que minimiza a
      covariancia (nao publicada; tratada como diagonal, dito assim);
  (b) ExoClock III (Kokori et al. 2023, ApJS 265, 4; J/ApJS/265/4 table7):
      efemeride linear revisada (solo + TESS), idem;
  (c) WASP-4 b: Bouma et al. 2020 (ApJL 893, L29; arXiv:2004.00637) Tabela 3,
      modelo com derivada de periodo: t0 = 2456180,558872(31) BJD_TDB, P =
      1,338231502(24) d, dP/dt = -2,74(40) x 10^-10 (Eq. 3 do texto: -8,64 +-
      1,26 ms/ano; o texto cita +-0,28 apos reescalar por sqrt(chi2_red) e a
      tabela +-0,40 - usa-se a tabela, o valor maior); t(E) = t0 + P E + 1/2
      (dP/dE) E^2 com dP/dE = P dP/dt; covariancia nao publicada ("small",
      E = 0 na media ponderada), tratada como diagonal.
Por HJ e efemeride: E = round((t_SW - T0)/P), O-C = t_SW - t(E), barra da
predicao propagada (sigma_T0^2 + E^2 sigma_P^2 [+ (E^2/2)^2 sigma_dPdE^2]),
barra da epoca = max(formal, std/sqrt n das temporadas) como na cadeia (run_hj_calib),
total = hypot. Media ponderada dos 4, chi2 em torno da media, comparacao
com -2,14 +- 0,78 e com 0. Tambem por temporada. Saidas:
vies_hj_indep.parquet, vies_hj_indep.json. Se uma efemeride nao estiver
acessivel, diz-se; nada e substituido.

EXPECTATIVA (escrita e commitada ANTES de rodar):
  - As epocas SuperWASP reproduzem o controle do Passo 4 (t0g a < 0,01 min).
  - IW22 e ExoClock III concordam entre si em 2006-08 a < 0,5 min nos 4;
    barra da predicao < 0,3 min nos 4 (Ndata grandes, epocas de referencia
    a 2-8 anos das temporadas).
  - WASP-4 b: o O-C contra a quadratica de Bouma 2020 fica +2 a +3,5 min
    ACIMA do da cadeia (-2,55): a efemeride TESS-only de 2019-26 propagada a
    2007 carrega o decaimento (~-3 min); contra as lineares de IW22/ExoClock
    (que absorvem o decaimento em P) fica entre os dois.
  - WASP-18 b permanece em -2,5 +- 1,1 (efemeride mais apertada; nenhum
    decaimento publicado ao nivel de 1 min em 2007).
  - Media ponderada dos 4 (barras da cadeia): **-1,6 +- 0,8** (faixa -2,2 a
    -1,0): compativel com -2,14 a 1 sigma e com 0 a ~2 sigma; chi2 dos 4 em
    torno da media 2-6 (3 gl). Leitura esperada: parte do -2,14 e o
    decaimento do WASP-4 b e o resto e WASP-18 b sozinho a 2,5 sigma - a
    correcao como "vies da cadeia" nao se sustenta com 4 calibradores.
Desvio em qualquer direcao vai ao relatorio.

RESULTADO (2026-09-17, contra a expectativa de 5520329): epocas SuperWASP
reproduzem o controle do Passo 4 a 0,003-0,02 min (WASP-19 0,067: temporada
ruidosa, alinhamento da grade). IW22 e ExoClock III concordam a 0,03-0,40
min; predicao 0,06-0,29 min. O-C: WASP-18 b -2,38/-2,41 +- 1,05 (cadeia
-2,85); WASP-4 b -1,28/-1,45 lineares e -0,81 +- 1,45 com Bouma 2020
(cadeia -2,55: +1,74 min, esperado +2 a +3,5 - um pouco menos); WASP-19 b
+0,46/+0,20 +- 2,59; WASP-6 b +0,40/0,00 +- 1,85. Media ponderada dos 4:
-1,26 (IW22 + Bouma) a -1,56 (ExoClock linear) +- 0,75; com Bouma no WASP-4:
-1,26 +- 0,75 (IW22) e -1,39 +- 0,74 (ExoClock): 1,7-1,9 sigma de 0 e 0,7-0,8
sigma de -2,14; chi2 2,0-2,4 / 3. Sem WASP-18 b: -0,2 a -0,4 +- 1,0. Uma
primeira rodada usou std em vez de std/sqrt n na barra global (WASP-19
4,48 em vez de 2,59) e deu -1,36 a -1,66: corrigido para a barra da cadeia
antes de gravar. Leitura: o -2,14 tinha 0,5-0,9 min de propagacao TESS-only
(decaimento do WASP-4 b e o P do WASP-18 b); o que sobra, -1,3 a -1,6 +-
0,75, e um efeito a ~2 sigma sustentado por um calibrador.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cadeia_oc as co  # noqa: E402
import config  # noqa: E402
import superwasp as sw  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
HJ = (("WASP-18 b", "WASP-018", "WASP-18b"), ("WASP-4 b", "WASP-004", "WASP-4b"), ("WASP-19 b", "WASP-019", "WASP-19b"), ("WASP-6 b", "WASP-006", "WASP-6b"))
BOUMA2020 = {"t0": 2456180.558872, "e_t0": 31e-6, "P": 1.338231502, "e_P": 24e-9, "dPdt": -2.74e-10, "e_dPdt": 0.40e-10,
             "fonte": "Bouma et al. 2020, ApJL 893, L29, Tabela 3 (modelo com derivada de periodo); covariancia nao publicada"}
EXPECTATIVA = "5520329"


def efemerides(nome_iw, nome_ek):
    iw = pd.read_parquet(config.CATALOGS / "ivshina_winn_2022_hj4.parquet").set_index("Sys").loc[nome_iw]
    ek = pd.read_parquet(config.CATALOGS / "exoclock3_hj4.parquet").set_index("Planet").loc[nome_ek]
    return {"IW22": {"T0": float(iw.T0), "e_T0": float(iw.e_T0), "P": float(iw.Per), "e_P": float(iw.e_Per), "dPdE": 0.0, "e_dPdE": 0.0, "Ndata": int(iw.Ndata), "n_T0": str(iw.n_T0)},
            "ExoClock3": {"T0": float(ek.T0), "e_T0": float(ek.e_T0), "P": float(ek.Per), "e_P": float(ek.e_Per), "dPdE": 0.0, "e_dPdE": 0.0}}


def predizer(ef, t):
    E = np.round((t - ef["T0"]) / ef["P"])
    tp = ef["T0"] + ef["P"] * E + 0.5 * ef["dPdE"] * E ** 2
    var = ef["e_T0"] ** 2 + (E * ef["e_P"]) ** 2 + (0.5 * E ** 2 * ef["e_dPdE"]) ** 2
    return E, tp, np.sqrt(var)


if __name__ == "__main__":
    cal = pd.read_parquet(config.DATA / "cache" / "superwasp" / "hj_calibradores.parquet").set_index("pl_name")
    res = pd.read_parquet(config.RESULTS / "hj_calibracao.parquet").set_index("pl_name")
    linhas, temporadas = [], []
    for nome, n_iw, n_ek in HJ:
        c = cal.loc[nome]; r = res.loc[nome]
        efs = efemerides(n_iw, n_ek)
        if nome == "WASP-4 b":
            b = BOUMA2020
            efs["Bouma2020_quad"] = {"T0": b["t0"], "e_T0": b["e_t0"], "P": b["P"], "e_P": b["e_P"], "dPdE": b["P"] * b["dPdt"], "e_dPdE": b["P"] * b["e_dPdt"]}
        t, f, _ = sw.carregar(c.swasp_id, float(c.ra), float(c.dec))
        ef0 = efs["ExoClock3"]
        t0_prev = ef0["T0"] + np.round((np.median(t) - ef0["T0"]) / ef0["P"]) * ef0["P"]
        t0g, sigf, ssub, subs, st = sw.epoca_por_subconjunto(t, f, ef0["P"], t0_prev, t14_h=float(c.dur_h))
        disp = float(ssub * np.sqrt(len(subs))) if np.isfinite(ssub) else np.nan
        # a barra da epoca GLOBAL e a da cadeia (run_hj_calib): max(formal, std/sqrt n) - o erro da media das
        # temporadas, que e a quantidade certa para a epoca combinada (nao std, que e a barra de UMA temporada)
        sig_ep = max(float(sigf), float(ssub)) if np.isfinite(ssub) else float(sigf)
        ctrl = json.loads((BASE / "vies_forma_26" / f"info__{nome.replace(' ', '_')}__transito.json").read_text(encoding="utf-8"))["real"]
        dif_ctrl = abs(t0g - ctrl["t0g"]) * 1440.0
        for rot, ef in efs.items():
            E, tp, sp = predizer(ef, np.array([t0g]))
            oc = float((t0g - tp[0]) * 1440.0)
            linhas.append({"hj": nome, "efemeride": rot, "t0g_sw": float(t0g), "E": int(E[0]), "oc_min": oc, "sig_pred_min": float(sp[0] * 1440.0),
                           "sig_epoca_min": sig_ep, "sig_formal_min": float(sigf), "disp_temporadas_min": disp, "sig_total_min": float(np.hypot(sig_ep, sp[0] * 1440.0)),
                           "oc_cadeia_min": float(r.oc_min), "sig_cadeia_min": float(r.sig_total_min), "n_temporadas": len(subs), "dif_controle_min": dif_ctrl, "status": st})
            for s in subs:
                Es, tps, sps = predizer(ef, np.array([s["t0"]]))
                temporadas.append({"hj": nome, "efemeride": rot, "t_medio": s["t_medio"], "t0": s["t0"], "E": int(Es[0]), "oc_min": float((s["t0"] - tps[0]) * 1440.0),
                                   "sig_formal_min": s["sigma_formal_min"], "sig_pred_min": float(sps[0] * 1440.0)})
        print(f"{nome}: t0g {t0g:.5f} (controle dif {dif_ctrl:.4f} min), formal {sigf:.2f}, disp {disp:.2f}, {len(subs)} temporadas | " +
              " | ".join(f"{rot}: O-C {l['oc_min']:+.2f} +- {l['sig_total_min']:.2f} (pred {l['sig_pred_min']:.2f})" for rot, l in zip(efs, linhas[-len(efs):]))
              + f" | cadeia {r.oc_min:+.2f} +- {r.sig_total_min:.2f}", flush=True)
    R = pd.DataFrame(linhas); T = pd.DataFrame(temporadas)
    R.to_parquet(BASE / "vies_hj_indep.parquet", index=False); T.to_parquet(BASE / "vies_hj_indep_temporadas.parquet", index=False)
    out = {"expectativa_commit": EXPECTATIVA, "bouma2020": BOUMA2020, "conjuntos": {}}
    print("\n== media ponderada dos 4 (pesos 1/sig_total^2, barras da epoca como na cadeia):")
    conjuntos = {"IW22 (linear, os 4)": {n: "IW22" for n, _, _ in HJ}, "ExoClock3 (linear, os 4)": {n: "ExoClock3" for n, _, _ in HJ},
                 "IW22 + Bouma2020 quad no WASP-4": {n: ("Bouma2020_quad" if n == "WASP-4 b" else "IW22") for n, _, _ in HJ},
                 "ExoClock3 + Bouma2020 quad no WASP-4": {n: ("Bouma2020_quad" if n == "WASP-4 b" else "ExoClock3") for n, _, _ in HJ},
                 "cadeia (TESS-only, run_hj_calib)": None}
    for rot, sel in conjuntos.items():
        if sel is None:
            oc = np.array([res.loc[n, "oc_min"] for n, _, _ in HJ]); s = np.array([res.loc[n, "sig_total_min"] for n, _, _ in HJ])
        else:
            rows = [R[(R.hj == n) & (R.efemeride == sel[n])].iloc[0] for n, _, _ in HJ]
            oc = np.array([x.oc_min for x in rows]); s = np.array([x.sig_total_min for x in rows])
        w = 1 / s ** 2; m = float(np.sum(w * oc) / w.sum()); sm = float(1 / np.sqrt(w.sum()))
        chi2 = float(np.sum(((oc - m) / s) ** 2)); p = float(stats.chi2.sf(chi2, len(oc) - 1))
        out["conjuntos"][rot] = {"oc": {n: float(v) for (n, _, _), v in zip(HJ, oc)}, "sig": {n: float(v) for (n, _, _), v in zip(HJ, s)}, "media_min": m, "sigma_min": sm,
                                 "chi2": chi2, "dof": len(oc) - 1, "p": p, "z_vs_zero": m / sm, "z_vs_cadeia": (m - co.VIES_CADEIA_MIN) / np.hypot(sm, co.VIES_CADEIA_SIG),
                                 "sem_WASP4_media": float(np.sum((w * oc)[[0, 2, 3]]) / w[[0, 2, 3]].sum()), "sem_WASP4_sigma": float(1 / np.sqrt(w[[0, 2, 3]].sum())),
                                 "sem_WASP18_media": float(np.sum((w * oc)[1:]) / w[1:].sum()), "sem_WASP18_sigma": float(1 / np.sqrt(w[1:].sum()))}
        o = out["conjuntos"][rot]
        print(f"  {rot:40s}: O-C {', '.join(f'{v:+.2f}' for v in oc)} | media {m:+.2f} +- {sm:.2f} min | chi2 {chi2:.1f}/{len(oc) - 1} (p {p:.2f}) | "
              f"vs 0: {m / sm:+.1f} sigma | vs -2,14 +- 0,78: {o['z_vs_cadeia']:+.1f} sigma | sem WASP-4: {o['sem_WASP4_media']:+.2f} +- {o['sem_WASP4_sigma']:.2f}; sem WASP-18: {o['sem_WASP18_media']:+.2f} +- {o['sem_WASP18_sigma']:.2f}")
    # concordancia entre as efemerides publicadas nas epocas SuperWASP
    conc = {}
    for n, _, _ in HJ:
        a = R[(R.hj == n) & (R.efemeride == "IW22")].iloc[0]; b = R[(R.hj == n) & (R.efemeride == "ExoClock3")].iloc[0]
        conc[n] = {"IW22_menos_ExoClock3_min": float(a.oc_min - b.oc_min) * -1.0 * -1.0, "sig_pred_IW22": float(a.sig_pred_min), "sig_pred_EK3": float(b.sig_pred_min)}
    out["concordancia_publicadas"] = conc
    print("  IW22 - ExoClock3 na epoca SuperWASP (min): " + ", ".join(f"{n} {v['IW22_menos_ExoClock3_min']:+.2f} (pred {v['sig_pred_IW22']:.2f}/{v['sig_pred_EK3']:.2f})" for n, v in conc.items()))
    print("  por temporada (O-C min):"); print(T.pivot_table(index=["hj", "t_medio"], columns="efemeride", values="oc_min").round(2).to_string())
    (BASE / "vies_hj_indep.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(f"  -> {BASE / 'vies_hj_indep.parquet'}, _temporadas.parquet, .json")
