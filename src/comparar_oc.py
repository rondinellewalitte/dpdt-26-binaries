# -*- coding: utf-8 -*-
"""Compara a cadeia rodada com a expectativa por componentes. Qualquer campo
fora da tolerancia, em qualquer direcao, e ESCALAR."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config


def comparar(exp, obs):
    d = []

    def num(campo, e, o, tol):
        if o is None or (isinstance(o, float) and o != o):
            d.append((campo, e, "<ausente>")); return
        if abs(e - o) > tol:
            d.append((campo, e, round(o, 5)))

    T, t = exp["tess"], obs["tess"]
    num("tess.n_aglomerados", T["n_aglomerados"], t["n_aglomerados"], 0)
    ag = {a["setor"]: a for a in t["aglomerados"]}
    for s in ("s4", "s31", "s97", "s105"):
        k = int(s[1:])
        num(f"tess.{s}.epoca_min", T[s]["epoca_btjd"] * 1440, ag[k]["t0"] * 1440, T[s]["tol_min"])
        num(f"tess.{s}.sigma_min", T[s]["sigma_min"], ag[k]["sig_min"], T[s]["tol_sigma"])
    for i, s in enumerate((4, 31, 97, 105)):
        num(f"tess.linear.residuo[{s}]", T["linear"]["residuos_min"][i], ag[s]["res_linear_min"], T["linear"]["tol_min"])
        num(f"tess.parabola.residuo[{s}]", T["parabola"]["residuos_min"][i], ag[s]["res_parabola_min"], T["parabola"]["tol_min"])
    num("tess.linear.chi2", T["linear"]["chi2"], t["linear"]["chi2"], T["linear"]["tol_chi2"])
    num("tess.parabola.chi2", T["parabola"]["chi2"], t["parabola"]["chi2"], T["parabola"]["tol_chi2"])
    num("tess.delta_chi2", T["delta_chi2"], t["delta_chi2"], T["tol_delta_chi2"])
    num("tess.Q", T["Q_d_por_ciclo2"], t["parabola"]["Q_d_por_ciclo2"], T["Q_d_por_ciclo2"] * T["tol_Q_rel"])

    S, s = exp["superwasp"], obs["superwasp"]
    if S["sourceid"] != s["sourceid"]:
        d.append(("superwasp.sourceid", S["sourceid"], s["sourceid"]))
    num("superwasp.npts", S["npts"], s["npts"], S["tol_npts"])
    num("superwasp.n_temporadas", S["n_temporadas"], s["n_temporadas"], 0)
    if not s["status_cobertura_por_bloco"].startswith(S["status_cobertura_por_bloco"]):
        d.append(("superwasp.status", S["status_cobertura_por_bloco"], s["status_cobertura_por_bloco"]))
    num("superwasp.epoca_min", S["epoca_bjd_tdb"] * 1440, s["epoca_bjd_tdb"] * 1440, S["tol_min"])
    num("superwasp.E", S["E"], s["E"], 0)
    num("superwasp.barra_interna_min", S["barra_interna_min"], s["barra_interna_min"], S["tol_barra"])
    num("superwasp.barra_total_min", S["barra_total_min"], s["barra_total_min"], S["tol_barra_total"])

    X, x = exp["teste"], obs["teste"]
    num("teste.residuo_no_linear_min", X["residuo_no_linear_min"], x["residuo_no_linear_min"], X["tol_min"])
    num("teste.sigmas_linear", X["sigmas_linear"], x["sigmas_linear"], X["tol_sigmas"])
    num("teste.residuo_na_parabola_min", X["residuo_na_parabola_min"], x["residuo_na_parabola_min"], X["tol_min_par"])
    num("teste.sigmas_parabola", X["sigmas_parabola"], x["sigmas_parabola"], X["tol_sigmas_par"])
    num("teste.parabola_prevista_min", X["parabola_prevista_2007_min"], x["parabola_prevista_min"], X["tol_prev"])
    num("teste.graus_de_liberdade", X["graus_de_liberdade_com_swasp"], x["graus_de_liberdade_com_swasp"], 0)
    return d


if __name__ == "__main__":
    exp = json.load(open(config.DATA / "orquestra" / "expectativa_oc_142874476.json", encoding="utf-8"))
    obs = json.load(open(config.DATA / "orquestra" / "cadeia_oc__142874476.json", encoding="utf-8"))
    d = comparar(exp, obs)
    n = 4 + 8 + 8 + 4 + 8 + 6   # campos comparados
    print(f"campos comparados: {n} | divergencias: {len(d)}")
    for campo, e, o in d:
        print(f"  ESCALAR  {campo:<36} esperado={e!r:<28} observado={o!r}")
    print("\nstatus:", "REPRODUZIU" if not d else "ESCALAR")
    print("adversarial (primeiro registro):", obs["adversarial"])
