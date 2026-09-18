# -*- coding: utf-8 -*-
"""Estado B -> estado D: o que a barra TESS |a - b| / 2 mudou nos registros dos 88.

Le os registros do estado B direto do git (commit 12b11bd, o ultimo antes da
troca) e os do estado D em data/orquestra/oc_lote/, e grava
`comparacao_estado_B_D.parquet` (uma linha por alvo medido nos dois estados)
com dP/dt, sigma, chi2_red, p_curv, p_adv, p_gof e classe nos dois estados,
o fator de barra TESS por epoca (mediana e minimo por alvo), o deslocamento
maximo das epocas SuperWASP (min) e a diferenca contra o refit
`barra_tess_metades_26.parquet` (estado B refeito com a regra /2 sem a
cadeia), em unidades de sigma.

Tambem confere o que NAO pode mudar: o conjunto dos 26 medidos, a Tabela 2
(o primeiro guard de cada um dos 62), o E de cada epoca, e as barras
SuperWASP (nao dependem da barra TESS a nao ser pelo P da escada).
"""
import glob
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
COMMIT_B = "12b11bd"   # o ultimo commit do ESTADO B (nota rev34): barra TESS |a - b| / sqrt 2; o estado D comeca em 98a3d69


def registro_B(tic):
    raw = subprocess.run(["git", "show", f"{COMMIT_B}:data/orquestra/oc_lote/oc__{tic}.json"],
                         capture_output=True, check=True, cwd=config.DATA.parent).stdout
    return json.loads(raw.decode("utf-8"))


def registro_D(tic):
    return json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))


def erro_curto(j):
    e = str(j.get("erro", ""))
    for k in ("pontos validos", "cobertura por bloco", "veio vazio", "escada nao fechou", "T14 do catalogo inaplicavel"):
        if k in e:
            return k
    return e[:40]


if __name__ == "__main__":
    tics = sorted(int(Path(f).stem.split("__")[1]) for f in glob.glob(str(BASE / "oc__*.json")))
    assert len(tics) == 88, len(tics)
    okB = {t for t in tics if registro_B(t).get("status") == "ok"}
    okD = {t for t in tics if registro_D(t).get("status") == "ok"}
    print(f"medidos: B {len(okB)}, D {len(okD)}; entram {sorted(okD - okB)}, saem {sorted(okB - okD)}")
    assert okB == okD, "o conjunto dos medidos mudou - escalar antes de qualquer outra coisa"
    # Tabela 2: o primeiro guard de cada um dos 62
    dif_guard = [(t, erro_curto(registro_B(t)), erro_curto(registro_D(t))) for t in tics if t not in okB
                 if erro_curto(registro_B(t)) != erro_curto(registro_D(t))]
    print(f"Tabela 2: {len(dif_guard)} alvos mudaram de guard {dif_guard}")

    refit = pd.read_parquet(BASE / "barra_tess_metades_26.parquet").set_index("tic")
    linhas = []
    for t in sorted(okB):
        b, d = registro_B(t), registro_D(t)
        pb = pd.DataFrame(b["pontos"]).sort_values("t0").reset_index(drop=True)
        pdd = pd.DataFrame(d["pontos"]).sort_values("t0").reset_index(drop=True)
        assert len(pb) == len(pdd) and (pb.fonte.values == pdd.fonte.values).all(), t
        assert (pb.E.values == pdd.E.values).all(), (t, "E mudou")
        tess = pb.fonte.str.startswith("TESS").values
        fator = pdd.sig_min.values[tess] / pb.sig_min.values[tess]
        desl_sw = np.abs(pdd.t0.values[~tess] - pb.t0.values[~tess]) * 1440.0
        dbar_sw = np.abs(pdd.sig_min.values[~tess] - pb.sig_min.values[~tess])
        pa, pdp = b.get("parabola", {}), d.get("parabola", {})
        r = {"tic": t, "n_tess": int(tess.sum()), "fator_tess_mediana": float(np.median(fator)), "fator_tess_min": float(fator.min()),
             "n_tess_encolhem": int((fator < 0.999).sum()), "desl_sw_max_min": float(desl_sw.max()), "dbar_sw_max_min": float(dbar_sw.max()),
             "P_escada_B": b["P_escada_d"], "P_escada_D": d["P_escada_d"], "sP_B": b["sP_escada_d"], "sP_D": d["sP_escada_d"],
             "dPdt_B": pa.get("dPdt_s_por_ano"), "dPdt_D": pdp.get("dPdt_s_por_ano"), "s_B": pa.get("sdPdt_s_por_ano"), "s_D": pdp.get("sdPdt_s_por_ano"),
             "chi2red_B": pa["chi2"] / pa["dof"] if pa else np.nan, "chi2red_D": pdp["chi2"] / pdp["dof"] if pdp else np.nan, "dof": d.get("dof"),
             "p_B": b.get("p_curvatura"), "p_D": d.get("p_curvatura"),
             "padv_B": b.get("adversarial", {}).get("p_adversarial"), "padv_D": d.get("adversarial", {}).get("p_adversarial"),
             "classe_B": b.get("classe"), "classe_D": d.get("classe")}
        from scipy import stats
        r["pgof_B"] = float(stats.chi2.sf(pa["chi2"], pa["dof"])) if pa else np.nan
        r["pgof_D"] = float(stats.chi2.sf(pdp["chi2"], pdp["dof"])) if pdp else np.nan
        if t in refit.index:
            r["refit_dPdt"] = float(refit.loc[t, "dPdt"]); r["refit_s"] = float(refit.loc[t, "s"])
            r["dif_refit_sigma"] = abs(r["dPdt_D"] - r["refit_dPdt"]) / r["refit_s"]
            r["refit_p_adv"] = float(refit.loc[t, "p_adv"])
        linhas.append(r)
    R = pd.DataFrame(linhas)
    R.to_parquet(BASE / "comparacao_estado_B_D.parquet", index=False)
    det = lambda c: (R[f"p_{c}"] < 0.05) & (R[f"padv_{c}"] < 0.05)  # noqa: E731
    print(f"\n== barras TESS: {int(R.n_tess_encolhem.sum())} de {int(R.n_tess.sum())} encolhem; fator mediano nas que encolhem "
          f"{np.median(np.concatenate([[r.fator_tess_min] for r in R.itertuples()])):.2f} (min por alvo)")
    print(f"== epocas SuperWASP: deslocamento maximo {R.desl_sw_max_min.max():.3f} min; barra SW muda no maximo {R.dbar_sw_max_min.max():.3f} min")
    print(f"== curvatura + adversarial: B {int(det('B').sum())} -> D {int(det('D').sum())}; entram {sorted(R[det('D') & ~det('B')].tic)}; saem {sorted(R[det('B') & ~det('D')].tic)}")
    print(f"== p_gof < 0,05: B {int((R.pgof_B < 0.05).sum())} {sorted(R[R.pgof_B < 0.05].tic)} -> D {int((R.pgof_D < 0.05).sum())} {sorted(R[R.pgof_D < 0.05].tic)}")
    print(f"== chi2_red mediano: B {R.chi2red_B.median():.3f} -> D {R.chi2red_D.median():.3f}; agregado B "
          f"{(R.chi2red_B * R.dof).sum() / R.dof.sum():.3f} -> D {(R.chi2red_D * R.dof).sum() / R.dof.sum():.3f} ({int(R.dof.sum())} gl)")
    print(f"== dP/dt: mediana |dD - dB| / sB = {(abs(R.dPdt_D - R.dPdt_B) / R.s_B).median():.3f}, max {(abs(R.dPdt_D - R.dPdt_B) / R.s_B).max():.3f}; "
          f"razao s_D / s_B mediana {(R.s_D / R.s_B).median():.3f}, min {(R.s_D / R.s_B).min():.3f}")
    if "dif_refit_sigma" in R:
        print(f"== contra o refit barra_tess_metades_26: |dP/dt_D - refit| / s: mediana {R.dif_refit_sigma.median():.4f}, max {R.dif_refit_sigma.max():.4f}; "
              f"p_adv D vs refit: max |dif| {(R.padv_D - R.refit_p_adv).abs().max():.4f}")
    print(f"== max dP/dt D {R.dPdt_D.max():+.4f}, max |dP/dt| {R.dPdt_D.abs().max():.4f}; sinais +/-: {int((R.dPdt_D > 0).sum())}/{int((R.dPdt_D < 0).sum())}")
    print(R[["tic", "fator_tess_min", "dPdt_B", "dPdt_D", "s_B", "s_D", "p_B", "p_D", "padv_B", "padv_D", "pgof_B", "pgof_D"]].round(4).to_string(index=False))
    print("gravado:", BASE / "comparacao_estado_B_D.parquet")
