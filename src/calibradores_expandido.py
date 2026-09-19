# -*- coding: utf-8 -*-
"""Rodada vinte, G: a amostra de calibradores do vies arquival, alargada ate onde o arquivo deixa.

O QUE ESTA EM JOGO. O vies da cadeia, -1,26 +- 0,75 min, e a unica quantidade do artigo medida
fora da amostra, e esta medida em QUATRO objetos. Tirar o mais preciso (WASP-18 b) leva a media a
-0,22 +- 1,04 - isto e, a afirmacao repousa num calibrador. A objecao e obvia e a resposta nao e
retorica: ou existem mais calibradores elegiveis, e eles entram, ou nao existem, e isso tambem e
resultado.

COMO OS DOZE FORAM ESCOLHIDOS (item 1 do pedido, respondido antes de medir). Ha uma tabela em
cache, `data/cache/superwasp/hj_calibradores.parquet`, com 35 planetas em transito que tem fonte
SuperWASP no arquivo CERIT-SC, ordenada por uma estimativa de SNR do transito na propria curva
arquival. `run_hj_calib.py` le essa tabela e faz `.head(N_ALVOS)` com `N_ALVOS = 12`. Ou seja: os
doze sao os doze de maior SNR de trinta e cinco, e o corte em doze e de CONVENIENCIA - custo de
baixar curvas TESS na epoca -, nao criterio fisico nem estatistico. Nao ha script versionado que
construa os 35: a tabela foi montada em sessao, do cruzamento da lista de planetas confirmados
(`catalogs.fetch_confirmed`) com o arquivo SuperWASP. Isto fica escrito no artigo como esta aqui.

O QUE ESTE SCRIPT FAZ. Refaz a elegibilidade em codigo e mede todo mundo que passa:
  (i)  fonte SuperWASP com curva que carrega (a tabela dos 35 ja traz `swasp_id`, `npts`);
  (ii) efemeride publicada que inclui timings de solo - Ivshina & Winn (2022), VizieR
       J/ApJS/259/62 table3, e ExoClock III (Kokori et al. 2023), J/ApJS/265/4 table7;
  (iii) epoca da era SuperWASP medida pela MESMA cadeia (`superwasp.epoca_por_subconjunto`, modo
       temporada, cobertura exigida), com a MESMA guarda de dispersao entre temporadas: barra
       entre temporadas acima de CAP_DISPERSAO_MIN = 5 min reprova o calibrador.
O TESS nao entra: o O-C independente e a epoca arquival contra a efemeride publicada.

EXPECTATIVA (escrita e commitada ANTES de rodar)

(a) ELEGIVEIS. Os 35 sao quase todos WASP, e IW22 (382 sistemas) e ExoClock III (450 planetas)
    cobrem bem essa classe. Espero de 20 a 32 dos 35 com pelo menos uma efemeride publicada.
(b) QUE RENDEM EPOCA. Espero 60 a 90% dos elegiveis: a curva carrega e a cobertura de fase passa.
(c) QUE PASSAM A GUARDA. Nos doze, 8 renderam epoca e 4 passaram - metade. Espero 30 a 60% dos que
    rendem epoca, ou seja, algo entre 6 e 15 calibradores no total, contra os 4 de hoje.
(d) O VIES. Espero a media ponderada entre -2,0 e -0,5 min, com barra MENOR que 0,75, e chi2/dof
    entre 0,5 e 2,5. Espero tambem que a sensibilidade a tirar o calibrador mais preciso caia
    abaixo de 0,5 min - e esse o ganho que justifica a rodada.
(e) SE NAO HOUVER MAIS NENHUM elegivel alem dos quatro, isso e o resultado e fecha a objecao por
    esgotamento do arquivo, nao por escolha.

PORTAO. Se a media ponderada nova cair FORA de -1,26 +- 0,75 (isto e, fora de [-2,01, -0,51]),
PARO e reporto antes de tocar na cadeia: seria um estado (G) e mexeria em todos os coeficientes.

Uso:
    python src/calibradores_expandido.py
"""
import json
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402
import superwasp as sw  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="backslashreplace")
except (AttributeError, ValueError):
    pass

BASE = config.DATA / "orquestra" / "oc_lote"
SAIDA = BASE / "calibradores_expandido.parquet"
SAIDA_JSON = BASE / "calibradores_expandido.json"
CAP_DISPERSAO_MIN = 5.0          # a mesma guarda dos quatro: acima disto o calibrador nao serve
JA_MEDIDOS = ["WASP-18 b", "WASP-4 b", "WASP-19 b", "WASP-6 b"]


def chave(nome):
    """'WASP-018', 'WASP-18b' e 'WASP-18 b' na mesma chave, que e a do SISTEMA.

    A letra do planeta e DESCARTADA de proposito: IW22 indexa por sistema e o ExoClock por
    planeta, e guardar a letra fazia IW22 nunca casar - a primeira passada desta rodada mediu com
    um catalogo so e nao disse nada. Quem chama confere que dois candidatos nao caem na mesma
    chave (nenhum sistema desta lista tem dois transitantes)."""
    n = str(nome).strip().upper().replace(" ", "")
    m = re.match(r"^([A-Z0-9\-]+?)-?0*(\d+)\s*[A-Z]?\s*[A-Z]?$", n)
    if m:
        return f"{m.group(1)}-{int(m.group(2))}"
    return n


def efemerides_publicadas():
    """IW22 e ExoClock III inteiras, do VizieR, indexadas pela chave normalizada."""
    from astroquery.vizier import Vizier
    V = Vizier(columns=["**"], row_limit=-1)
    iw = V.get_catalogs("J/ApJS/259/62")["J/ApJS/259/62/table3"].to_pandas()
    ek = V.get_catalogs("J/ApJS/265/4")["J/ApJS/265/4/table7"].to_pandas()
    print(f"  IW22 table3: {len(iw)} sistemas | ExoClock III table7: {len(ek)} planetas")
    out = {}
    for r in iw.itertuples():
        k = chave(r.Sys)
        out.setdefault(k, {})["IW22"] = {"T0": float(r.T0), "e_T0": float(r.e_T0), "P": float(r.Per),
                                         "e_P": float(r.e_Per), "dPdE": 0.0, "e_dPdE": 0.0}
    for r in ek.itertuples():
        k = chave(r.Planet)
        out.setdefault(k, {})["ExoClock3"] = {"T0": float(r.T0), "e_T0": float(r.e_T0), "P": float(r.Per),
                                              "e_P": float(r.e_Per), "dPdE": 0.0, "e_dPdE": 0.0}
    return out


def predizer(ef, t):
    E = np.round((t - ef["T0"]) / ef["P"])
    tp = ef["T0"] + ef["P"] * E + 0.5 * ef["dPdE"] * E ** 2
    var = ef["e_T0"] ** 2 + (E * ef["e_P"]) ** 2
    return int(E), float(tp), float(np.sqrt(var))


def media_ponderada(d, coluna="oc_min", barra="sig_total_min"):
    w = 1.0 / d[barra].values ** 2
    m = float((w * d[coluna].values).sum() / w.sum())
    s = float(1.0 / np.sqrt(w.sum()))
    chi2 = float((w * (d[coluna].values - m) ** 2).sum())
    return {"media_min": m, "sigma_min": s, "chi2": chi2, "dof": int(len(d) - 1),
            "chi2_dof": chi2 / max(len(d) - 1, 1), "n": int(len(d))}


def main():
    cal = pd.read_parquet(config.DATA / "cache" / "superwasp" / "hj_calibradores.parquet")
    print(f"== candidatos com fonte SuperWASP: {len(cal)} (a tabela que `run_hj_calib` corta em 12 por SNR)")
    pub = efemerides_publicadas()
    cal["chave"] = [chave(n) for n in cal.pl_name]
    dup = cal.chave[cal.chave.duplicated()].tolist()
    assert not dup, f"dois candidatos no mesmo sistema, a chave nao serve: {dup}"
    cal["tem_efemeride"] = [k in pub for k in cal.chave]
    eleg = cal[cal.tem_efemeride].reset_index(drop=True)
    print(f"== com efemeride publicada (IW22 e/ou ExoClock III): {len(eleg)} de {len(cal)}")
    print("   sem efemeride:", ", ".join(sorted(cal.pl_name[~cal.tem_efemeride])) or "nenhum")

    linhas = []
    for c in eleg.itertuples():
        efs = pub[c.chave]
        rot0 = "ExoClock3" if "ExoClock3" in efs else "IW22"
        try:
            t, f, _ = sw.carregar(c.swasp_id, float(c.ra), float(c.dec))
        except Exception as e:  # noqa: BLE001
            print(f"  {c.pl_name:14s} curva nao carrega: {type(e).__name__} {str(e)[:70]}")
            linhas.append({"pl_name": c.pl_name, "status": "sem curva", "motivo": str(e)[:120]})
            continue
        ef0 = efs[rot0]
        t0_prev = ef0["T0"] + np.round((np.median(t) - ef0["T0"]) / ef0["P"]) * ef0["P"]
        try:
            t0g, sigf, ssub, subs, st = sw.epoca_por_subconjunto(
                t, f, ef0["P"], t0_prev, t14_h=float(c.dur_h), n_sub=2, modo="temporada", exigir_cobertura=True)
        except Exception as e:  # noqa: BLE001
            print(f"  {c.pl_name:14s} epoca falhou: {type(e).__name__} {str(e)[:70]}")
            linhas.append({"pl_name": c.pl_name, "status": "sem epoca", "motivo": str(e)[:120]})
            continue
        if not np.isfinite(t0g) or not str(st).startswith("testado"):
            print(f"  {c.pl_name:14s} sem epoca utilizavel ({st})")
            linhas.append({"pl_name": c.pl_name, "status": "sem epoca", "motivo": str(st)[:120]})
            continue
        disp = float(ssub * np.sqrt(len(subs))) if np.isfinite(ssub) else np.nan
        sig_ep = max(float(sigf), float(ssub)) if np.isfinite(ssub) else float(sigf)
        passa = np.isfinite(disp) and disp <= CAP_DISPERSAO_MIN
        for rot, ef in efs.items():
            E, tp, sp = predizer(ef, t0g)
            oc = float((t0g - tp) * 1440.0)
            linhas.append({"pl_name": c.pl_name, "status": "medido" if passa else "reprovado na dispersao",
                           "efemeride": rot, "t0g_sw": float(t0g), "E": E, "oc_min": oc,
                           "sig_pred_min": sp * 1440.0, "sig_epoca_min": sig_ep, "sig_formal_min": float(sigf),
                           "disp_temporadas_min": disp, "sig_total_min": float(np.hypot(sig_ep, sp * 1440.0)),
                           "n_temporadas": len(subs), "ja_medido": c.pl_name in JA_MEDIDOS, "snr": float(c.snr)})
        print(f"  {c.pl_name:14s} {len(subs)} temporadas | formal {sigf:5.2f} | disp {disp:8.2f} | "
              f"{'PASSA' if passa else 'reprova'} | " +
              " | ".join(f"{r}: O-C {l['oc_min']:+7.2f} ± {l['sig_total_min']:5.2f}"
                         for r, l in zip(efs, linhas[-len(efs):])))

    R = pd.DataFrame(linhas)
    R.to_parquet(SAIDA, index=False)
    med = R[(R.status == "medido")]
    print(f"\n== rendem epoca: {R[R.status.isin(['medido', 'reprovado na dispersao'])].pl_name.nunique()} de {len(eleg)}")
    print(f"== passam a guarda de dispersao ({CAP_DISPERSAO_MIN:.0f} min): {med.pl_name.nunique()}")
    res = {"n_candidatos": int(len(cal)), "n_elegiveis": int(len(eleg)),
           "n_rendem_epoca": int(R[R.status.isin(["medido", "reprovado na dispersao"])].pl_name.nunique()),
           "n_passam": int(med.pl_name.nunique()), "cap_dispersao_min": CAP_DISPERSAO_MIN,
           "sem_efemeride": sorted(cal.pl_name[~cal.tem_efemeride]), "conjuntos": {}}
    for rot in ("IW22", "ExoClock3"):
        d = med[med.efemeride == rot]
        if len(d) < 2:
            continue
        m = media_ponderada(d)
        # sensibilidade: tirar o calibrador de menor barra total
        i = d.sig_total_min.idxmin()
        m["sem_o_mais_preciso"] = media_ponderada(d.drop(index=i))
        m["mais_preciso"] = str(d.loc[i, "pl_name"])
        m["alvos"] = sorted(d.pl_name)
        res["conjuntos"][rot] = m
        print(f"\n== {rot}: media ponderada {m['media_min']:+.2f} ± {m['sigma_min']:.2f} min "
              f"(n={m['n']}, χ² {m['chi2']:.1f}/{m['dof']} = {m['chi2_dof']:.2f})")
        s = m["sem_o_mais_preciso"]
        print(f"   sem {m['mais_preciso']} (o mais preciso): {s['media_min']:+.2f} ± {s['sigma_min']:.2f} "
              f"(desloca {abs(s['media_min'] - m['media_min']):.2f} min)")
    SAIDA_JSON.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    print(f"\n-> {SAIDA} e {SAIDA_JSON}")


if __name__ == "__main__":
    main()
