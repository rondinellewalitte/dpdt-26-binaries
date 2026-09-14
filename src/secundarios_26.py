# -*- coding: utf-8 -*-
"""O-C do MINIMO SECUNDARIO nos 26, nas curvas TESS, com o mesmo estimador e
a mesma barra do primario. Ponto 1.2 da revisao: movimento apsidal em
sistema excentrico produz quadratica espuria SO no primario (primario e
secundario andam em sentidos opostos); dP/dt secular move os dois juntos.
E o segundo confundidor da familia do V527 Dra, e nunca foi testado.

Como: `oc_lote.epocas_2min` com a semente do catalogo deslocada de P/2 -
nada mais muda (ancora no setor mais recente, trapezio com o T14 do
catalogo, barra = max(formal, |t0a - t0b|/sqrt 2 das metades)). O primario e
re-extraido pelo mesmo caminho e tem de reproduzir o `pontos` do JSON do
lote - controle embutido.

Por alvo e setor: d = (t_sec - t_prim - P/2) em minutos, modulo P, com
barra hypot(sig_prim, sig_sec). Secundario "mensuravel" = profundidade
ajustada > 0 e barra finita; onde nao for, conta-se e diz-se. Por alvo: d
mediano, e a deriva entre o primeiro e o ultimo setor (min/ano) com barra;
"andam juntos" = deriva compativel com zero em 2 sigma.

EXPECTATIVA (antes de rodar): (a) nos 21 com P < 2 d (contato ou quase),
secundario mensuravel na maioria, d compativel com zero em 2 sigma e deriva
compativel com zero; (b) nos 5 com P entre 2 e 4,4 d (126763885, 139256217,
237116051, 329248002, 359552377), d pode ser != 0 (orbita excentrica) e a
deriva e a medida que interessa - espero 0 a 2 com deriva significativa;
(c) alguns secundarios rasos nao mensuraveis - espero <= 5 alvos. Desvio em
qualquer direcao vai para a nota como esta.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402
import oc_lote as L  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


def medir_alvo(tic, P0, t14_h, t0_cat, man, pontos_json):
    prim = L.epocas_2min(tic, P0, t14_h, man, t0_cat)
    sec = L.epocas_2min(tic, P0, t14_h, man, t0_cat + 0.5 * P0)
    # controle embutido: o primario re-extraido reproduz o JSON do lote
    ref = {int(f.split("s")[1]): t for f, t in zip(pontos_json.fonte, pontos_json.t0) if f.startswith("TESS")}
    for r in prim.itertuples():
        if r.setor in ref and abs(r.t0_btjd - ref[r.setor]) * 1440 > 0.01:
            raise RuntimeError(f"TIC {tic} s{r.setor}: primario re-extraido difere do JSON em "
                               f"{(r.t0_btjd - ref[r.setor]) * 1440:.3f} min")
    m = prim.merge(sec, on="setor", suffixes=("_p", "_s"))
    m["mensuravel"] = (m.prof_ajustada_ppm_s > 0) & np.isfinite(m.sigma_min_s) & (m.sigma_min_s > 0)
    d = (m.t0_btjd_s - m.t0_btjd_p - 0.5 * P0) % P0
    d = np.where(d > 0.5 * P0, d - P0, d) * 1440.0
    m["d_min"] = d
    m["sd_min"] = np.hypot(m.sigma_min_p, m.sigma_min_s)
    m["tic"] = tic
    return m


def resumo(tic, P0, m):
    ok = m[m.mensuravel].sort_values("setor")
    out = {"tic": tic, "P_d": P0, "n_setores": len(m), "n_sec_mensuravel": int(len(ok))}
    if len(ok) == 0:
        out["veredito"] = "secundario nao mensuravel"
        return out
    out["d_mediano_min"] = float(ok.d_min.median())
    out["sd_mediano_min"] = float(ok.sd_min.median())
    out["d_sobre_sigma"] = float(ok.d_min.median() / ok.sd_min.median())
    if len(ok) < 2:
        out["veredito"] = "1 setor: deslocamento medido, deriva nao testavel"
        return out
    a, b = ok.iloc[0], ok.iloc[-1]
    anos = (b.t0_btjd_p - a.t0_btjd_p) / 365.25
    if anos < 0.5:
        out["veredito"] = f"setores a {anos * 12:.1f} meses: deriva nao testavel"
        return out
    deriva = (b.d_min - a.d_min) / anos
    s_deriva = np.hypot(a.sd_min, b.sd_min) / anos
    out.update({"alavanca_anos": float(anos), "deriva_min_por_ano": float(deriva),
                "s_deriva_min_por_ano": float(s_deriva), "deriva_sobre_sigma": float(deriva / s_deriva)})
    out["veredito"] = ("andam juntos (deriva compativel com 0 em 2 sigma)" if abs(deriva / s_deriva) <= 2
                       else "deriva do secundario significativa: primario e secundario NAO andam juntos")
    return out


if __name__ == "__main__":
    if not L.controles():
        sys.exit("controles nao reproduziram: nada roda")
    base = config.DATA / "orquestra" / "oc_lote"
    t26 = pd.read_parquet(base / "tabela_26.parquet")
    alvos = pd.read_parquet(config.DATA / "orquestra" / "alvos_oc_88.parquet").set_index("tic")
    man = L.manifestos()
    setores, resumos = [], []
    for tic in sorted(t26.TIC):
        a = alvos.loc[tic]
        j = json.loads((base / f"oc__{tic}.json").read_text(encoding="utf-8"))
        pontos = pd.DataFrame(j["pontos"])
        try:
            m = medir_alvo(int(tic), float(j["P_ref_d"]), float(a.t14_h), float(a.t0_btjd), man, pontos)
        except Exception as e:  # noqa: BLE001 - contado, nao escondido
            resumos.append({"tic": int(tic), "P_d": float(j["P_ref_d"]), "veredito": f"falhou: {str(e)[:120]}"})
            print(f"TIC {tic}: falhou: {str(e)[:120]}", flush=True)
            continue
        setores.append(m)
        r = resumo(int(tic), float(j["P_ref_d"]), m)
        resumos.append(r)
        print(f"TIC {tic} P {r['P_d']:.4f}: {r['n_sec_mensuravel']}/{r['n_setores']} secundarios mensuraveis"
              + (f"; d {r['d_mediano_min']:+.2f} +- {r['sd_mediano_min']:.2f} min" if "d_mediano_min" in r else "")
              + (f"; deriva {r['deriva_min_por_ano']:+.2f} +- {r['s_deriva_min_por_ano']:.2f} min/ano" if "deriva_min_por_ano" in r else "")
              + f" -> {r['veredito']}", flush=True)
    S = pd.concat(setores, ignore_index=True) if setores else pd.DataFrame()
    R = pd.DataFrame(resumos)
    S.to_parquet(base / "secundarios_26_setores.parquet", index=False)
    R.to_parquet(base / "secundarios_26.parquet", index=False)
    print("\n== resumo (expectativa: maioria mensuravel e andando junto; 0-2 com deriva; <= 5 nao mensuraveis)")
    print(R.veredito.value_counts().to_string())
    if "d_sobre_sigma" in R:
        print(f"  |d|/sigma > 2 (secundario fora da fase 0,5): {(R.d_sobre_sigma.abs() > 2).sum()} de {R.d_sobre_sigma.notna().sum()}")
