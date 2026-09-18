# -*- coding: utf-8 -*-
"""Rodada vinte e um: o segundo template propagado NO ESPACO DE EPOCAS QUE A CADEIA USA.

O DEFEITO QUE ISTO CONSERTA. A 3.4 comparou o deslocamento da epoca arquival GLOBAL (todas as
temporadas juntas) com a coluna de offset minimo da Tabela 3. Mas a 3.3 diz - e o codigo faz - que
cada temporada do SuperWASP entra no ajuste como EPOCA PROPRIA. Para os tres alvos de temporada
unica os dois numeros coincidem; para os demais nao, e ali a cadeia enfrenta os deslocamentos por
temporada (mediana 0,64 min, p90 4,91, maximo 10,0), nao o da combinada (0,54 e 3,27). Havia uma
defesa - um offset COMUM e a forma maximamente degenerada com curvatura, e deslocamentos de sinal
variavel sao bem menos degenerados - mas ela e um argumento, e as epocas alternativas ja estao
medidas. Entao em vez de argumentar, reajusta-se.

O QUE FAZ. Para cada um dos 11 do D11: pega os `pontos` do proprio JSON da cadeia (as epocas TESS
aglomeradas e as temporadas arquivais, com barra, vies e piso ja aplicados), desloca SO as
arquivais pelo dt medido da sua temporada (`estimador_alt_d11.parquet`), e reajusta com a mesma
maquinaria - `cadeia_oc.ajustar` para a linear e a parabola, `cadeia_oc.adversarial` para a porta.
Uma deteccao sobrevive se p_curv < 0,05 E p_adv < 0,05, que e a definicao do D11.

Onde a cadeia partiu uma temporada unica em duas metades (tres alvos), as duas metades recebem o
dt daquela temporada: a largura do template e uma propriedade da temporada, nao da metade. Isso e
declarado porque nesses tres o deslocamento fica COMUM as duas epocas arquivais - o caso mais
degenerado com curvatura, e portanto o mais desfavoravel.

EXPECTATIVA (escrita e commitada ANTES de rodar)

(a) SOBREVIVENTES. Espero 9 a 11 dos 11 sobreviverem. O deslocamento por temporada e maior que o
    da combinada, mas tem sinal variavel dentro do mesmo alvo (198408416: +0,33/+1,20/+0,64;
    329246824: +0,06/-10,02/-0,29), e um deslocamento de sinal variavel nao imita curvatura.
    Menos de 9 e desvio: paro e reporto antes de tocar em qualquer texto.

(b) OS DOIS EM RISCO. 329246824, que tem uma temporada com dt = -10,02 min contra um offset comum
    limiar de 7 min, e 232634196, cujas tres temporadas andam no MESMO sentido (+5,02, +3,78,
    +8,17) - e deslocamento de mesmo sinal em todas as arquivais e justamente o offset comum, que
    ali precisa passar de 10 min para derrubar a deteccao. Espero 232634196 sobreviver (o
    deslocamento comum efetivo fica em ~5 min, metade do limiar) e 329246824 ser o mais ameacado.

(c) TAMANHO DO EFEITO. Espero |dPdt(novo) - dPdt(velho)| / sigma com mediana abaixo de 0,5 e
    maximo abaixo de 3. Acima de 3 em qualquer alvo vai ao relatorio nomeado.

(d) AS DUAS TEMPORADAS SEM MEDIDA. 144194304 (temporada 0) e 198388252 (temporada 1) nao tem
    largura livre medivel. Elas entram no reajuste com dt = 0 - que e o valor da cadeia atual - e
    os dois alvos ficam marcados como cobertura PARCIAL. Nenhum dos dois pode ser contado como
    "verificado" sem ressalva, e e por isso que o "11 de 11" da 3.4 precisa de qualificacao.

Uso:
    python src/estimador_alt_refit.py
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

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="backslashreplace")
except (AttributeError, ValueError):
    pass

BASE = config.DATA / "orquestra" / "oc_lote"
SAIDA = BASE / "estimador_alt_refit.parquet"


def ajusta(pt, P):
    """(dPdt, sigma, p_curv, p_adv) com a mesma maquinaria da cadeia."""
    sig_d = pt.sig_min.values / 1440.0
    _, _, chi2l, _ = co.ajustar(pt.E.values, pt.t0.values, sig_d, 1)
    cp, _, chi2p, cov = co.ajustar(pt.E.values, pt.t0.values, sig_d, 2)
    Q, sQ = float(cp[0]), float(np.sqrt(cov[0, 0]))
    dPdt = 2 * Q / P * 365.25 * 86400
    s = 2 * sQ / P * 365.25 * 86400
    p_curv = float(stats.chi2.sf(chi2l - chi2p, 1))
    ag = pt.copy()
    ag["setor"] = range(len(pt))
    adv = co.adversarial(ag, P, float(pt.t0.iloc[-1]))
    return dPdt, s, p_curv, float(adv["p_adversarial"])


def main():
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    S = pd.read_parquet(BASE / "estimador_alt_d11.parquet")
    d11 = sorted(int(t) for t in t26.index[t26.curvatura])
    linhas = []
    for tic in d11:
        j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
        P = float(j["P_escada_d"])
        pt = pd.DataFrame(j["pontos"])
        temporadas = S[S.tic == tic]
        # cada ponto arquival recebe o dt da temporada mais proxima no tempo (t_fixa em BJD)
        desloc, parcial = [], False
        for r in pt.itertuples():
            if not str(r.fonte).startswith("SuperWASP"):
                desloc.append(0.0)
                continue
            t_ponto = float(r.t0) + 2457000.0
            k = (temporadas.t0_fixa_bjd - t_ponto).abs().idxmin()
            dt = float(temporadas.loc[k, "dt_min"])
            if not np.isfinite(dt) or bool(temporadas.loc[k, "na_borda"]):
                dt, parcial = 0.0, True
            desloc.append(dt)
        pt2 = pt.copy()
        pt2["t0"] = pt.t0.values + np.array(desloc) / 1440.0
        a = ajusta(pt, P)
        b = ajusta(pt2, P)
        sobrevive_a = a[2] < 0.05 and a[3] < 0.05
        sobrevive_b = b[2] < 0.05 and b[3] < 0.05
        linhas.append({"tic": tic, "n_pontos": len(pt), "n_arquivais": int(sum(1 for x in desloc if x != 0.0)),
                       "desloc_min": [round(x, 3) for x in desloc if x != 0.0],
                       "desloc_max_min": max((abs(x) for x in desloc), default=0.0),
                       "cobertura_parcial": parcial,
                       "dPdt_cadeia": a[0], "s_cadeia": a[1], "p_curv_cadeia": a[2], "p_adv_cadeia": a[3],
                       "dPdt_alt": b[0], "s_alt": b[1], "p_curv_alt": b[2], "p_adv_alt": b[3],
                       "desloc_sigma": abs(b[0] - a[0]) / a[1],
                       "sobrevive_cadeia": sobrevive_a, "sobrevive_alt": sobrevive_b})
        print(f"  TIC {tic}: dP/dt {a[0]:+.4f} ± {a[1]:.4f} -> {b[0]:+.4f} ± {b[1]:.4f} "
              f"({abs(b[0] - a[0]) / a[1]:.2f}σ) | p_curv {a[2]:.1e} -> {b[2]:.1e} | "
              f"p_adv {a[3]:.3f} -> {b[3]:.3f} | {'SOBREVIVE' if sobrevive_b else 'CAI'}"
              f"{' (cobertura parcial)' if parcial else ''}")
    R = pd.DataFrame(linhas)
    R.to_parquet(SAIDA, index=False)
    print(f"\n  sobrevivem com as epocas da cadeia: {int(R.sobrevive_cadeia.sum())} de {len(R)}")
    print(f"  sobrevivem com as epocas do template alternativo: {int(R.sobrevive_alt.sum())} de {len(R)}")
    print(f"  deslocamento do coeficiente: mediana {R.desloc_sigma.median():.2f}σ, maximo "
          f"{R.desloc_sigma.max():.2f}σ (TIC {int(R.loc[R.desloc_sigma.idxmax(), 'tic'])})")
    print(f"  cobertura parcial (temporada sem largura livre): {sorted(int(t) for t in R.tic[R.cobertura_parcial])}")
    print(f"\n-> {SAIDA}")


if __name__ == "__main__":
    main()
