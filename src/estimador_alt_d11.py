# -*- coding: utf-8 -*-
"""Rodada vinte: o estimador ALTERNATIVO de timing, com a largura do perfil LIVRE.

A PERGUNTA. Toda a maquinaria do vies de forma - o estado (F), a varredura de offset comum, a
coluna da Tabela 3, um dos dois cortes da 5.3.2 - existe porque as epocas arquivais saem de um
trapezio de largura FIXA: T14 e medida uma vez na curva inteira (`superwasp.largura_do_minimo`) e
depois mantida em todas as temporadas. Isso e um erro de template, e ate aqui ele nunca foi
confrontado com outro metodo. A verificacao de primeira ordem que falta e simples: refazer as
mesmas epocas com a largura LIVRE, nos mesmos dados e nas mesmas temporadas, e reportar a
diferenca de epoca.

O QUE MUDA E O QUE NAO MUDA. Muda um parametro: T14 deixa de ser imposta e passa a ser ajustada
junto com t0 e a profundidade, numa grade 2-D (t0 x T14). Nao muda nada mais - a mesma curva
limpa, o mesmo corte em temporadas (gap > 30 d), a mesma semente, a mesma janela de meio periodo,
o mesmo criterio de cobertura de fase. Assim a diferenca de epoca que sair e atribuivel ao
template, e nao a outra coisa.

ALVOS. Os onze do D11, que sao os alvos sobre os quais o artigo afirma alguma coisa. O custo e o
de reajustar ~2-4 temporadas por alvo numa grade 2-D.

EXPECTATIVA (escrita e commitada ANTES de rodar)

(a) LARGURA. Espero a largura livre dentro de 0,8 a 1,3 vezes a fixa na mediana das temporadas -
    `largura_do_minimo` mede na propria curva, entao nao deve estar longe. Espero tambem que uma
    fracao das temporadas NAO seja medivel com largura livre (a largura corre para a borda da
    grade porque o eclipse esta mal coberto): 10 a 30% delas. Essas entram como nao-medida
    declarada, nunca como valor cortado na borda.

(b) DIFERENCA DE EPOCA. Espero |dt| = |t0(livre) - t0(fixa)| com mediana abaixo de 1,5 min e
    percentil 90 abaixo de 4 min sobre as temporadas medidas. A escala de comparacao e a propria
    varredura de offset: a Tabela 3 diz, por alvo, o menor offset comum que derruba a deteccao
    (0,75 a >10 min). Se |dt| ficar sistematicamente abaixo desses valores, o erro de template
    esta coberto pela varredura que o artigo ja publica; se ficar acima em alvos cujo offset
    minimo e pequeno, a varredura NAO cobre e o texto tem de dizer isso.

(c) EFEITO NAS DETECCOES. Espero que NENHUMA deteccao com offset minimo >= 3 min mude de classe,
    e que os alvos com offset minimo de 0,75 a 2 min (390021728 entre eles) sejam os unicos em
    risco. Espero no maximo 2 dos 11 mudando de classe. Se mudarem 3 ou mais, ou se a mediana de
    |dt| passar de 3 min, PARO e reporto antes de tocar em qualquer texto.

(d) SINAL. Nao espero vies de sinal em dt (o template errado desloca para os dois lados conforme
    a assimetria da cobertura): espero a mediana de dt compativel com zero dentro da dispersao,
    e e por isso que a comparacao util e |dt| e nao dt.

Uso:
    python src/estimador_alt_d11.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402
import superwasp as sw  # noqa: E402
import epocas_142874476 as ep  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="backslashreplace")
except (AttributeError, ValueError):
    pass

BASE = config.DATA / "orquestra" / "oc_lote"
SAIDA = BASE / "estimador_alt_d11.parquet"
FATORES = np.geomspace(0.45, 2.2, 25)     # T14 livre: a fixa vezes isto


def epoca_largura_livre(t, f, P, t0_prev, t14_fix):
    """(t0, sigma, T14 ajustada, encostou_na_borda) com a largura LIVRE.

    Mesma `medir_epoca` do estimador da cadeia, varrida sobre uma grade de larguras: o minimo
    global em (t0, T14) e a medida. A borda da grade e reportada, nunca cortada em silencio."""
    melhor = None
    for fator in FATORES:
        t14 = t14_fix * fator
        t0, sig, prof = sw.medir_epoca(t, f, P, t0_prev, t14_h=t14)
        if not np.isfinite(t0) or not np.isfinite(prof) or prof <= 0:
            continue
        m = ep.sim.modelo(t, t0, P=P, t14_h=t14, prof=prof) - 1.0
        chi2 = float(np.sum((f - 1.0 - m) ** 2))
        if melhor is None or chi2 < melhor[0]:
            melhor = (chi2, t0, sig, t14, fator)
    if melhor is None:
        return np.nan, np.nan, np.nan, True
    _, t0, sig, t14, fator = melhor
    borda = bool(fator <= FATORES[0] * 1.001 or fator >= FATORES[-1] * 0.999)
    return t0, sig, t14, borda


def main():
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    alvos = pd.read_parquet(config.DATA / "orquestra" / "alvos_oc_26.parquet").set_index("tic")
    d11 = sorted(int(t) for t in t26.index[t26.curvatura])
    print(f"== estimador alternativo (largura livre) nos {len(d11)} alvos do D11")
    linhas = []
    for tic in d11:
        j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
        sws = j["superwasp"]
        P = float(j["P_escada_d"])
        # a largura FIXA da cadeia e a do CATALOGO (oc_lote passa r.t14_h): e este o template
        t14_fix = float(alvos.loc[tic, "t14_h"])
        t, f, _ = sw.carregar(sws["sourceid"], float(alvos.loc[tic, "ra"]), float(alvos.loc[tic, "dec"]))
        o = np.argsort(t)
        t, f = t[o], f[o]
        corte = np.where(np.diff(t) > 30.0)[0]
        grupos = [g for g in np.split(np.arange(len(t)), corte + 1) if len(g) > 200]
        t0_prev = float(sws["temporadas"][0]["t0"])
        for k, g in enumerate(grupos):
            t0_fix, sig_fix, prof = sw.medir_epoca(t[g], f[g], P, t0_prev, t14_h=t14_fix)
            t0_liv, sig_liv, t14_liv, borda = epoca_largura_livre(t[g], f[g], P, t0_prev, t14_fix)
            dt = (t0_liv - t0_fix) * 1440.0 if np.isfinite(t0_liv) and np.isfinite(t0_fix) else np.nan
            linhas.append({"tic": tic, "temporada": k, "n_pontos": len(g),
                           "t0_fixa_bjd": t0_fix, "sig_fixa_min": sig_fix, "prof_fixa": prof,
                           "t0_livre_bjd": t0_liv, "sig_livre_min": sig_liv,
                           "t14_fixa_h": t14_fix, "t14_livre_h": t14_liv,
                           "razao_t14": t14_liv / t14_fix if np.isfinite(t14_liv) else np.nan,
                           "na_borda": borda, "dt_min": dt})
            print(f"  TIC {tic} temporada {k} (n={len(g):5d}): T14 {t14_fix:.2f} -> "
                  f"{t14_liv:.2f} h ({t14_liv / t14_fix:.2f}x){' BORDA' if borda else ''} | "
                  f"dt = {dt:+.2f} min" if np.isfinite(dt) else
                  f"  TIC {tic} temporada {k} (n={len(g):5d}): sem medida")
    R = pd.DataFrame(linhas)
    R.to_parquet(SAIDA, index=False)
    ok = R[np.isfinite(R.dt_min) & ~R.na_borda]
    print(f"\n  temporadas: {len(R)} | medidas com largura livre fora da borda: {len(ok)} "
          f"({len(ok) / len(R):.0%}) | na borda: {int(R.na_borda.sum())}")
    if len(ok):
        print(f"  razao T14 livre/fixa: mediana {ok.razao_t14.median():.2f} "
              f"(min {ok.razao_t14.min():.2f}, max {ok.razao_t14.max():.2f})")
        print(f"  |dt|: mediana {ok.dt_min.abs().median():.2f} min, p90 {ok.dt_min.abs().quantile(0.9):.2f}, "
              f"max {ok.dt_min.abs().max():.2f} (TIC {int(ok.loc[ok.dt_min.abs().idxmax(), 'tic'])})")
        print(f"  dt com sinal: mediana {ok.dt_min.median():+.2f} min")
    print(f"\n-> {SAIDA}")


if __name__ == "__main__":
    main()
