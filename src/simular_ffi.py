"""Simulacao: que precisao de epoca uma FFI de 30 min entrega?

A PERGUNTA QUE ISTO DECIDE. O TIC 142874476 tem tres aglomerados de epocas -
SuperWASP (2006-2008), setor 31 (2020) e setor 105 (2026) - e tres parametros
numa parabola. Ajuste sem grau de liberdade nao e testado: a parabola passa
pelos tres pontos por construcao, e uma senoide de tempo de luz por terceiro
corpo tambem passa. So epoca INTERMEDIARIA separa as duas hipoteses, e as
candidatas sao FFI: setor 4 (2018, entre o SuperWASP e o s31) e setor 97 (2025,
entre o s31 e o s105).

MAS A CADENCIA DE 30 MIN PODE NAO ENTREGAR. O eclipse dura 1,967 h = 118 min,
entao cabem ~4 pontos nele, e cada ponto e uma INTEGRAL de 30 min sobre a curva -
o borramento nao e ruido, e distorcao deterministica da forma. Antes de extrair
qualquer coisa do s4 e preciso saber que precisao o metodo recupera no minimo
MEDIO do aglomerado. Se der pior que ~2 min, o s4 nao separa parabola de
senoide e so o s97 conta.

O QUE A SIMULACAO PRECISA TER, para nao repetir o erro de gerador que ja
apareceu tres vezes neste projeto: a integracao de 30 min explicita (nao
amostragem pontual), o gap de downlink no meio do setor, ruido correlacionado
alem do branco, e fase inicial sorteada - o minimo nao cai em fase zero no dado
real. Sem esses incomodos a simulacao devolve a precisao de um caso que nao
existe.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import config
import ids as ids_mod

import numpy as np
import pandas as pd

OUT = config.RESULTS / "simular_ffi.parquet"

# TIC 142874476, medidos aqui e na reanalise com o setor 31.
P_D = 1.4740475
T14_H = 1.967
PROF = 0.0316            # 31.600 ppm
CAD_MIN = 30.0
SETOR_D = 27.4
GAP_D = 1.0              # downlink no meio do setor


def modelo(t, t0, P=P_D, t14_h=T14_H, prof=PROF, forma_borda=0.25):
    """Trapezio periodico. `forma_borda` = fracao de T14 em cada borda."""
    t14 = t14_h / 24.0
    fase = ((t - t0 + 0.5 * P) % P) - 0.5 * P     # -P/2 .. P/2, minimo em 0
    x = np.abs(fase) / (t14 / 2)                  # 1 no contato externo
    plano = 1.0 - 2 * forma_borda                 # meia-largura do fundo plano
    y = np.ones_like(t)
    dentro = x < 1
    # fundo plano
    y[dentro & (x <= plano)] = 1 - prof
    # bordas lineares
    b = dentro & (x > plano)
    y[b] = 1 - prof * (1 - x[b]) / (1 - plano)
    return y


def curva_ffi(t0, dias=SETOR_D, cad_min=CAD_MIN, ruido_ppm=1500,
              ruido_vermelho_ppm=500, gap_d=GAP_D, n_sub=30, rng=None):
    """Curva FFI simulada, com a INTEGRACAO da exposicao feita de verdade.

    Cada ponto e a media do modelo sobre a janela de `cad_min`, amostrada em
    `n_sub` subpontos. Amostrar o modelo no centro da janela subestimaria o
    borramento, que e o efeito dominante aqui - e a simulacao devolveria uma
    precisao que o dado real nao tem.
    """
    rng = rng or np.random.default_rng()
    cad = cad_min / 1440.0
    t = np.arange(0, dias, cad)
    meio = dias / 2
    t = t[(t < meio - gap_d / 2) | (t > meio + gap_d / 2)]
    desl = (np.arange(n_sub) + 0.5) / n_sub * cad - cad / 2
    y = np.mean([modelo(t + d, t0) for d in desl], axis=0)
    ruido = rng.normal(0, ruido_ppm / 1e6, len(t))
    if ruido_vermelho_ppm > 0:
        # ruido correlacionado: passeio suavizado na escala de horas
        w = rng.normal(0, 1, len(t))
        k = max(3, int(0.25 / cad))               # ~6 h
        nucleo = np.ones(k) / k
        vermelho = np.convolve(w, nucleo, mode="same")
        vermelho *= (ruido_vermelho_ppm / 1e6) / (np.std(vermelho) or 1)
        ruido = ruido + vermelho
    return t, y + ruido


def medir_epoca(t, f, P=P_D, t0_chute=None, span_min=90.0, n=901):
    """Epoca MEDIA do aglomerado, por minimos quadrados contra o proprio modelo.

    Devolve (t0, sigma_min). O sigma vem do delta-chi2 = 1 na curva de chi2, que
    e a incerteza estatistica correta para um parametro - e nao de std/sqrt(n),
    que subestima com ruido correlacionado.
    """
    if t0_chute is None:
        t0_chute = t[int(np.argmin(f))]
    grade = t0_chute + np.linspace(-span_min, span_min, n) / 1440.0
    var = np.var(f)
    chi2 = np.array([np.sum((f - modelo(t, g)) ** 2) for g in grade]) / var
    i = int(np.argmin(chi2))
    t0 = float(grade[i])
    # largura onde chi2 sobe 1 acima do minimo
    alvo = chi2[i] + 1.0
    esq = grade[:i + 1][chi2[:i + 1] <= alvo]
    dir_ = grade[i:][chi2[i:] <= alvo]
    if len(esq) and len(dir_):
        sig = (dir_.max() - esq.min()) / 2 * 1440.0
    else:
        sig = np.nan
    return t0, float(sig)


def rodar(n_real=200, cenarios=None, semente=20260907):
    """Varre cenarios de ruido e cadencia. Fase sorteada em cada realizacao."""
    if cenarios is None:
        cenarios = [
            ("s4 FFI 30 min, ruido tipico", 30.0, 1500, 500),
            ("s4 FFI 30 min, ruido pessimista", 30.0, 3000, 1000),
            ("s4 FFI 30 min, so ruido branco", 30.0, 1500, 0),
            ("s97 FFI 200 s", 200 / 60.0, 900, 300),
            ("controle: 2 min como o s105", 2.0, 700, 200),
            ("controle negativo: 2 h", 120.0, 1500, 500),
        ]
    rng = np.random.default_rng(semente)
    linhas = []
    for nome, cad, br, verm in cenarios:
        erros, sigmas = [], []
        for _ in range(n_real):
            t0v = float(rng.uniform(0, P_D))      # fase sorteada
            t, f = curva_ffi(t0v, cad_min=cad, ruido_ppm=br,
                             ruido_vermelho_ppm=verm, rng=rng)
            t0m, sig = medir_epoca(t, f)
            # traz para o ciclo mais proximo do valor verdadeiro
            dt = (t0m - t0v + 0.5 * P_D) % P_D - 0.5 * P_D
            erros.append(dt * 1440.0)
            sigmas.append(sig)
        e = np.array(erros)
        linhas.append({
            "cenario": nome, "cadencia_min": cad, "ruido_branco_ppm": br,
            "ruido_vermelho_ppm": verm, "n": n_real,
            "vies_min": float(np.median(e)),
            "espalhamento_min": float(np.std(e)),
            "mad_min": float(np.median(np.abs(e - np.median(e))) * 1.4826),
            "sigma_formal_mediano_min": float(np.nanmedian(sigmas)),
            "p95_abs_min": float(np.percentile(np.abs(e), 95)),
        })
        print(f"  {nome:36s} espalhamento {linhas[-1]['espalhamento_min']:6.2f} min"
              f"   vies {linhas[-1]['vies_min']:+6.2f}", flush=True)
    d = pd.DataFrame(linhas)
    d.to_parquet(OUT, index=False)
    return d


if __name__ == "__main__":
    print(f"TIC 142874476: P = {P_D} d, T14 = {T14_H} h, profundidade "
          f"{PROF*1e6:.0f} ppm")
    print(f"setor de {SETOR_D} d -> {SETOR_D/P_D:.1f} ciclos; "
          f"eclipse com {T14_H*60/CAD_MIN:.1f} pontos em cadencia de "
          f"{CAD_MIN:.0f} min\n")
    d = rodar()
    print()
    print("=" * 96)
    print("PRECISAO DA EPOCA MEDIA DO AGLOMERADO")
    print("=" * 96)
    print(d[["cenario", "cadencia_min", "espalhamento_min", "mad_min",
             "vies_min", "sigma_formal_mediano_min"]].round(3).to_string(index=False))
    s4 = d[d["cenario"].str.startswith("s4 FFI 30 min, ruido tipico")].iloc[0]
    print(f"\n  CRITERIO: o s4 so separa parabola de senoide se a epoca sair "
          f"melhor que ~2 min.")
    print(f"  medido no cenario tipico: {s4['espalhamento_min']:.2f} min "
          f"-> {'SERVE' if s4['espalhamento_min'] < 2 else 'NAO SERVE'}")


# -----------------------------------------------------------------------------
# ONDE ESTA SIMULACAO ACERTOU E ONDE ELA FALHOU
#
# ACERTOU A PRECISAO. Previu espalhamento de 0,47 min para a epoca media do
# aglomerado no s4; as tres reducoes independentes (QLP, eleanor, TGLC) sairam
# com espalhamento de 0,51 min. Previsao feita ANTES da extracao.
#
# ERROU O VIES, E POR CONSTRUCAO. Previu +0,01 min. `src/controle_ffi.py`
# degradou a curva SPOC de 2 min do proprio objeto para bins de 30 min, onde a
# verdade e conhecida, e mediu vies de -0,5 min (s31) e -1,3 min (s105). A 10
# min o vies e +0,1 e 0,0 - desprezivel.
#
# O MOTIVO E O MODELO DESTE ARQUIVO: `modelo()` e um trapezio SIMETRICO sobre
# linha de base PLANA. Integrar uma janela simetrica sobre uma funcao simetrica
# devolve vies zero - o resultado estava embutido na premissa. A estrela real
# tem modulacao de mare e um secundario em volta do eclipse, e e essa assimetria
# LOCAL que desloca o minimo ajustado quando a janela chega a 30 min, um quarto
# da duracao do eclipse.
#
# Quarta ocorrencia de "o gerador de teste precisa do incomodo real", e a mais
# instrutiva: aqui o gerador nao deixou o teste passar por engano - ele mediu
# corretamente a quantidade que sabia medir (precisao) e ficou cego para a que
# nao modelava (acuracia). A licao operacional e que simulacao valida precisao;
# so dado com verdade conhecida valida acuracia.
#
# EFEITO NO RESULTADO. A correcao e +0,9 +- 0,4 min na epoca do s4, e ela EMPURRA
# A CURVATURA PARA CIMA:
#
#   sem correcao            delta chi2  5,93   p = 0,015
#   com +0,9 min            delta chi2 12,23   p = 0,0005
#
# Direcao favoravel foi sorte, nao desenho: um vies de mesmo modulo e sinal
# oposto teria apagado o resultado. O que muda de verdade e que o vies agora e
# MEDIDO em vez de suposto zero.
