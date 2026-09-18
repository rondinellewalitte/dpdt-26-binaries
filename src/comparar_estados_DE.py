# -*- coding: utf-8 -*-
"""Estado D -> estado E: o que a troca do vies da cadeia (-2,14 +- 0,78 -> -1,26 +- 0,75 min)
mudou nos registros dos 88 - e o que NAO pode ter mudado.

Estado E = a constante VIES_CADEIA em `cadeia_oc.py` trocada para a media ponderada do
Passo E da rodada seis (`vies_hj_indep.py`: as mesmas temporadas SuperWASP dos 4 Jupiteres
quentes contra Ivshina & Winn 2022 + Bouma et al. 2020 no WASP-4 b: -1,26 +- 0,75). O
estado D (6b1854e) fica no git; este script le os registros dos dois estados e grava
`comparacao_estado_D_E.parquet` (uma linha por alvo medido) com dP/dt, sigma, p_curv,
p_adv, p_gof, classe e veredito nos dois estados, o deslocamento das epocas SuperWASP e
a diferenca contra a PREVISAO ALGEBRICA (abaixo), em sigmas.

Como o vies entra: `oc_lote` grava cada temporada SuperWASP com t0 = bruto - VIES/1440 e
barra hypot(max(formal, dispersao), VIES_SIG); a escada, as epocas TESS, o gate de
cobertura e a contagem de ciclos nao dependem dele. O GLS e linear nos dados: para barras
fixas, dP/dt e EXATAMENTE linear em delta (conferido no sensib_vies_26 do estado D: nao-
linearidade 8e-9 sigma). Logo o estado E e previsivel por algebra a partir do estado D
antes de rodar, e e isso que se testa.

EXPECTATIVA (escrita e commitada ANTES de rodar; cada desvio vai ao humano):
  Controle TIC 142874476: ANTES (estado D, 6b1854e) A e B REPRODUZIRAM (log
    logs_E/00_controle_antes_estadoD.txt). DEPOIS: contra a v3 (= v2 com os quatro campos
    que carregam a constante propagados por algebra: O-C no linear 22,72 -> 21,84; sigmas
    10,16 -> 9,81; O-C na parabola -5,48 -> -6,36; sigmas 2,45 -> 2,85; barra_total 2,237
    -> 2,226) A e B REPRODUZEM; contra a v2, EXATAMENTE esses quatro campos escalam.
  Os 88: os mesmos 26 medidos; Tabela 2 (primeiro guard dos 62) identica; E identico em
    toda epoca; epocas e barras TESS identicas ao digito (< 1e-6 min); epocas SuperWASP
    deslocadas -0,880 +- 0,001 min (t0 = bruto + 1,26 em vez de bruto + 2,14); barras
    SuperWASP diminuem entre 0,000 e 0,030 min (hypot com 0,75 em vez de 0,78).
  dP/dt dos 26: igual a interpolacao linear do sensib_vies_26 do estado D (sigma_v 0,78,
    delta 0 e -2,14) em delta = -1,26, dentro de 0,05 sigma (a unica diferenca e o piso
    0,75, que muda os pesos SuperWASP em <= 2%). Todos os 26 movem para o lado NEGATIVO;
    mediana 0,27 sigma, maximo 0,95 sigma (230386284); nenhum muda de sinal.
  Vereditos: D11 = os mesmos 11; D9 = 9; p_gof < 0,05 os mesmos 3 (198388252, 229914020,
    237116051); nucleo {198408416, 232634196, 329246824, 392536812} intacto; nenhuma nova;
    classes identicas nos 26. chi2_red mediano 0,70 -> 0,65-0,75; agregado 1,47 -> 1,40-1,55.
  Secao 5 (cada script contra o proprio estado D):
    qual_barra_26: TESS 1,59 -> 1,50-1,70; SuperWASP 1,36 -> 1,25-1,45; 23 adequados 0,68 /
      0,74 -> +- 0,05.
    nulo_barras_reestimadas_26: SW 0,77, TESS 0,98, condicionado 0,77 / 0,81 -> +- 0,02 (o
      nulo e simulado em torno do modelo; so o piso 0,75 entra). nulo_pisos_26: +- 0,02.
    vies_alavanca_26: imparcial < 2%. estrutura_swasp_26: f = 0 reproduz o nulo +- 0,02.
    epsilon_26: Sigma eps 6,4 +- 0,1; eps3 4,25 +- 0,1 (depende das barras, nao do vies).
    completude_secular_26: A50 0,014 +- 0,001; --fator 0,81-0,85 +- 0,03.
    gls_vies_26: nenhum veredito muda; sigma +2% mediana.
    secundarios_26: 6 derivas (com 115244268), valores +- 0,02 min/ano (TESS-only).
    reinflar_tess_26: 1 de 11 em 1,0 (232634196), 0 em 2,6.
    v527_ltte: 93% / 84% +- 2. comparar_brno_26: V564 Dra z -4,8 -> -5,0 a -5,2 (dP/dt
      -0,0207); CV Dra -4,4 -> -4,6 a -4,8; concordancias na janela inalteradas.
    caso_232634196: temporadas -134,6 / -98,6 / -108,3 +- 19,8 (as tres -0,88); primario
      -2,66 e secundario -11,64 iguais (TESS-only). previsao_232634196: separacao 16,9 ->
      16,8-17,2 min; janela +- 2,2.
    tabela_nota: Tabela 4 com as mesmas particoes (23 adequados; 15 com sigma < 0,01; 11 de
      26 sobrevivem); |dP/dt| mediana 0,016 -> 0,015-0,017; teste de sinal igual.
    jitter_232634196: sigma_j 1,18 identico (TESS-only); refit 18: dP/dt -0,212 -> -0,213
      a -0,214 +- 0,020, p_gof 0,25-0,35. jitter_primarios_26 --resumo/--quad: sigma_j por
      alvo, variograma e q identicos; inflacao A 5 / B 4 / C 6 de 11 (+- 1; 232634196 sempre).
    deriva_gp_26: hiperparametros identicos (Matern A 0,835, tau_c 68 d; exp 0,824 / 91 d);
      desenho 7 de 11 - o unico em risco e 198408416 (p_adv 0,028 em -2,14 e 0,073 em 0 ->
      ~0,045 em -1,26: fica por pouco); --por-alvo 6 +- 1; 232634196 com 18 epocas: p_gof
      0,02-0,04, separacao 14,0-14,6 min (6 sigma).
    epsilon_longo_26: controles reproduzem o epsilon_26 novo; Sigma eps 5,79 / 16,50 / 10,81
      +- 0,1; N 1,88 +- 0,03; --gp: 3,88 / 13,66 / 7,18 +- 0,15, N 1,25 +- 0,05, FA 0.
    sensib_vies_26 (grade 0 / -1,26 / -1,56 / -2,14 / -4,28; sigma_v 0,75 / 2,14): ATUAL
      reproduz a Tabela 3 (assercao); as linhas com sigma_v 2,14 reproduzem o estado D ao
      digito nos tres deltas antigos (mesmos pontos brutos); ExoClock -1,56: D11 11/11,
      nucleo intacto na diagonal, GP 7 (198408416 p_adv ~0,038); delta 0: 11/11; -4,28:
      9/11; GP 6 / 7 / 8 em 0 / -2,14 / -4,28.
  Nao rodados (nao dependem da constante ou sao a origem dela): vies_hj_indep (E),
    checa_tempo_swasp (A), vies_forma_swasp (Passo 4), jitter_primarios_26 --medir (799
    setores TESS), caso_232634196_completo (F: refeito dentro do Passo H para o D9).
    Nao gerados (texto/figuras/CDS): tabela_nota --inserir, figuras_nota, epocas_cds_26,
    md_para_tex.

RESULTADO (2026-09-17, contra a expectativa acima; logs em scratchpad/logs_E):
  Controle 142874476 ANTES (estado D, 6b1854e) e DEPOIS (estado E, contra a v3):
  A e B REPRODUZIRAM nos dois. Contra a v2, escalam EXATAMENTE os quatro campos
  previstos (O-C no linear 22,72 -> 21,84; sigmas 10,16 -> 9,81; O-C na parabola
  -5,48 -> -6,36; sigmas 2,45 -> 2,85), com a barra_total 2,237 -> 2,226.
  Os 88: mesmos 26 medidos; Tabela 2 identica (0 alvos mudam de guard); E
  identico; TESS identico ao digito (max |dif| 0,00 min em t0 e barra); SW
  deslocadas -0,8800 min em TODAS (desvio maximo 0,0000); barra SW -0,021 a
  -0,000 min; P da escada identico. dP/dt = a previsao algebrica dentro de
  0,0001 sigma (mediana; max 0,017 em 237116051, cuja barra SW e a menor do
  conjunto): a cadeia e linear no vies, como previsto. Todos os 26 movem para o
  lado negativo (26 de 26), mediana -0,285 sigma, max -0,954 (230386284);
  nenhum muda de sinal; s_E/s_D mediana 0,9989.
  Vereditos: D11 11 -> 11 (mesmos), D9 9, p_gof < 0,05 os mesmos 3, nucleo
  intacto, 0 classes mudam; chi2_red mediano 0,703 -> 0,702; agregado 1,468 ->
  1,444. Tabela 4: |dP/dt| mediana 0,0157 -> 0,0137 (26), 0,0181 -> 0,0175 (23
  adequados), 0,0091 -> 0,0085 (16 com sigma < 0,01); testes de sinal iguais;
  curvatura marginal (0,01 < p < 0,05) troca de 126945917 para 199688409 e
  229461186, nenhuma sobrevive ao adversarial.
  Seccao 5 (todos dentro do esperado): qual_barra SW 1,36 -> 1,31 (23 adequados
  0,74 -> 0,72), TESS 1,59 (0,68); nulos identicos (SW 0,77, TESS 0,98; com
  corte 0,77 / 0,81; percentil do SW medido 32% -> 24%); nulo_pisos 0,78/0,95;
  vies_alavanca imparcial; epsilon 6,38 / eps3 4,25; estrutura f = 0 reproduz o
  nulo; completude A50 0,014 e supressao iguais; --fator 0,85 (TESS parado);
  gls: nenhum veredito muda, sigma +2%; secundarios 6 derivas; reinflar 1 de 11
  em d = 1,0 e 0 em 2,6; v527 93% / 84%; Brno: V564 Dra z -4,82 -> -5,08 no
  registro inteiro e -1,57 na janela, CV Dra -4,37 -> -4,70 (-0,39 na janela),
  |z| <= 2 na janela 4 de 4 (a mediana da diferenca dado-a-dado no SuperWASP cai
  de +1,65 para +0,77 min: o vies novo aproxima as nossas epocas SW das de
  Brno); jitter_232634196 refit 18 -0,2135 +- 0,0202, p_gof 0,288, frase (b);
  GP: hiperparametros identicos por construcao (TESS-only), desenho 7 de 11 (os
  mesmos), 232634196 com 18 epocas p_gof 0,023 e separacao 14,5 min;
  epsilon_longo 5,81 / 16,50 / 10,83, N 1,88; --gp N 1,25 [0,80, 2,46], FA 0;
  sensib: ATUAL reproduz a Tabela 3 (assercao), D11 11/11 em -1,26, 10/11 em
  -1,56 com sigma_v 2,14, nucleo intacto na diagonal em todas.
  DESVIOS (tres, todos para o lado de MAIS sobreviventes):
  1. 377253090 passa a sobreviver a TRES testes que reprovava no estado D: a
     inflacao branca A (5 -> 6 de 11), a quadratica pos-hoc C (6 -> 7) e o GP
     por alvo (6 -> 7). Mecanismo: p_adv 0,0014 -> 0,00013 (o dP/dt vai de
     -0,0057 para -0,0063, +0,59 sigma), o que o poe acima do limiar nos tres.
     Nao estava previsto - a expectativa dizia "A 5 / B 4 / C 6 (+- 1)".
  2. No sensib com sigma_v 2,14, delta -1,56 (ExoClock) derruba 390021728
     (10 de 11) enquanto -1,26 mantem 11: a sensibilidade dos vereditos a
     ESCOLHA entre as duas calibracoes independentes existe, mas so com a barra
     do vies inflada a 2,14 min.
  3. O GP agrupado com sigma_v 2,14 derruba 198408416 tambem em -1,26 e -1,56
     (6 de 11, nucleo quebra), o que no estado D so acontecia em delta 0. O
     nucleo continua intacto no GP com a barra 0,75 usada pela cadeia.
"""
import glob
import io
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cadeia_oc as co  # noqa: E402
import config  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
COMMIT_D = "6b1854e"   # o ultimo commit do ESTADO D (vies -2,14 +- 0,78)
VIES_D, SIG_D = -2.14, 0.78
TOL_INTERP_SIGMA = 0.05


def git_bytes(caminho):
    return subprocess.run(["git", "show", f"{COMMIT_D}:{caminho}"], capture_output=True, check=True, cwd=config.DATA.parent).stdout


def registro_D(tic):
    return json.loads(git_bytes(f"data/orquestra/oc_lote/oc__{tic}.json").decode("utf-8"))


def registro_E(tic):
    return json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))


def erro_curto(j):
    e = str(j.get("erro", ""))
    for k in ("pontos validos", "cobertura por bloco", "veio vazio", "escada nao fechou", "T14 do catalogo inaplicavel"):
        if k in e:
            return k
    return e[:40]


def veredito(j):
    pa = j.get("parabola", {})
    p_gof = float(stats.chi2.sf(pa["chi2"], pa["dof"])) if pa else np.nan
    p_c = j.get("p_curvatura"); p_a = j.get("adversarial", {}).get("p_adversarial")
    det = (p_c is not None) and (p_a is not None) and p_c < 0.05 and p_a < 0.05
    return pa.get("dPdt_s_por_ano"), pa.get("sdPdt_s_por_ano"), p_c, p_a, p_gof, det, j.get("classe")


if __name__ == "__main__":
    assert abs(co.VIES_CADEIA_MIN - (-1.26)) < 1e-9 and abs(co.VIES_CADEIA_SIG - 0.75) < 1e-9, "a cadeia nao esta no estado E"
    tics = sorted(int(Path(f).stem.split("__")[1]) for f in glob.glob(str(BASE / "oc__*.json")))
    assert len(tics) == 88, len(tics)
    okD = {t for t in tics if registro_D(t).get("status") == "ok"}
    okE = {t for t in tics if registro_E(t).get("status") == "ok"}
    print(f"medidos: D {len(okD)}, E {len(okE)}; entram {sorted(okE - okD)}, saem {sorted(okD - okE)}")
    assert okD == okE, "o conjunto dos medidos mudou - escalar antes de qualquer outra coisa"
    dif_guard = [(t, erro_curto(registro_D(t)), erro_curto(registro_E(t))) for t in tics if t not in okD
                 if erro_curto(registro_D(t)) != erro_curto(registro_E(t))]
    print(f"Tabela 2: {len(dif_guard)} alvos mudaram de guard {dif_guard}")

    # a previsao algebrica: interpolacao do sensib_vies_26 do estado D em delta = -1,26
    Sd = pd.read_parquet(io.BytesIO(git_bytes("data/orquestra/oc_lote/sensib_vies_26.parquet")))
    d = Sd[(Sd.metodo == "diagonal") & (Sd.sigma_v_min == SIG_D)].pivot(index="tic", columns="delta_min", values="dPdt")
    prev = d[VIES_D] + ((co.VIES_CADEIA_MIN - VIES_D) / (0.0 - VIES_D)) * (d[0.0] - d[VIES_D])

    linhas = []
    for t in sorted(okD):
        b, e = registro_D(t), registro_E(t)
        pb = pd.DataFrame(b["pontos"]).sort_values("t0").reset_index(drop=True)
        pe = pd.DataFrame(e["pontos"]).sort_values("t0").reset_index(drop=True)
        assert len(pb) == len(pe) and (pb.fonte.values == pe.fonte.values).all(), t
        assert (pb.E.values == pe.E.values).all(), (t, "E mudou")
        tess = pb.fonte.str.startswith("TESS").values
        d_tess = np.abs(pe.t0.values[tess] - pb.t0.values[tess]) * 1440.0
        d_btess = np.abs(pe.sig_min.values[tess] - pb.sig_min.values[tess])
        desl_sw = (pe.t0.values[~tess] - pb.t0.values[~tess]) * 1440.0
        dbar_sw = pe.sig_min.values[~tess] - pb.sig_min.values[~tess]
        dD, sD, pcD, paD, pgD, detD, clD = veredito(b)
        dE, sE, pcE, paE, pgE, detE, clE = veredito(e)
        r = {"tic": t, "n_tess": int(tess.sum()), "n_sw": int((~tess).sum()),
             "max_dif_t0_tess_min": float(d_tess.max()), "max_dif_barra_tess_min": float(d_btess.max()),
             "desl_sw_min_media": float(desl_sw.mean()), "desl_sw_min_max_desvio": float(np.abs(desl_sw - desl_sw.mean()).max()),
             "dbar_sw_min_min": float(dbar_sw.min()), "dbar_sw_min_max": float(dbar_sw.max()),
             "P_escada_D": b["P_escada_d"], "P_escada_E": e["P_escada_d"],
             "dPdt_D": dD, "dPdt_E": dE, "s_D": sD, "s_E": sE, "dPdt_E_previsto": float(prev.loc[t]) if t in prev.index else np.nan,
             "p_D": pcD, "p_E": pcE, "padv_D": paD, "padv_E": paE, "pgof_D": pgD, "pgof_E": pgE, "detecta_D": detD, "detecta_E": detE,
             "classe_D": clD, "classe_E": clE, "dof": e.get("dof"),
             "chi2red_D": b["parabola"]["chi2"] / b["parabola"]["dof"] if b.get("parabola") else np.nan,
             "chi2red_E": e["parabola"]["chi2"] / e["parabola"]["dof"] if e.get("parabola") else np.nan}
        if dD is not None and dE is not None:
            r["dif_E_menos_D_sigma"] = (dE - dD) / sD
            r["dif_contra_previsao_sigma"] = (dE - r["dPdt_E_previsto"]) / sD if np.isfinite(r["dPdt_E_previsto"]) else np.nan
        linhas.append(r)
    R = pd.DataFrame(linhas)
    R.to_parquet(BASE / "comparacao_estado_D_E.parquet", index=False)

    print(f"\n== TESS: max |dif t0| {R.max_dif_t0_tess_min.max():.2e} min, max |dif barra| {R.max_dif_barra_tess_min.max():.2e} min (esperado 0)")
    print(f"== SuperWASP: deslocamento medio {R.desl_sw_min_media.mean():+.4f} min (esperado -0,880), max desvio de uma epoca do valor comum {R.desl_sw_min_max_desvio.max():.4f}; "
          f"barra muda entre {R.dbar_sw_min_min.min():+.4f} e {R.dbar_sw_min_max.max():+.4f} min (esperado -0,030 a 0)")
    print(f"== P da escada: max |dif| {(R.P_escada_E - R.P_escada_D).abs().max():.2e} d (esperado 0)")
    q = R.dropna(subset=["dif_contra_previsao_sigma"])
    print(f"== dP/dt contra a previsao algebrica (interpolacao do estado D em -1,26): |dif| mediana {q.dif_contra_previsao_sigma.abs().median():.4f} sigma, "
          f"max {q.dif_contra_previsao_sigma.abs().max():.4f} ({int(q.loc[q.dif_contra_previsao_sigma.abs().idxmax(), 'tic'])}); tolerancia {TOL_INTERP_SIGMA}: "
          f"{'DENTRO' if q.dif_contra_previsao_sigma.abs().max() <= TOL_INTERP_SIGMA else 'FORA - escalar'}")
    print(f"== dP/dt E - D em sigma_D: mediana {q.dif_E_menos_D_sigma.median():+.3f}, min {q.dif_E_menos_D_sigma.min():+.3f}, max {q.dif_E_menos_D_sigma.max():+.3f}; "
          f"negativos {int((q.dif_E_menos_D_sigma < 0).sum())} de {len(q)}; mudam de sinal: {sorted(int(t) for t in q[np.sign(q.dPdt_D) != np.sign(q.dPdt_E)].tic)}; "
          f"razao s_E / s_D mediana {(q.s_E / q.s_D).median():.4f}")
    print(f"== curvatura + adversarial: D {int(R.detecta_D.sum())} -> E {int(R.detecta_E.sum())}; entram {sorted(int(t) for t in R[R.detecta_E & ~R.detecta_D].tic)}; "
          f"saem {sorted(int(t) for t in R[R.detecta_D & ~R.detecta_E].tic)}")
    print(f"== p_gof < 0,05: D {sorted(int(t) for t in R[R.pgof_D < 0.05].tic)} -> E {sorted(int(t) for t in R[R.pgof_E < 0.05].tic)}")
    print(f"== classes que mudam: {[(int(r.tic), r.classe_D, r.classe_E) for r in R.itertuples() if r.classe_D != r.classe_E]}")
    print(f"== chi2_red mediano: D {R.chi2red_D.median():.3f} -> E {R.chi2red_E.median():.3f}; agregado D {(R.chi2red_D * R.dof).sum() / R.dof.sum():.3f} -> E {(R.chi2red_E * R.dof).sum() / R.dof.sum():.3f}")
    nuc = (198408416, 232634196, 329246824, 392536812)
    print(f"== nucleo: " + ", ".join(f"{t} {'ok' if bool(R.set_index('tic').loc[t, 'detecta_E']) else 'CAI'}" for t in nuc))
    print(R[["tic", "dPdt_D", "dPdt_E", "dPdt_E_previsto", "s_D", "dif_E_menos_D_sigma", "dif_contra_previsao_sigma", "p_E", "padv_E", "pgof_E", "detecta_D", "detecta_E"]].round(5).to_string(index=False))
    print("gravado:", BASE / "comparacao_estado_D_E.parquet")
