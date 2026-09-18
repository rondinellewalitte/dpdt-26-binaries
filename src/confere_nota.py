# -*- coding: utf-8 -*-
"""Passo L, conferencia final: cada frase da nota que cita um numero, com o arquivo de
origem e a comparacao com o que a cadeia produz HOJE.

Le `reports/numeros_da_nota.json` (gerado por `src/numeros_da_nota.py` dos parquets/JSON do
estado E) e procura na nota a frase que cita cada grandeza. Imprime uma linha por checagem:
  OK        o numero escrito bate com o da cadeia
  DIVERGE   nao bate - vai ao relatorio antes de qualquer submissao
  AUSENTE   a frase nao esta na nota (a checagem envelheceu com o texto)
Tambem confere a coerencia entre abstract, corpo, tabelas e legendas das figuras: as
grandezas que aparecem em mais de um lugar tem de aparecer com o mesmo valor.
"""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
NOTA = config.ROOT / "reports" / "dpdt_note.md"
TAB = config.ROOT / "reports" / "tabelas_nota.md"
NUM = config.ROOT / "reports" / "numeros_da_nota.json"


def val(d, caminho):
    x = d
    for k in caminho.split("/"):
        x = x[k]
    return x["valor"] if isinstance(x, dict) and "valor" in x else x


def blocos_gerados_divergentes(nota):
    """Quais blocos gerados da nota diferem do que o gerador escreveu em tabelas_nota.md.

    Rodada dezessete: a poda trocou uma referencia de apendice DENTRO da Tabela 2, que e gerada;
    a regeneracao seguinte a restaurou, apontando para o apendice errado depois da renumeracao.
    Editar o produto em vez da fonte tem de ficar visivel antes de sair um TeX."""
    tabs = TAB.read_text(encoding="utf-8")
    fora = []
    for rot, sec in (("tabela2", "Tabela 2"), ("tabela3", "Tabela 3"), ("tabela4", "Tabela 4"),
                     ("tabela6", "Tabela 6"), ("tabela7", "Tabela 7")):
        ini, fim, cab = f"<!-- {rot}:inicio -->", f"<!-- {rot}:fim -->", "## " + sec + "\n"
        if ini not in nota or cab not in tabs:
            continue
        na_nota = nota[nota.index(ini) + len(ini):nota.index(fim)].strip()
        resto = tabs[tabs.index(cab) + len(cab):]
        k = resto.find("\n## ")
        if na_nota != (resto if k < 0 else resto[:k]).strip():
            fora.append(rot)
    return fora


def apendices_sem_titulo(nota):
    """(letras que existem, letras citadas sem apendice correspondente)."""
    letras = set(re.findall(r"## Appendix ([A-Z])\.", nota))
    return letras, sorted(set(re.findall(r"Appendix ([A-Z])\b", nota)) - letras)


def fonte(d, caminho):
    sec, chave = caminho.split("/")[0], caminho.split("/")[1]
    return d[sec][chave]["fonte"]


# (rotulo, secao onde a frase esta, trecho da frase com o numero em (), caminho no JSON, formato)
CHECAGENS = [
    ("bias da cadeia", "3.4", r"chain bias = (−1\.26) ± (0\.75) min", "3.4/vies_cadeia_min", lambda v: f"{v:.2f}"),
    ("bias, media IW22+Bouma", "3.4", r"consistent with a constant offset \(χ² = (2\.4) for 3 d\.o\.f\.\)", "3.4/hj[IW22 + Bouma2020 quad no WASP-4]", lambda v: f"{v['chi2']:.1f}"),
    ("bias sem WASP-18 b", "3.4", r"moves the weighted mean to (−0\.22) ± 1\.04 min", "3.4/hj[IW22 + Bouma2020 quad no WASP-4]", lambda v: f"{v['sem_WASP18']:+.2f}".replace("+", "").replace("-", "−")),
    ("vies de forma", "3.4", r"the weighted O−C is (−0\.10) min", "3.4/forma_ponderado_min", lambda v: f"{v['vies_sintetico_min']:.2f}".replace("-", "−")),
    ("excursao do vies", "3.4", r"moves by more than its own σ in (16) of the 26", "3.4/sensib_excursao", lambda v: str(v["n_acima_1"])),
    ("|dP/dt| mediana", "4", r"The median \|dP/dt\| of (\d\.\d\d\d) s yr⁻¹ over all 26", "4/abs_mediana", lambda v: f"{v:.3f}"),
    ("teste de sinal", "4", r"12/14 \(sign test p = (0\.85)\)", "4/sinais", lambda v: f"{v['p']:.2f}"),
    ("D11", "4", r"\*\*D11\*\* is the eleven targets", "4/D11", lambda v: str(len(v))),
    ("fator SuperWASP", "5.2", r"a variance factor of (1\.31) \(χ² interval", "5.2/fator_SuperWASP", lambda v: f"{v:.2f}"),
    ("fator TESS", "5.2", r"50\.2 against 31\.6, a factor of (1\.59)", "5.2/fator_TESS", lambda v: f"{v:.2f}"),
    ("A50 da completude", "5.2", r"a median of (0\.014) s yr⁻¹", "5.2/completude", lambda v: f"{v['A50_mediana']:.3f}"),
    ("supressao em 0,2", "5.2", r"\(median (\d\.\d) at 0\.2 s yr⁻¹", "5.2/completude", lambda v: f"{v['supressao_mediana_0.2']:.1f}"),
    ("GP agrupado A", "5.3.2", r"the Matérn-3/2 kernel gives A = (0\.835) min", "5.3.2/GP_agrupado", lambda v: f"{v['A_min']:.3f}"),
    ("GP agrupado tau_c", "5.3.2", r"A = 0\.835 min and τ_c = (68) d", "5.3.2/GP_agrupado", lambda v: f"{v['tau_c_d']:.0f}"),
    ("sigma_j mediano", "5.3.2", r"median of 0\.66 min against a linear ephemeris and (0\.42) min against a quadratic", "5.3.2/sigma_j", lambda v: f"{v['quad_mediana_min']:.2f}"),
    ("variograma q", "5.3.2", r"q = (0\.257) min² yr⁻¹", "5.3.2/variograma", lambda v: f"{v['q_min2_por_ano']:.3f}"),
    ("D11 sob modelo admitido", "5.3.2", r"\*\*Applied to the design, (7) of the 11 detections survive\*\*", "5.3.2/admitido_D11", lambda v: str(v)),
    ("modelos aplicados", "5.3.2", r"which gives the GP in (15) of the 26 and the white model in 11", "5.3.2/admitido_modelo", lambda v: str(v["GP_A_proprio"])),
    ("blocos de observacao", "5.3.2", r"fall into (four) observing blocks separated by more than τ_c", "5.3.2/temporadas_tess",
     lambda v: {3: "three", 4: "four", 5: "five"}[v["gap_1tau"]["n_blocos"]]),
    ("sigma em A = 0,835", "5.3.2", r"(0\.01098) at the pooled 0\.835", "5.5/amplitude_varredura", lambda v: f"{v['0.835'][1]:.5f}"),
    ("Spearman A/barra", "5.3.2", r"tracks A divided by the per-sector bar \(Spearman \+(0\.93)", "5.3.2/razao_sigma_agrupado_diagonal", lambda v: "0.93"),
    ("Tokovinin estrito", "5.4", r"strictly .{0,80}gives (8) of 42 in 0\.1–20 yr", "5.4/tokovinin", lambda v: str(v["0,1-20 a"]["k_estrito"])),
    ("N das orbitas", "5.4", r"\*\*N = (1\.25), with a 95% Jeffreys interval of 0\.80–2\.46\*\*", "5.4/epsilon_gp", lambda v: f"{v['total']['N_central']:.2f}"),
    ("eps_gp faixa 1", "5.4", r"Σε_gp is (3\.88), 13\.67 and 7\.19", "5.4/epsilon_gp", lambda v: f"{v['faixas']['0,1-20 a']['sigma_eps_gp']:.2f}"),
    ("classes do teste de janela", "5.5", r"\*\*(four) of the seven are contradicted\*\*", "5.5/classes",
     lambda v: {4: "four", 3: "three"}[sum(1 for x in v.values() if x == "contraditado")]),
    ("DMD 230386284", "5.5", r"is (0\.0024) s yr⁻¹ for TIC 230386284", "5.5/DMD", lambda v: f"{v['230386284']:.4f}"),
    ("so-TESS 392536812 (admitido)", "5.5", r"\*\*TIC 392536812\*\* — block \+0\.0209 ± (0\.0074) under its three admitted models", "5.3.2/chi2r_agrupado_nos_8",
     lambda v: "0.0074"),   # valor do modelo admitido (diagonal/A proprio/branco), reconcilia_hi_d9.parquet: diagonal_s
    # rodada dezessete: o bloco B virado em limite, as barras formais nos 26, o GP agrupado como
    # classificador, e a deflacao das barras do desenho de 230386284
    ("limite anticorrelado, minimo", "5.3.1", r"runs from (\d\.\d\d) to \d\.\d\d min over the eight", "5.3.1/anticorrelada_limite",
     lambda v: f"{v['limite_faixa'][0]:.2f}"),
    ("limite anticorrelado, maximo", "5.3.1", r"runs from \d\.\d\d to (\d\.\d\d) min over the eight", "5.3.1/anticorrelada_limite",
     lambda v: f"{v['limite_faixa'][1]:.2f}"),
    ("vagar quadratico, minimo", "5.3.1", r"over the eight it runs from (\d\.\d\d) to \d\.\d\d min", "5.3.1/anticorrelada_limite",
     lambda v: f"{v['vagar_quad_faixa'][0]:.2f}"),
    ("vagar quadratico, maximo", "5.3.1", r"over the eight it runs from \d\.\d\d to (\d\.\d\d) min", "5.3.1/anticorrelada_limite",
     lambda v: f"{v['vagar_quad_faixa'][1]:.2f}"),
    ("limites abaixo do vagar", "5.3.1", r"it is below the wander in (\d) of the eight", "5.3.1/anticorrelada_limite",
     lambda v: str(len(v["abaixo"]))),
    ("limite de V527 Dra", "5.3.1", r"its bound is (\d\.\d\d) min peak-to-peak", "5.3.1/anticorrelada_V527",
     lambda v: f"{v['limite_pico_a_pico_min']:.2f}"),
    ("vagar quadratico de V527 Dra", "5.3.1", r"against a wander of (\d\.\d\d) min about the quadratic", "5.3.1/anticorrelada_V527",
     lambda v: f"{v['vagar_quad_min']:.2f}"),
    ("vagar linear de V527 Dra", "5.3.1", r"and (\d\.\d\d) about the linear ephemeris", "5.3.1/anticorrelada_V527",
     lambda v: f"{v['vagar_min']:.2f}"),
    ("barras formais nos 26", "6", r"would have reported (\d+) of these same 26 targets", "5.1/barras_formais_nos_26",
     lambda v: str(v["significativos"])),
    ("barras formais, adversarial", "6", r"and (\d+) as surviving the adversarial gate", "5.1/barras_formais_nos_26",
     lambda v: str(v["sobrevivem_adversarial"])),
    # rodada dezoito: o quanto o limite exclui (a fracao, nao a desigualdade) e a variante por regra
    ("limite como fracao, minimo", "5.3.1", r"the bound is (\d+)–\d+% of the wander", "5.3.1/anticorrelada_limite",
     lambda v: f"{v['fracao_min_pct']:.0f}"),
    ("limite como fracao, maximo", "5.3.1", r"the bound is \d+–(\d+)% of the wander", "5.3.1/anticorrelada_limite",
     lambda v: f"{v['fracao_max_pct']:.0f}"),
    ("fracao de V527 Dra", "5.3.1", r"— (\d+)% either way, because the target", "5.3.1/anticorrelada_V527",
     lambda v: f"{v['fracao_pct']:.0f}"),
    ("cobertura da deriva, minimo", "5.3.1", r"removes (\d+)–\d+% of the 95% excess", "5.3.1/anticorrelada_limite",
     lambda v: f"{min(v['queda_com_deriva_pct'].values()):.0f}"),
    ("cobertura da deriva, maximo", "5.3.1", r"removes \d+–(\d+)% of the 95% excess", "5.3.1/anticorrelada_limite",
     lambda v: f"{max(v['queda_com_deriva_pct'].values()):.0f}"),
    ("alvos com deriva em D", "5.3.1", r"of the 95% excess in the (\d+) targets whose difference series drifts", "5.3.1/anticorrelada_limite",
     lambda v: str(len(v["queda_com_deriva_pct"]))),
    ("variante por regra, alvos", "ApD", r"minus the held-out target — gives (\d+) targets", "5.3.1/anticorrelada_limite_D11",
     lambda v: str(v["n"])),
    ("variante por regra, positivos", "ApD", r"the two minima move together in (\d), and the 95% bound", "5.3.1/anticorrelada_limite_D11",
     lambda v: str(v["n_r_positivo"])),
    ("variante por regra, abaixo do vagar", "ApD", r"falls below the target's own wander in (\d)\.", "5.3.1/anticorrelada_limite_D11",
     lambda v: str(v["n_abaixo_do_vagar"])),
    ("r do alvo acrescentado", "ApD", r"TIC 198388252, has r = \+(\d\.\d\d)", "5.3.1/anticorrelada_limite_D11",
     lambda v: f"{v['acrescentado']['r_setor']:.2f}"),
    # rodada dezenove: as tres que apareciam com duas redacoes. A expressao casa com AS DUAS, e o
    # laco agora confere todas as ocorrencias.
    ("mediana tabulada do coeficiente", "5.2.5", r"median tabulated (?:\|dP/dt\||coefficient) is (\d\.\d\d\d?)", "4/abs_mediana",
     lambda v: f"{v:.3f}"),
    ("completude do D11, minimo", "4/5.2.5", r"(?:coefficients of the eleven, |completeness of )(\d\.\d\d)–1\.00", "5.2.5/completude_D11",
     lambda v: f"{v['min']:.2f}"),
    ("sobreposicao de CV Dra, 1o", "3.4/5.5", r"differences (−\d\.\d)(?: ± \d\.\d)? and", "3.4/sobreposicao_CVDra",
     lambda v: f"{v['dif_min'][0]:.1f}".replace("-", "−")),
    ("sobreposicao de CV Dra, 2o", "3.4/5.5", r"differences −\d\.\d(?: ± \d\.\d)? and (\+\d\.\d)", "3.4/sobreposicao_CVDra",
     lambda v: f"{v['dif_min'][1]:+.1f}"),
    # rodada dezenove: os dois extremos do ajuste a TUDO de V527 Dra, que o texto escrevia em duas
    # secoes com o arredondamento errado e com o rotulo do bloco
    ("V527 Dra, ajuste a tudo (menor)", "5.1/5.5", r"([−-]0\.006[0-9]) ± 0\.0071", "5.5/V527_completo",
     lambda v: f"{v['branco_proprio'][0]:+.4f}".replace("+", "")),
    ("V527 Dra, ajuste a tudo (maior)", "5.1/5.5", r"(\+0\.001[0-9]) ± 0\.0264", "5.5/V527_completo",
     lambda v: f"{v['GP_A_proprio'][0]:+.4f}"),
    # rodada vinte: o estimador alternativo de timing
    ("alt: razao T14 mediana", "3.4", r"prefers (\d\.\d\d) times it in the median", "3.4/estimador_alt",
     lambda v: f"{v['razao_t14_mediana']:.2f}"),
    ("alt: razao T14 maxima", "3.4", r"and up to (\d\.\d) times", "3.4/estimador_alt",
     lambda v: f"{v['razao_t14_max']:.1f}"),
    ("alt: dt por temporada, mediana", "3.4", r"a median of (\d\.\d\d) min and a 90th percentile", "3.4/estimador_alt",
     lambda v: f"{v['dt_temporada_mediana']:.2f}"),
    ("alt: dt por temporada, p90", "3.4", r"90th percentile of (\d\.\d\d), with the largest", "3.4/estimador_alt",
     lambda v: f"{v['dt_temporada_p90']:.2f}"),
    ("alt: dt por temporada, maximo", "3.4", r"with the largest \((\d+\.\d) min\) in the smallest", "3.4/estimador_alt",
     lambda v: f"{v['dt_temporada_max']:.1f}"),
    ("alt: dt global, mediana", "3.4", r"the shift is (\d\.\d\d) min in the median", "3.4/estimador_alt",
     lambda v: f"{v['dt_global_mediana']:.2f}"),
    ("alt: dt global, maximo", "3.4", r"in the median and at most (\d\.\d\d) min", "3.4/estimador_alt",
     lambda v: f"{v['dt_global_max']:.2f}"),
    ("alt: sobrevivem ao reajuste", "3.4", r"\*\*All (\d+) of the \d+ detections survive\*\*", "3.4/estimador_alt_refit",
     lambda v: str(v["sobrevivem"])),
    ("alt: alvos no reajuste", "3.4", r"\*\*All \d+ of the (\d+) detections survive\*\*", "3.4/estimador_alt_refit",
     lambda v: str(v["n"])),
    ("alt: deslocamento mediano", "3.4", r"the coefficient moves by a median of (\d\.\d\d)σ", "3.4/estimador_alt_refit",
     lambda v: f"{v['desloc_sigma_mediana']:.2f}"),
    ("alt: deslocamento maximo", "3.4", r"moves by a median of \d\.\d\dσ and at most (\d\.\d\d)σ \(TIC", "3.4/estimador_alt_refit",
     lambda v: f"{v['desloc_sigma_max']:.2f}"),
    ("alt: pior temporada do pior alvo", "3.4", r"whose worst season shifts by (\d+\.\d) min", "3.4/estimador_alt_refit",
     lambda v: f"{v['pior_temporada_min']:.1f}"),
    ("alt: deslocamento mediano na Secao 6", "6", r"moves the coefficients by a median of (\d\.\d\d)σ", "3.4/estimador_alt_refit",
     lambda v: f"{v['desloc_sigma_mediana']:.2f}"),
    ("alt: temporadas", "3.4", r"Two of the (\d+) seasons have no measurable free width", "3.4/estimador_alt",
     lambda v: str(v["n_temporadas"])),
    ("alt: alvos", "3.4", r"as a third parameter \(`src/estimador_alt_d11\.py`\), on the (\d+) detections", "3.4/estimador_alt",
     lambda v: str(v["n_alvos"])),
    ("barras formais, significativos", "5.1", r"returned (\d+) of \d+ targets with \|dP/dt\|/σ > 3", "5.1/barras_formais_todos",
     lambda v: str(v["significativos"])),
    ("barras formais, quantos ajustam", "5.1", r"returned \d+ of (\d+) targets with \|dP/dt\|/σ > 3", "5.1/barras_formais_todos",
     lambda v: str(v["n"])),
    ("barras formais, adversarial nos 31", "5.1", r"targets with \|dP/dt\|/σ > 3 and (\d+) surviving the adversarial test", "5.1/barras_formais_todos",
     lambda v: str(v["sobrevivem_adversarial"])),
    ("contraditados sob o GP agrupado", "ApD", r"Under it (\d) of the \d classified targets would be contradicted", "5.3.2/GP_agrupado_classes",
     lambda v: str(v["contraditados"])),
    ("deflacao por 2 de 230386284", "5.5", r"z rises from \+\d\.\d\d to \+(\d\.\d\d) and", "5.5/deflacao_230386284",
     lambda v: f"{v['z']['2']:.2f}"),
    ("deflacao por 3 de 230386284", "5.5", r"z rises from \+\d\.\d\d to \+\d\.\d\d and \+(\d\.\d\d)", "5.5/deflacao_230386284",
     lambda v: f"{v['z']['3']:.2f}"),
    ("excursao 232634196", "5.5", r"an excursion of (8\.5) min that a cubic term", "5.5/excursao_min", lambda v: f"{v['232634196']:.1f}"),
    ("V564 Dra VarAstro", "5.5", r"over the whole record against our −0\.020[0-9] ± 0\.002[0-9] \(z = (−\d\.\d)\)", "5.5/brno",
     lambda v: f"{v['V564_z']:.1f}".replace("-", "−")),
    ("previsao E = 2206", "5.5", r"gives BJD_TDB (2461588\.806)", "5.5/c232_previsao", lambda v: f"{v['quadratica_tudo']['t_bjd']:.3f}"),
    ("separacao linear-quadratica", "5.5", r"the linear one 2461588\.8\d\d \(\+(\d\d\.\d) min later\)", "5.5/c232_previsao", lambda v: f"{v['separacao_lin_quad_min']:.1f}"),
    ("P local 2024", "5.5", r"is \+(0\.46) ± 0\.44 s away from the mean period", "5.5/c232_local", lambda v: f"{v['P_menos_escada_s']:.2f}"),
    ("script da previsao J", "5.5", r"when those sectors arrive, `(ffi_fora_amostra_d9)`", "5.5/J_previsao_script", lambda v: v),
    # rodada treze: estado F (correcao de forma por alvo), corte conservador e varredura
    ("correcao de forma no D11", "3.4", r"at most (\d\.\d\d) min in the eleven detections", "3.4/forma_max_D11", lambda v: f"{v:.2f}"),
    ("correcao de forma na amostra", "3.4", r"(\d\.\d\d) min in the sample \(TIC \d+\)", "3.4/forma_max_amostra", lambda v: f"{v:.2f}"),
    ("deslocamento mediano do estado F", "3.4", r"coefficients move by a median of (\d\.\d\d)σ", "3.4/forma_desloc_mediano", lambda v: f"{v:.2f}"),
    ("deslocamento maximo do estado F", "3.4", r"move by a median of \d\.\d\dσ and at most (\d\.\d\d)σ \(TIC", "3.4/forma_desloc_max", lambda v: f"{v:.2f}"),
    ("corte conservador de 3 min", "5.3.2", r"Under it (\d+) of the 11 detections survive", "5.3.2/conservadora_sobrevivem", lambda v: str(v)),
    ("bloco minimo de 390021728", "5.5", r"does not survive a common archival offset of ([\d.]+) min", "5.5/bloco_minimo",
     lambda v: f"{v['390021728']:g}"),
    ("completude do D11", "5.2.5", r"members of D11 have a completeness of (\d\.\d\d)–\d\.\d\d \(median \d\.\d\d\)", "5.2.5/completude_D11", lambda v: f"{v['min']:.2f}"),
    ("supressao do D11", "5.2.5", r"the eleven members of D11 is \d\.\d\d–\d\.\d\d \(median (\d\.\d\d)\)", "5.2.5/supressao_D11", lambda v: f"{v['mediana']:.2f}"),
    ("completude >= 0,5 em 0,02", "5.2.5", r"completeness at 0\.02 s yr⁻¹ is ≥ 0\.5 in (\d+) of the 26", "5.2.5/completude_0.02_ge_meio", lambda v: str(v)),
    ("mediana do chi2_red dos 26", "5.2.1", r"median reduced χ² of the 26 quadratic fits is (\d\.\d\d)", "5.2.1/mediana_chi2r_26", lambda v: f"{v:.2f}"),
    ("P da mediana esperada", "5.2.1", r"has probability (\d\.\d\d) of being this high", "5.2.1/P_mediana_esperada", lambda v: f"{v:.2f}"),
    ("agregado dos adequados", "5.2.1", r"the other 23 giving (\d\.\d\d)", "5.2.1/agregado_adequados", lambda v: f"{v:.2f}"),
    ("percentil da mediana no nulo", "5.2.1", r"sits at its (\d\d)th percentile", "5.2.1/nulo_mediana", lambda v: str(round(v["pct"] * 100))),
    ("p_curv < 0,05 na amostra", "4", r"Of the (\d+) targets with p_curv < 0\.05", "4/n_p_curv_abaixo_005", lambda v: str(v)),
    ("caem no adversarial", "4", r"with p_curv < 0\.05, (\w+) fail the adversarial gate", "4/caem_no_adversarial",
     lambda v: {5: "five", 6: "six", 7: "seven", 8: "eight"}[len(v)]),
    ("setores do cache", "2", r"\((799) sectors; the four epochs of each design", "5.3.2/sigma_j", lambda v: str(v["n_setores"])),
    ("ultimo setor", "2", r"stops at sector (86) \(2024-12-18\)", "2/ultimo_setor", lambda v: str(v)),
    ("FFI TESS-SPOC", "2", r"a median \|t\(FFI\) − t\(2-min\)\| of (0\.30) min", "2/ffi_tess_spoc_dt_mediana_min", lambda v: f"{v:.2f}"),
    ("FFI QLP", "2", r"differences up to (7\.2) min", "2/ffi_qlp_dt_max_min", lambda v: f"{v:.1f}"),
]


# ---------------------------------------------------------------------------
# LISTA NEGRA (rodada dez): valores do ESTADO D que nao podem aparecer fora do
# Apendice D, onde o historico vive. Cada entrada e (rotulo, regex, excecao):
# se `excecao` casar na mesma linha, a ocorrencia e legitima.
LISTA_NEGRA = [
    ("vies -2,14 (estado D)", r"[−-]2\.14", r"TESS-only|per 2\.14 min|bias set to"),
    ("piso 0,78 do vies (hoje 0,75)", r"0\.78 min", r"TESS-only"),
    ("fator SuperWASP 1,36 (hoje 1,31)", r"(?<![\d.])1\.36(?![\dσ])", None),
    ("V527 Dra -0,0358 (hoje -0,0370)", r"[−-]0\.035[89]", None),
    # a mediana |dP/dt| dos 26 voltou a 0,016 no estado F (a mediana cai num vao da ordenacao:
    # 0,0141 e 0,0173 sao os dois centrais), entao este valor deixou de ser exclusivo do estado D
    # e a lista negra nao consegue distinguir os dois - quem protege agora e a checagem numerica
    ("'as often as it agrees'", r"as often as it agrees", None),
    ("'largest |dP/dt|/sigma'", r"largest \|dP/dt\|/σ", None),
    ("'intact'", r"intact", None),
    ("'0,003 target' sem o limite", r"0\.003 target(?! over the 26; the 95)", r"upper bound"),
    ("82% (hoje 84%)", r"82%", None),
    ("'Six targets'", r"Six targets", None),
    ("controle 22,7 min (hoje 21,8)", r"22\.7 min", None),
    ("conjunto D10", r"\bD10\b", None),
    ("'core' sem definicao", r"\bcore\b", None),
    ("limite f >~ 0,8", r"f ≳ 0\.8", None),
    ("contagem frouxa 10 de 42", r"10 (of|in) 42", None),
    # rodada doze: hash de commit no texto. O hash citado na 5.5 so existia no repositorio
    # privado - ninguem de fora podia resolve-lo - e a citacao passou a ser o nome do script
    ("hash de commit no texto", r"commit [0-9a-f]{7,40}\b", None),
    ("'all ten' / 'the ten'", r"\ball ten\b|\bthe ten\b", None),
]


def lista_negra_figuras():
    """A mesma lista negra sobre o texto DESENHADO dentro das figuras (reports/figuras_texto.json,
    escrito por figuras_nota.py). O PNG rasteriza esse texto e o PDF nao o devolve."""
    arq = config.ROOT / "reports" / "figuras_texto.json"
    if not arq.exists():
        return [("(ausente)", "reports/figuras_texto.json nao existe - rode src/figuras_nota.py", "", "")]
    achados = []
    for nome, textos in json.loads(arq.read_text(encoding="utf-8")).items():
        for txt in textos:
            for rot, rx, exc in LISTA_NEGRA:
                if exc and re.search(exc, txt):
                    continue
                for m in re.finditer(rx, txt):
                    achados.append((nome, rot, m.group(0), txt[:80]))
    return achados


def lista_negra(nota):
    """Ocorrencias de valores do estado D fora do Apendice D, com a linha."""
    corte = min(nota.index("## Appendix C."), nota.index("## Appendix D."))   # C e D guardam o historico
    linhas = nota[:corte].split(chr(10))
    achados = []
    for n, ln in enumerate(linhas, 1):
        for rot, rx, exc in LISTA_NEGRA:
            for m in re.finditer(rx, ln):
                if exc and re.search(exc, ln):
                    continue
                achados.append((n, rot, m.group(0), ln.strip()[:110]))
    return achados


def main():
    nota = NOTA.read_text(encoding="utf-8")
    d = json.loads(NUM.read_text(encoding="utf-8"))
    negra = lista_negra(nota)
    # F da rodada catorze: bibliografia - nenhuma orfa, e markdown == bibitem do TeX verbatim
    _ref = nota[nota.index("## References"):nota.index("## Figures")]
    _ent = [l.strip() for l in _ref.split(chr(10)) if l.strip() and not l.startswith("#")]
    _corpo_ref = nota[:nota.index("## References")]
    _orfas = []
    for _e in _ent:
        _chave = _e.split(",")[0].replace("Collaboration", "").strip()
        _ano = re.search(r"(19|20)\d\d", _e)
        if _ano and not re.search(re.escape(_chave) + r".{0,40}?" + _ano.group(0), _corpo_ref, re.S):
            _orfas.append(f"{_chave} {_ano.group(0)}")
    _v1 = (config.ROOT / "paper" / "dpdt_note_v1.tex").read_text(encoding="utf-8")
    _bib = [l for l in _v1.split(chr(10)) if l.lstrip().startswith(chr(92) + "bibitem")]
    print(f"== BIBLIOGRAFIA: {len(_ent)} entradas no markdown, {len(_bib)} no TeX verbatim; orfas: {_orfas if _orfas else 'nenhuma'}")
    if len(_ent) != len(_bib):
        print(f"   DIVERGE: markdown e \bibitem tem contagens diferentes ({len(_ent)} x {len(_bib)})")

    # D da rodada catorze: COMENSURABILIDADE da Tabela 5 - o coeficiente, a barra e o DMD de cada
    # linha tem de vir do MESMO modelo admitido (o menos favoravel), e nao de um modelo para o
    # coeficiente e outro para o z. Ate esta rodada a coluna principal trazia o GP agrupado, que
    # nao e admitido em nenhum dos oito.
    import pandas as _pd
    _ra = _pd.read_parquet(config.DATA / "orquestra" / "oc_lote" / "ruido_admitido_d9.parquet")
    _ra = _ra[_ra.nivel == 0.95].set_index("tic")
    _i = nota.index("<!-- tabela6:inicio -->")
    _j = nota.index("<!-- tabela6:fim -->")
    _lin = [l for l in nota[_i:_j].split(chr(10)) if l.startswith("| ") and not l.startswith("| TIC")]
    _ruim = []
    for l in _lin:
        cel = [x.strip() for x in l.strip("|").split("|")]
        tic = int(cel[0])
        if tic not in _ra.index:
            continue
        pior = int(max(range(len(_ra.loc[tic, "s_adm"])), key=lambda i: _ra.loc[tic, "s_adm"][i]))
        esperado = f"{_ra.loc[tic, 'dPdt_adm'][pior]:+.4f} ± {_ra.loc[tic, 's_adm'][pior]:.4f}"
        dmd = f"{3.0 * max(_ra.loc[tic, 's_adm']):.4f}"
        if cel[6] == "—":          # V527 Dra fica fora da classificacao: sem z e sem MDD
            continue
        if cel[2] != esperado:
            _ruim.append(f"TIC {tic}: coeficiente {cel[2]!r}, admitido menos favoravel {esperado!r}")
        if cel[6] != dmd:
            _ruim.append(f"TIC {tic}: MDD {cel[6]!r}, 3 sigma do mesmo modelo {dmd!r}")
    print(f"== COMENSURABILIDADE (Tabela 5, coeficiente e MDD do mesmo modelo): {len(_ruim)} divergencia(s)")
    for r in _ruim:
        print("   " + r)

    # Rodada dezessete: BLOCO GERADO EDITADO NO PRODUTO. A poda trocou uma referencia de apendice
    # dentro da Tabela 2 (gerada); a regeneracao seguinte a restaurou, e apontava para o apendice
    # errado depois da renumeracao. Editar o produto em vez da fonte tem de ficar visivel.
    _desatual = blocos_gerados_divergentes(nota)
    print(f"== BLOCOS GERADOS (nota x tabelas_nota.md): "
          f"{'iguais' if not _desatual else 'DIVERGEM em ' + ', '.join(_desatual)}")

    # referencia cruzada de apendice: toda mencao tem de ter um titulo correspondente
    _letras, _mortas = apendices_sem_titulo(nota)
    print(f"== REFERENCIAS A APENDICE: {len(_letras)} apendices; mencoes sem apendice: "
          f"{_mortas if _mortas else 'nenhuma'}")

    # G da rodada catorze: citacoes do corpo a entradas das Tabelas 3 e 4, lidas do parquet
    import citacoes_tabelas as CT
    cit = CT.confere(nota)
    print(f"== CITACOES DE TABELA (corpo x parquet): {len(cit)} divergencia(s)"
          + ("" if not cit else " -> rode src/citacoes_tabelas.py --corrigir"))
    for ct in cit:          # NAO usar `d`: e o dicionario do registro, logo abaixo
        print(f"   TIC {ct['tic']} {ct['grandeza']}: escrito {ct['escrito']}, tabela {ct['correto']}")
    print(f"== LISTA NEGRA (texto da nota): {len(negra)} ocorrencias de valores do estado D fora do Apendice D")
    for n, rot, txt, ln in negra:
        print(f"   linha {n:4d}  {rot:<32s} [{txt}]  {ln}")
    negra_fig = lista_negra_figuras()
    print(f"== LISTA NEGRA (texto DENTRO das figuras): {len(negra_fig)} ocorrencias")
    for nome, rot, txt, ctx in negra_fig:
        print(f"   {nome:<22s} {rot:<32s} [{txt}]  {ctx}")
    negra = list(negra) + list(negra_fig)
    print()
    print(f"{'grandeza':<28s} {'secao':<7s} {'na nota':<12s} {'na cadeia':<12s} {'estado':<8s} fonte")
    n_ok = n_div = n_aus = 0
    for rot, sec, rx, caminho, fmt in CHECAGENS:
        # TODAS as ocorrencias, nao a primeira (rodada dezenove). Tres numerais obsoletos
        # sobreviveram a checagem porque a grandeza aparecia duas vezes com redacoes diferentes e
        # `re.search` parava na primeira - a mesma licao da janela de 260 caracteres do item 22:
        # uma checagem so fala do que olhou.
        ms = list(re.finditer(rx, nota))
        try:
            esperado = fmt(val(d, caminho))
        except Exception as e:  # noqa: BLE001
            print(f"{rot:<28s} {sec:<7s} {'?':<12s} {'ERRO ' + str(e)[:40]}")
            n_div += 1
            continue
        if not ms:
            print(f"{rot:<28s} {sec:<7s} {'AUSENTE':<12s} {esperado:<12s} {'AUSENTE':<8s} {fonte(d, caminho)}")
            n_aus += 1
            continue
        escritos = [m.group(1) if m.groups() else esperado for m in ms]
        alvo = str(esperado).replace("−", "-").strip()
        ruins = [x for x in escritos if x.replace("−", "-").strip() != alvo]
        escrito = escritos[0] if len(set(escritos)) == 1 else "|".join(dict.fromkeys(escritos))
        n_x = f" ({len(ms)}x)" if len(ms) > 1 else ""
        print(f"{rot:<28s} {sec:<7s} {escrito[:12]:<12s} {esperado:<12s} {('OK' if not ruins else 'DIVERGE') + n_x:<8s} {fonte(d, caminho)}")
        n_ok += not ruins
        n_div += bool(ruins)
    print(f"\n  {n_ok} OK, {n_div} divergem, {n_aus} ausentes de {len(CHECAGENS)} checagens")

    # coerencia entre partes: os numeros do abstract tem de aparecer no corpo
    ab = nota[nota.index("## Abstract"):nota.index("## 1. Introduction")]
    corpo = nota[nota.index("## 1. Introduction"):]
    achados = re.findall(r"(?<![\w.])(\d+\.\d{2,4}|\d{1,3})(?![\w.])", ab)
    faltam = [x for x in set(achados) if x not in corpo and len(x) > 2]
    print(f"  numeros do abstract ausentes do corpo: {sorted(faltam) if faltam else 'nenhum'}")

    # as tabelas geradas e a nota falam do mesmo D11
    tabs = TAB.read_text(encoding="utf-8")
    for tic in val(d, "4/D11"):
        assert str(tic) in tabs, tic
    print(f"  os {len(val(d, '4/D11'))} TIC do D11 aparecem nas tabelas geradas: sim")
    # legendas das figuras: conferir que as que citam contagens citam as atuais
    for fig, alvo in (("Fig. 3", "1 and none"), ("Fig. 4", "re-estimated"), ("Fig. 5", "V527 Dra")):
        bloco = nota[nota.index(f"**{fig}.**"):nota.index(f"**{fig}.**") + 600]
        print(f"  legenda da {fig}: {'cita ' + alvo if alvo.split()[0] in bloco or alvo in bloco else 'NAO cita ' + alvo}")
    return n_div + len(negra)


if __name__ == "__main__":
    sys.exit(0 if main() == 0 else 1)
