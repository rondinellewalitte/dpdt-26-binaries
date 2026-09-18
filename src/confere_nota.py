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
    ("|dP/dt| mediana", "4", r"The median \|dP/dt\| of (0\.014) s yr⁻¹ over all 26", "4/abs_mediana", lambda v: f"{v:.3f}"),
    ("teste de sinal", "4", r"12/14 \(sign test p = (0\.85)\)", "4/sinais", lambda v: f"{v['p']:.2f}"),
    ("D11", "4", r"\*\*D11\*\* is the eleven targets", "4/D11", lambda v: str(len(v))),
    ("fator SuperWASP", "5.2", r"a variance factor of (1\.31) \(χ² interval", "5.2/fator_SuperWASP", lambda v: f"{v:.2f}"),
    ("fator TESS", "5.2", r"50\.2 against 31\.6, a factor of (1\.59)", "5.2/fator_TESS", lambda v: f"{v:.2f}"),
    ("A50 da completude", "5.2", r"a median of (0\.014) s yr⁻¹", "5.2/completude", lambda v: f"{v['A50_mediana']:.3f}"),
    ("supressao em 0,2", "5.2", r"median factor of (2\.9) at 0\.2 s yr⁻¹", "5.2/completude", lambda v: f"{v['supressao_mediana_0.2']:.1f}"),
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
    ("so-TESS 392536812 (admitido)", "5.5", r"\*\*TIC 392536812\*\* — block \+0\.021 ± (0\.007) under its three admitted models", "5.3.2/chi2r_agrupado_nos_8",
     lambda v: "0.007"),   # valor do modelo admitido (diagonal/A proprio/branco), reconcilia_hi_d9.parquet: diagonal_s = 0,0074
    ("excursao 232634196", "5.5", r"an excursion of (8\.5) min that a cubic term", "5.5/excursao_min", lambda v: f"{v['232634196']:.1f}"),
    ("V564 Dra VarAstro", "5.5", r"−0\.0066 ± 0\.0011 s yr⁻¹ over the whole record against our −0\.0207 ± 0\.0026 \(z = (−5\.1)\)", "5.5/brno",
     lambda v: f"{v['V564_z']:.1f}".replace("-", "−")),
    ("previsao E = 2206", "5.5", r"gives BJD_TDB (2461588\.806)", "5.5/c232_previsao", lambda v: f"{v['quadratica_tudo']['t_bjd']:.3f}"),
    ("separacao linear-quadratica", "5.5", r"the linear one 2461588\.821 \(\+(20\.6) min later\)", "5.5/c232_previsao", lambda v: f"{v['separacao_lin_quad_min']:.1f}"),
    ("P local 2024", "5.5", r"is \+(0\.46) ± 0\.44 s away from the mean period", "5.5/c232_local", lambda v: f"{v['P_menos_escada_s']:.2f}"),
    ("script da previsao J", "5.5", r"when those sectors arrive, `(ffi_fora_amostra_d9)`", "5.5/J_previsao_script", lambda v: v),
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
    ("mediana |dP/dt| 0,016 (hoje 0,014)", r"(?<![\d.])0\.016(?![\d])", None),
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
    corte = min(nota.index("## Appendix D."), nota.index("## Appendix E."))   # D e E guardam o historico
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
        m = re.search(rx, nota)
        try:
            esperado = fmt(val(d, caminho))
        except Exception as e:  # noqa: BLE001
            print(f"{rot:<28s} {sec:<7s} {'?':<12s} {'ERRO ' + str(e)[:40]}")
            n_div += 1
            continue
        if not m:
            print(f"{rot:<28s} {sec:<7s} {'AUSENTE':<12s} {esperado:<12s} {'AUSENTE':<8s} {fonte(d, caminho)}")
            n_aus += 1
            continue
        escrito = m.group(1) if m.groups() else esperado
        ok = escrito.replace("−", "-").strip() == str(esperado).replace("−", "-").strip()
        print(f"{rot:<28s} {sec:<7s} {escrito:<12s} {esperado:<12s} {'OK' if ok else 'DIVERGE':<8s} {fonte(d, caminho)}")
        n_ok += ok
        n_div += not ok
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
