# -*- coding: utf-8 -*-
"""Converte reports/dpdt_note.md em paper/dpdt_note.tex (classe aa) preservando
as decisoes de producao ja tomadas no TeX anterior (paper/dpdt_note_v1.tex):
preambulo, autor/instituto, keywords, as Tabelas 1 e 2, os ambientes das
Tabelas 3-5 (caption, rotulo, cabecalho: duas sidewaystable* e uma table*),
os quatro ambientes de figura, e a bibliografia em \\bibitem autor-ano. Esses
blocos sao copiados VERBATIM do TeX anterior. O CORPO das Tabelas 3, 4 e 5
(cada linha de dados) e gerado aqui a partir das tabelas markdown que
tabela_nota.py escreve entre os marcadores da nota - nenhum numero e
transcrito a mao; os codigos de classe da Tabela 3 (C, --, a, dagger) e o
asterisco da Tabela 5 sao derivados das colunas do markdown. O texto corrido,
os titulos de secao, o abstract estruturado, as captions das figuras, os
agradecimentos e os apendices vem do markdown, convertidos aqui. As
referencias cruzadas "Section X.Y", "Table N", "Fig. N", "Appendix X" viram
\\ref pelos rotulos do TeX anterior.

Uso: python src/md_para_tex.py  (escreve paper/dpdt_note.tex)
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402

MD = config.ROOT / "reports" / "dpdt_note.md"
TEX_V1 = config.ROOT / "paper" / "dpdt_note_v1.tex"
TEX = config.ROOT / "paper" / "dpdt_note.tex"

# ---------- rotulos ----------
SECOES = {  # titulo do markdown -> (comando, rotulo)
    "1. Introduction": ("section", "s:intro"),
    "2. The sample as a vehicle: 26 binaries with data sufficient for the tests": ("section", "s:sel"),
    "3. Method": ("section", "s:met"),
    "3.1 Eclipse epochs": ("subsection", "ss:epocas"),
    "3.2 Period ladder and cycle count": ("subsection", "ss:escada"),
    "3.3 Timing uncertainties from the data": ("subsection", "ss:barras"),
    "3.4 SuperWASP timing bias": ("subsection", "ss:vies"),
    "3.5 Fit and tests": ("subsection", "ss:testes"),
    "4. The chain applied to the 26": ("section", "s:res"),
    "5. Three mechanisms that can produce a period change, and what the data can test": ("section", "s:val"),
    "5.1 Formal timing errors": ("subsection", "ss:formais"),
    "5.2 The residual test and its null: what data-derived bars can be validated against": ("subsection", "ss:qualbarra"),
    "5.2.1 The data-derived bars against the expected reduced χ²": ("subsubsection", "sss:esperado"),
    "5.2.2 The leverage-corrected decomposition by data set, and the composition error it exposed": ("subsubsection", "sss:decomp"),
    "5.2.3 The null of the estimator: bars re-estimated from two or three samples": ("subsubsection", "sss:nulo"),
    "5.2.4 What the residuals leave for orbits, with the same estimator": ("subsubsection", "sss:orbitas"),
    "5.2.5 A secular term through the same estimator: completeness and self-suppression": ("subsubsection", "sss:secular"),
    "5.3 Between-sector timing noise: the secondaries and the TESS sensitivity ceiling": ("subsection", "ss:setores"),
    "5.3.1 Secondary minima: the differential drift between sectors": ("subsubsection", "sss:secundarios"),
    "5.3.2 TESS bars at the between-sector scale: a sensitivity test with a ceiling": ("subsubsection", "sss:teto"),
    "5.4 Undersampled third-body orbits": ("subsection", "ss:orb"),
    "5.5 External records: published work, archival minima, and the one detection carried by the lever": ("subsection", "ss:externo"),
    "6. What is and is not claimed": ("section", "s:claim"),
    "Data and code": ("section", "s:data"),
    "Data availability": ("section", "s:avail"),
    "Appendix A. Outcome of the recovery pass, per target": ("appendix", "a:A"),
    "Appendix B. Control target and period-ladder audit": ("appendix", "a:B"),
    "Appendix C. The orbit population, the efficiency, and the occurrence numbers": ("appendix", "a:C"),
    "Appendix D. Record of readings withdrawn between versions": ("appendix", "a:D"),
}
REF_SECAO = {"3.1": "ss:epocas", "3.2": "ss:escada", "3.3": "ss:barras", "3.4": "ss:vies", "3.5": "ss:testes",
             "5.1": "ss:formais", "5.2": "ss:qualbarra", "5.3": "ss:setores", "5.4": "ss:orb", "5.5": "ss:externo",
             "5.2.1": "sss:esperado", "5.2.2": "sss:decomp", "5.2.3": "sss:nulo", "5.2.4": "sss:orbitas", "5.2.5": "sss:secular",
             "5.3.1": "sss:secundarios", "5.3.2": "sss:teto",
             "1": "s:intro", "2": "s:sel", "3": "s:met", "4": "s:res", "5": "s:val", "6": "s:claim"}
REF_TABELA = {"1": "t:funil", "2": "t:naomedidos", "3": "t:coef", "4": "t:dist", "5": "t:barras"}
REF_FIG = {"1": "f:dpdt", "2": "f:cvdra", "3": "f:barras", "4": "f:compl", "5": "f:v527"}   # CV Dra abre a Secao 5 (revisao externa, item 6)
REF_APX = {"A": "a:A", "B": "a:B", "C": "a:C", "D": "a:D"}


# ---------- conversao inline ----------
def esc(t):
    """escapa & % # _ fora de codigo e de matematica (o texto corrido nao tem $)."""
    return t.replace("\\", "\\textbackslash{}").replace("&", "\\&").replace("%", "\\%").replace("#", "\\#").replace("_", "\\_")


SIMBOLOS = [
    ("Σχ²", "$\\Sigma\\chi^2$"), ("Q·E²", "$Q\\cdot E^2$"), ("QE²", "$QE^2$"), ("0.78²", "$0.78^2$"), ("a₃", "$a_3$"), ("10⁻²¹", "$10^{-21}$"), ("10⁻¹⁴", "$10^{-14}$"),
    ("10⁻¹⁹", "$10^{-19}$"), ("10⁻¹²", "$10^{-12}$"), ("std(d)/√n", "${\\rm std}(d)/\\sqrt{n}$"), ("std(d)", "${\\rm std}(d)$"),
    ("√2–√3", "$\\sqrt{2}$--$\\sqrt{3}$"), ("|t_a − t_b|/√2", "$|t_a - t_b|/\\sqrt{2}$"), ("|d₁ − d₂|/√2", "$|d_1 - d_2|/\\sqrt{2}$"), ("|Δ|/√2", "$|\\Delta|/\\sqrt{2}$"), ("±Δ/2", "$\\pm\\Delta/2$"), ("difference Δ", "difference $\\Delta$"),
    ("±(d₁ − d₂)/2", "$\\pm(d_1 - d_2)/2$"), ("z² = r²/s²", "$z^2 = r^2/s^2$"), ("z² = 0.5", "$z^2 = 0.5$"), ("√n", "$\\sqrt{n}$"),
    ("F(1, ν)", "$F(1,\\nu)$"), ("ν/(ν − 2)", "$\\nu/(\\nu-2)$"), ("ν ≤ 2", "$\\nu \\le 2$"), ("ν = 1–2", "$\\nu$ = 1--2"),
    ("10⁻³⁰", "$10^{-30}$"), ("10⁻⁶", "$10^{-6}$"), ("10⁻³", "$10^{-3}$"), ("10⁻⁴", "$10^{-4}$"), ("10⁻⁷", "$10^{-7}$"),
    ("χ²_red", "$\\chi^2_{\\rm red}$"), ("Δχ²", "$\\Delta\\chi^2$"), ("ΔG", "$\\Delta G$"), ("χ²/ν", "$\\chi^2/\\nu$"), ("χ²", "$\\chi^2$"),
    ("Σε₃", "$\\Sigma\\varepsilon_3$"), ("Σε", "$\\Sigma\\varepsilon$"), ("ε₃", "$\\varepsilon_3$"), ("ε", "$\\varepsilon$"),
    ("|dP/dt|/σ", "$|\\mathrm{d}P/\\mathrm{d}t|/\\sigma$"), ("|dP/dt|", "$|\\mathrm{d}P/\\mathrm{d}t|$"),
    ("σ(dP/dt)", "$\\sigma(\\mathrm{d}P/\\mathrm{d}t)$"), ("dP/dt", "$\\mathrm{d}P/\\mathrm{d}t$"),
    ("σ_P × N", "$\\sigma_P \\times N$"), ("P₃", "$P_3$"), ("M₃", "$M_3$"), ("M_bin", "$M_{\\rm bin}$"), ("M☉", "$M_\\odot$"),
    ("p_gof", "$p_{\\rm gof}$"), ("p_curv", "$p_{\\rm curv}$"), ("p_adv", "$p_{\\rm adv}$"), ("T14", "$T_{14}$"),
    ("Σz²", "$\\Sigma z^2$"), ("Σ(1−h)", "$\\Sigma(1-h)$"), ("Σ (1 − h)", "$\\Sigma(1-h)$"), ("z²/(1−h)", "$z^2/(1-h)$"), ("E[r²/σ²] = 1 − h", "$E[r^2/\\sigma^2]=1-h$"),
    ("f × Σε", "$f\\times\\Sigma\\varepsilon$"), ("√(2.73 − 0.37)", "$\\sqrt{2.73-0.37}$"), ("√1.92", "$\\sqrt{1.92}$"),
    ("O−C", "O$-$C"), ("β", "$\\beta$"), ("σ", "$\\sigma$"), ("ν", "$\\nu$"), ("ρ", "$\\rho$"), ("τ", "$\\tau$"),
    ("≥", "$\\ge$"), ("≤", "$\\le$"), ("≈", "$\\approx$"), ("≳", "$\\gtrsim$"), ("≲", "$\\lesssim$"), ("≫", "$\\gg$"), ("≠", "$\\ne$"),
    ("×", "$\\times$"), ("±", "$\\pm$"), ("→", "$\\rightarrow$"), ("⁻¹", "$^{-1}$"), ("⁻²¹", "$^{-21}$"), ("⁻¹⁴", "$^{-14}$"),
    ("⁻³⁰", "$^{-30}$"), ("⁻⁶", "$^{-6}$"), ("⁻³", "$^{-3}$"), ("⁻⁴", "$^{-4}$"), ("⁻⁷", "$^{-7}$"),
    ("″", "\\arcsec"), ("′", "\\arcmin"), ("°", "$^\\circ$"), ("—", "---"), ("–", "--"), ("“", "``"), ("”", "''"), ("’", "'"),
    ("−", "$-$"), ("…", "\\ldots"), ("½", "$\\frac{1}{2}$"),
]


def inline(t):
    """markdown inline -> TeX; protege codigo e matematica antes de escapar."""
    partes = re.split(r"(`[^`]*`)", t)
    out = []
    for k, parte in enumerate(partes):
        if k % 2 == 1:                       # codigo
            # nomes de arquivo longos em \texttt nao quebram: ponto de quebra opcional depois de "src/"
            out.append("\\texttt{" + esc(parte[1:-1]).replace("src/", "src/\\allowbreak ") + "}")
            continue
        # simbolos primeiro (com placeholders, para o escape nao mexer neles), depois o escape
        s = parte
        subst = {}
        for n, (a, b) in enumerate(SIMBOLOS):
            if a in s:
                # terminador DIFERENTE do inicio: com "\x00N\x00", o fim de um placeholder mais os digitos do texto
                # seguinte mais o inicio do proximo formavam outro placeholder ("40–60°" virou "40\x0082\x0060\x00..." e
                # o "\x0060\x00" foi lido como o simbolo 60)
                chave = f"\x00{n}\x01"
                s = s.replace(a, chave)
                subst[chave] = b
        s = esc(s)
        for chave, b in subst.items():
            s = s.replace(chave, b)
        # numeros com espaco de milhar: "12 497" -> "12\,497"
        s = re.sub(r"(?<=\d) (?=\d{3}\b)", "\\\\,", s)
        # negrito e italico do markdown
        s = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", s)
        s = re.sub(r"(?<![\w\\])\*(?!\*)(.+?)(?<!\*)\*(?![\w])", r"\\textit{\1}", s)
        # aspas ASCII do markdown -> aspas TeX
        s = re.sub(r'"([^"\n]+)"', r"``\1''", s)
        # "~" do markdown (aprox.) -> $\sim$ quando precede numero
        s = re.sub(r"~(?=\d)", "$\\\\sim$", s)
        # no markdown as colunas das inflacoes estao na Tabela 3; no TeX elas sao a Tabela 5
        s = s.replace("The last row of Table 4 and the last two columns of Table 3 give", "The last row of Table 4 and Table 5 give")
        s = s.replace("(Table 3, last two columns; Table 4; Fig. 3)", "(Table 5; Table 4; Fig. 3)")
        # referencias cruzadas; "Table N of <Autor>" e tabela de outro artigo, nao vira \ref (pego pela guarda de
        # paragrafos: 'Table 5 of Tokovinin' saia como \ref{t:barras}, a NOSSA Tabela 5)
        externas = re.findall(r"Table \d of [A-Z]\w+", s)
        for n, ext in enumerate(externas):
            s = s.replace(ext, f"\x02{n}\x02")
        s = re.sub(r"Sections? (\d(?:\.\d){0,2}) and (\d(?:\.\d){0,2})", lambda m: f"Sects.~\\ref{{{REF_SECAO[m.group(1)]}}} and \\ref{{{REF_SECAO[m.group(2)]}}}", s)
        s = re.sub(r"Section (\d(?:\.\d){0,2})", lambda m: f"Sect.~\\ref{{{REF_SECAO[m.group(1)]}}}", s)
        s = re.sub(r"Tables (\d)--(\d)", lambda m: f"Tables~\\ref{{{REF_TABELA[m.group(1)]}}}--\\ref{{{REF_TABELA[m.group(2)]}}}", s)
        s = re.sub(r"Table (\d)", lambda m: f"Table~\\ref{{{REF_TABELA[m.group(1)]}}}", s)
        s = re.sub(r"Figs\. (\d), (\d) and (\d)", lambda m: "Figs.~" + ", ".join(f"\\ref{{{REF_FIG[g]}}}" for g in m.groups()[:2]) + f" and \\ref{{{REF_FIG[m.group(3)]}}}", s)
        s = re.sub(r"Figs\. (\d) and (\d)", lambda m: f"Figs.~\\ref{{{REF_FIG[m.group(1)]}}} and \\ref{{{REF_FIG[m.group(2)]}}}", s)
        s = re.sub(r"Fig\. (\d)", lambda m: f"Fig.~\\ref{{{REF_FIG[m.group(1)]}}}", s)
        s = re.sub(r"Appendix ([ABC])", lambda m: f"Appendix~\\ref{{{REF_APX[m.group(1)]}}}", s)
        for n, ext in enumerate(externas):
            s = s.replace(f"\x02{n}\x02", ext)
        # unidades
        s = s.replace(" s yr$^{-1}$", " s\\,yr$^{-1}$").replace(" min yr$^{-1}$", " min\\,yr$^{-1}$").replace(" d yr$^{-1}$", " d\\,yr$^{-1}$")
        # matematica adjacente: "$\ge$ 2" fica; junta "$a$$b$" que a substituicao possa ter criado
        s = s.replace("$$", "")
        out.append(s)
    return "".join(out)


# ---------- blocos verbatim do TeX anterior ----------
def bloco_tex(v1, ini, fim):
    i = v1.index(ini)
    j = v1.index(fim, i) + len(fim)
    return v1[i:j]


# ---------- corpo das Tabelas 3, 4 e 5, gerado das tabelas markdown de tabela_nota.py ----------
def linhas_md(md, marcador):
    """Linhas de uma tabela markdown entre <!-- marcador:inicio --> e :fim, como listas de celulas
    (cabecalho incluido; a linha de tracos nao). '\\|' dentro de celula e um pipe literal."""
    bloco = md[md.index(f"<!-- {marcador}:inicio -->"):md.index(f"<!-- {marcador}:fim -->")]
    out = []
    for l in bloco.splitlines():
        if not l.startswith("|") or re.match(r"\|\s*-", l):
            continue
        out.append([c.strip().replace("\\|", "|") for c in re.split(r"(?<!\\)\|", l.strip())[1:-1]])
    return out


def p_tex(v):
    """'< 10⁻³' -> $<10^{-3}$; os demais p ficam como estao; o '*' de 'ainda detectado' (posto por tabela_nota.py) e preservado."""
    estrela = "*" if v.endswith("*") else ""
    v = v.rstrip("*")
    return ("$<10^{-3}$" if v.startswith("<") else v) + estrela


def nome_tex(n):
    return esc(n).replace("UCAC4 ", "UCAC4\\,")


def codigo_classe(c):
    """coluna 'class' do markdown -> codigo da Tabela 3: C curvatura, -- nenhuma; a falha no adversarial; dagger parabola inadequada."""
    assert c.startswith("curvature") or c.startswith("no significant curvature"), c
    cod = "C" if c.startswith("curvature") else "--"
    if "fails adversarial" in c:
        cod += "a"
    if "parabola inadequate" in c:
        cod += "\\,$\\dagger$"
    return cod


def corpo_tabela3(md):
    linhas = linhas_md(md, "tabela3")
    cab = linhas[0]
    assert cab[:14] == ["TIC", "Name", "RA, Dec (J2000, deg)", "Tmag", "P (d)", "N (TESS+SW)", "d.o.f.", "span (yr)", "quadratic coefficient as dP/dt (s yr⁻¹)",
                        "χ²_red", "p_gof", "p_curv", "p_adv", "class"], cab[:14]
    t3, t5 = [], []
    for c in linhas[1:]:
        ra, dec = c[2].split(", ")
        t3.append(" & ".join([c[0], nome_tex(c[1]), ra, dec, c[3], c[4], c[5], c[6], c[7], c[8].replace("±", "$\\pm$"), c[9],
                              p_tex(c[10]), p_tex(c[11]), p_tex(c[12]), codigo_classe(c[13])]) + "\\\\")
        celulas = [c[0], nome_tex(c[1])]
        for par in c[14:16]:                      # "p_curv, p_adv*" por configuracao (TESS +1,0 / +2,6); * = ainda detectado
            pc, pa = [x.strip() for x in par.split(",")]
            celulas += [p_tex(pc), p_tex(pa)]
        t5.append(" & ".join(celulas) + "\\\\")
    assert len(t3) == 26, len(t3)
    return "\n".join(t3) + "\n", "\n".join(t5) + "\n"


def tabela4(md, caption_label):
    linhas = linhas_md(md, "tabela4")
    corpo = "\n".join(" & ".join(inline(c) for c in l) + "\\\\" for l in linhas[1:]) + "\n"
    cab = " & ".join(inline(c) for c in linhas[0]) + "\\\\\n"
    return (caption_label + "\\centering\n\\small\n\\begin{tabular}{p{0.30\\linewidth} p{0.19\\linewidth} p{0.20\\linewidth} p{0.20\\linewidth}}\n"
            "\\hline\\hline\n" + cab + "\\hline\n" + corpo + "\\hline\n\\end{tabular}\n\\end{table*}")


def trocar_corpo(env, corpo):
    """substitui as linhas de dados de um ambiente tabular (entre o \\hline do cabecalho e o \\hline final)."""
    a = env.index("\\hline\n", env.index("\\hline\\hline\n") + len("\\hline\\hline\n")) + len("\\hline\n")
    b = env.rindex("\\hline\n\\end{tabular}")
    return env[:a] + corpo + env[b:]


def main():
    md = MD.read_text(encoding="utf-8")
    v1 = TEX_V1.read_text(encoding="utf-8")
    preambulo = v1[:v1.index("\\begin{document}")]
    preambulo = re.sub(r"^%%.*\n(%%.*\n)*", "", preambulo)  # cabecalho de comentarios antigo
    # Tabela 4 (343 pt) + Fig. 2 (~290 pt) somam 0,90 da altura de texto (705 pt): abaixo do 0,95 que o aa.cls exige
    # para uma pagina so de floats, e acima do que cabe no topo com texto embaixo; sem isto a Fig. 2 vai da p. 6 para a p. 8
    preambulo += "\\renewcommand*\\dblfloatpagefraction{0.85}\n"
    cabecalho = ("%% How archival eclipse timing can produce period change\n%% A&A regular article. Requires aa.cls (in this folder).\n"
                 "%% Figures expected in ./figures/ : fig1_barras.png fig4_dpdt.png fig2_v527_injecao.png fig3_cvdra.png fig5_completude.png\n"
                 "%% Generated from reports/dpdt_note.md by src/md_para_tex.py. Table bodies (Tables 3-5) come from the markdown\n"
                 "%% tables written by src/tabela_nota.py; preamble, captions, figure environments and bibliography are carried\n"
                 "%% verbatim from dpdt_note_v1.tex. Compile: pdflatex dpdt_note ; pdflatex dpdt_note\n")
    autor = bloco_tex(v1, "\\author{", "\\date{Received ---; accepted ---}")
    keywords = bloco_tex(v1, "\\keywords{", "}")
    tabelas = {
        "t:funil": bloco_tex(v1, "\\begin{table}\n\\caption{Selection funnel.}", "\\end{table}"),
        "t:naomedidos": bloco_tex(v1, "\\begin{table*}\n\\caption{\\textbf{The 62 targets", "\\end{table*}"),
        "t:coef": bloco_tex(v1, "\\begin{sidewaystable*}\n\\caption{Quadratic ephemerides", "\\end{sidewaystable*}"),
        "t:dist": bloco_tex(v1, "\\begin{table*}\n\\caption{The distribution, in three partitions", "\\end{table*}"),
    }
    i = v1.index("\\end{table*}", v1.index("\\label{t:dist}")) + len("\\end{table*}")
    tabelas["t:barras"] = bloco_tex(v1[i:], "\\begin{sidewaystable*}", "\\end{sidewaystable*}")
    figuras = {}
    for rot, arq in (("f:barras", "fig1_barras"), ("f:dpdt", "fig4_dpdt"), ("f:v527", "fig2_v527_injecao"), ("f:cvdra", "fig3_cvdra")):
        j = v1.index("{figures/" + arq)               # o ambiente, nao a mencao no cabecalho
        a = v1.rfind("\\begin{figure", 0, j)
        b = v1.index("\\end{figure", j)
        b = v1.index("}", b) + 1
        figuras[rot] = v1[a:b]
    figuras["f:compl"] = ("\\begin{figure*}\n\\centering\n\\includegraphics[width=\\linewidth]{figures/fig5_completude.png}\n"
                         "\\caption{PLACEHOLDER}\n\\label{f:compl}\n\\end{figure*}")
    bib = bloco_tex(v1, "\\begin{thebibliography}{99}", "\\end{thebibliography}")
    # rotulos de secao antigos dentro das tabelas/figuras: s:orb (antiga secao 6) -> ss:orb
    for d in (tabelas, figuras):
        for k in d:
            d[k] = (d[k].replace("\\ref{s:orb}", "\\ref{ss:orb}").replace("(Section 6)", "(Sect.~\\ref{ss:orb})")
                    .replace("(Section 5, sensitivity ceiling)", "(Sect.~\\ref{ss:setores}, sensitivity ceiling)")
                    .replace("(Section 5)", "(Sect.~\\ref{s:val})"))
    # Tabela 3 (caption verbatim do v1): o O-C publicado de V527 Dra esta em 5.4; e o aviso de que nenhuma entrada e taxa
    chave = "V527 Dra (TIC 424461577) passes all three tests and has a published cyclic O$-$C (Sect.~\\ref{s:val})."
    assert tabelas["t:coef"].count(chave) == 1, chave
    tabelas["t:coef"] = tabelas["t:coef"].replace(chave,
        "No entry of this table is a period-change rate: each is a quadratic coefficient over 16--20\\,yr, and Sect.~\\ref{s:val} shows what such a "
        "coefficient can be made of; the table is not to be cited as measured d$P$/d$t$ (Sect.~\\ref{s:claim}). "
        "V527 Dra (TIC 424461577) passes all three tests and has a published cyclic O$-$C (Sect.~\\ref{ss:orb}).")
    chave = "d$P$/d$t$ in s\\,yr$^{-1}$."
    assert tabelas["t:coef"].count(chave) == 1, chave
    tabelas["t:coef"] = tabelas["t:coef"].replace(chave, "d$P$/d$t$ the quadratic coefficient expressed as a rate, in s\\,yr$^{-1}$.")
    # Tabela 4 (caption verbatim do v1): as particoes do estado B
    for a, b in (("the 22 whose quadratic fit is adequate", "the 23 whose quadratic fit is adequate"), ("and the 15 that are adequate", "and the 16 that are adequate")):
        assert tabelas["t:dist"].count(a) == 1, a
        tabelas["t:dist"] = tabelas["t:dist"].replace(a, b)
    # Tabela 5: gerada inteira aqui (so o teto TESS; a coluna SuperWASP x1,65 saiu em 2026-09-13 - Sect. 5.2.3)
    corpo3, corpo5 = corpo_tabela3(md)
    tabelas["t:coef"] = trocar_corpo(tabelas["t:coef"], corpo3)
    tabelas["t:barras"] = ("\\begin{table*}\n\\caption{Curvature with the TESS bars inflated in quadrature by 1.0 and by 2.6 min\\,yr$^{-1}$ times the "
                           "target lever, the sensitivity ceiling of Sect.~\\ref{sss:teto}, as $p_{\\rm curv}$, $p_{\\rm adv}$; the SuperWASP bars are those of "
                           "Sect.~\\ref{ss:barras}, uncorrected for the reason given in Sect.~\\ref{sss:nulo}. A detection requires both $p<0.05$ (marked $*$).}\n"
                           "\\label{t:barras}\n\\centering\n\\tiny\n\\begin{tabular}{l l r r r r}\n\\hline\\hline\n"
                           "TIC & Name & \\multicolumn{2}{c}{TESS $+$1.0 min\\,yr$^{-1}$} & \\multicolumn{2}{c}{TESS $+$2.6 min\\,yr$^{-1}$}\\\\\n"
                           " &  & $p_{\\rm curv}$ & $p_{\\rm adv}$ & $p_{\\rm curv}$ & $p_{\\rm adv}$\\\\\n\\hline\n"
                           + corpo5 + "\\hline\n\\end{tabular}\n\\end{table*}")
    tabelas["t:dist"] = tabela4(md, tabelas["t:dist"][:tabelas["t:dist"].index("\\centering")])

    # ---------- abstract estruturado ----------
    ab = md[md.index("## Abstract") + len("## Abstract"):md.index("## 1. Introduction")]
    partes = {}
    for chave in ("Context", "Aims", "Methods", "Results", "Conclusions"):
        m = re.search(r"\*" + chave + r"\.\*\s*(.+?)(?=\n\s*\n\*|\Z)", ab, flags=re.S)
        partes[chave] = inline(m.group(1).strip())
    abstract = "\\abstract{" + "}\n{".join(partes[k] for k in ("Context", "Aims", "Methods", "Results", "Conclusions")) + "}\n"

    # ---------- corpo ----------
    corpo = md[md.index("## 1. Introduction"):md.index("## References")]
    corpo = corpo[:corpo.index("## Acknowledgements")]          # agradecimentos vao no lugar da classe
    agradecimentos = md[md.index("## Acknowledgements") + len("## Acknowledgements"):md.index("## References")].strip()
    figs_md = md[md.index("## Figures"):md.index("## Appendix A.")]
    # as CAPTIONS das figuras vem do markdown (a fonte), nao do TeX anterior: so o
    # ambiente (posicao, largura, arquivo, rotulo) e reaproveitado
    for m in re.finditer(r"\*\*Fig\. (\d)\.\*\* (.+?)\s*\(`reports/figures/[^`]+`\)", figs_md):
        rot = REF_FIG[m.group(1)]
        legenda = inline(m.group(2).strip())
        figuras[rot] = re.sub(r"\\caption\{.*?\}\n\\label", lambda _: "\\caption{" + legenda + "}\n\\label", figuras[rot], count=1, flags=re.S)
        assert legenda[:30] in figuras[rot], f"caption da {rot} nao entrou"
    apendices = md[md.index("## Appendix A."):]
    saida = []
    tabela_pendente = None

    def emitir(secao_md):
        nonlocal saida
        linhas = secao_md.split("\n")
        par = []
        em_tabela = False

        lista = None                                              # "itemize" / "enumerate" aberta, ou None

        def fechar_lista():
            nonlocal lista
            if lista:
                saida.append("\\end{" + lista + "}\n\n")
                lista = None

        def flush():
            nonlocal par
            if par:
                texto = " ".join(x.strip() for x in par).strip()
                if texto:
                    # linha em branco depois de cada paragrafo: com "\n" simples o TeX EMENDA os paragrafos
                    # consecutivos num so (foi assim em todas as versoes ate 2026-09-13: cada secao saia como um paragrafo)
                    saida.append(inline(texto) + "\n\n")
                par = []

        for ln in linhas:
            if ln.startswith("<!-- tabela3:inicio"):
                flush(); saida.append(figuras["f:dpdt"] + "\n" + tabelas["t:coef"] + "\n"); em_tabela = True; continue   # Fig. 2 na secao 4
            if ln.startswith("<!-- tabela4:inicio"):
                # a Tabela 4 (distribuicao) e depois a Tabela 5 (as tres inflacoes, colunas que no
                # markdown estao na Tabela 3): numeracao 3, 4, 5 como no TeX anterior
                flush(); saida.append(tabelas["t:dist"] + "\n" + tabelas["t:barras"] + "\n"); em_tabela = True; continue
            if ln.startswith("<!-- tabela3:fim") or ln.startswith("<!-- tabela4:fim"):
                em_tabela = False; continue
            if em_tabela == "md":                                # tabela markdown: pula as linhas "|" (e as vazias) e sai na primeira que nao e
                if ln.startswith("|") or not ln.strip():
                    continue
                em_tabela = False
            if em_tabela or ln.startswith("<!--"):
                continue
            if ln.startswith("**Table 1. Selection funnel.**"):
                flush(); saida.append(tabelas["t:funil"] + "\n"); em_tabela = "md"; continue
            if ln.startswith("**Table 2. "):
                flush(); saida.append(tabelas["t:naomedidos"] + "\n"); em_tabela = "md"; continue
            if ln.startswith("**Table 3. ") or ln.startswith("**Table 4. "):
                continue                                          # as captions ja estao nos blocos verbatim
            if ln.startswith("    chain bias ="):
                flush(); saida.append("\\begin{equation}\n\\text{chain bias} = -2.14 \\pm 0.78\\ \\mathrm{min}.\n\\end{equation}\n"); continue
            m_item = re.match(r"(- |\d+\. )(.*)", ln)
            if m_item:                                            # listas: "- " -> itemize, "1. " -> enumerate
                flush()
                tipo = "itemize" if m_item.group(1) == "- " else "enumerate"
                if lista != tipo:
                    fechar_lista(); saida.append("\\begin{" + tipo + "}\n"); lista = tipo
                saida.append("\\item " + inline(m_item.group(2).strip()) + "\n")
                continue
            if not ln.strip():
                flush(); fechar_lista(); continue
            par.append(ln)
        flush(); fechar_lista()

    for m in re.finditer(r"^(##+) (.+)$", corpo, flags=re.M):
        pass
    partes_sec = re.split(r"^(##+ .+)$", corpo, flags=re.M)
    # partes_sec: ['', '## 1. Introduction', texto, '## 2...', texto, ...]
    for k in range(1, len(partes_sec), 2):
        titulo = partes_sec[k].lstrip("# ").strip()
        texto = partes_sec[k + 1]
        cmd, rot = SECOES[titulo]
        nome = re.sub(r"^\d(\.\d){0,2}\.? ", "", titulo)
        # matematica num titulo vai para o bookmark do PDF via \texorpdfstring (senao o hyperref avisa e descarta)
        nome_tex = re.sub(r"\$([^$]+)\$", lambda m: "\\texorpdfstring{$" + m.group(1) + "$}{" + m.group(1).replace("\\chi^2", "chi2") + "}", inline(nome))
        saida.append(f"\n\\{cmd}{{{nome_tex}}}\\label{{{rot}}}\n")
        if rot == "s:val":                                        # Fig. 2 (CV Dra) na ABERTURA da Secao 5: a demonstracao que enquadra os tres mecanismos
            saida.append(figuras["f:cvdra"] + "\n")
        if rot == "ss:qualbarra":                                 # Fig. 3 (eficiencia x barra) no INICIO da decomposicao, para flutuar cedo
            saida.append(figuras["f:barras"] + "\n")
        if rot in ("s:data", "s:avail"):        # caminhos longos em \texttt: sem overfull
            saida.append("\\begin{sloppypar}\n"); emitir(texto); saida.append("\\end{sloppypar}\n")
        else:
            emitir(texto)
        if rot == "ss:orb":
            saida.append(figuras["f:v527"] + "\n")
        if rot == "sss:secular":                                  # Fig. 5 (completude e supressao) na propria 5.2.5
            saida.append(figuras["f:compl"] + "\n")
    # apendices
    saida.append("\n\\appendix\n")
    partes_apx = re.split(r"^(## Appendix .+)$", apendices, flags=re.M)
    for k in range(1, len(partes_apx), 2):
        titulo = partes_apx[k].lstrip("# ").strip()
        cmd, rot = SECOES[titulo]
        nome = re.sub(r"^Appendix [A-Z]\. ", "", titulo)
        assert not nome.startswith("Appendix "), f"titulo de apendice nao despido (regex das letras desatualizada): {titulo!r}"
        saida.append(f"\n\\section{{{inline(nome)}}}\\label{{{rot}}}\n")
        emitir(partes_apx[k + 1])
    saida.append("\n\\begin{acknowledgements}\n" + inline(agradecimentos) + "\n\\end{acknowledgements}\n\n" + bib + "\n\n\\end{document}\n")

    titulo_md = md.split("\n", 1)[0].lstrip("# ").strip()
    t1, t2 = titulo_md.split(": ", 1)
    doc = (cabecalho + preambulo + "\\begin{document}\n\n"
           + f"\\title{{{inline(t1)}}}\n\\subtitle{{{inline(t2[0].upper() + t2[1:])}}}\n\n" + autor + "\n\n" + abstract + "\n" + keywords + "\n\n\\maketitle\n"
           + "".join(saida))
    TEX.write_text(doc, encoding="utf-8")
    print(f"-> {TEX} ({len(doc.split())} palavras); secoes: {len(partes_sec) // 2}; apendices: {len(partes_apx) // 2}")
    # sobras de markdown ou de unicode que a conversao nao cobriu
    # placeholder que escapou (um "\x00N\x01" nao restaurado) e um byte de controle no TeX: falha alta
    escapou = re.findall(r"[\x00-\x08\x0b\x0c\x0e-\x1f]", doc)
    assert not escapou, f"{len(escapou)} bytes de controle no TeX: placeholder da conversao inline nao restaurado"
    sobras = sorted(set(re.findall(r"[^\x00-\x7F]", doc)) - set("áàâãéêíóôõúçÁÉÍÓÚñüöäšřČŞıŠŁłČěřž°"))
    print("caracteres nao-ASCII restantes no TeX (fora de acentos):", "".join(sobras))
    # o que sobra e o que pdflatex ja aceitou (² via textcomp; ć ő na bibliografia); qualquer outro caractere e falha alta ANTES do pdflatex
    assert not (set(sobras) - set("²ćő")), f"caracteres novos fora da tabela de simbolos: {set(sobras) - set(chr(0xb2) + chr(0x107) + chr(0x151))}"
    print("marcadores markdown restantes: **:", doc.count("**"), "| ` :", doc.count("`") - doc.count("``"))
    # AUSENCIA, nao so presenca: nada de tabela markdown crua no corpo, e nenhum trecho
    # distintivo aparece mais de uma vez (a Secao 2 saiu duplicada uma vez por isso)
    corpo_sem_tab = re.sub(r"\\begin\{(sideways)?table\*?\}.*?\\end\{(sideways)?table\*?\}", "", doc, flags=re.S)
    pipes = [l for l in corpo_sem_tab.split("\n") if l.lstrip().startswith("|")]
    assert not pipes, f"tabela markdown crua no corpo: {pipes[0][:60]!r}"
    for frase in ("Union of the TESS eclipsing-binary catalogues", "SuperWASP curve with $<$ 500 valid points",
                  "The bias was measured rather than left as a caveat", "None of the three mechanisms is unknown",
                  "Consequence of the formal-error run"):
        n = doc.count(frase)
        assert n == 1, f"{frase!r} aparece {n}x no TeX (esperado 1)"
    print("ausencia conferida: 0 linhas de tabela crua fora de ambientes; 5 trechos distintivos aparecem exatamente 1x")
    # paragrafos: no corpo (fora de ambientes e da bibliografia), duas linhas de texto consecutivas sem linha em
    # branco entre elas sao dois paragrafos do markdown que o TeX vai emendar num so
    corpo_texto = re.sub(r"\\begin\{(thebibliography|figure\*?|abstract|acknowledgements)\}.*?\\end\{\1\}", "", corpo_sem_tab, flags=re.S)
    corpo_texto = corpo_texto[corpo_texto.index("\\maketitle"):]
    linhas_tex = corpo_texto.split("\n")

    def e_texto(l):   # linha de paragrafo: nao comeca por comando estrutural (os run-ins comecam por \textbf e SAO texto)
        l = l.lstrip()
        return bool(l) and (not l.startswith(("\\", "%", "&", "TIC ", "Quantity")) or l.startswith("\\text"))
    emendas = [(a[:50], b[:50]) for a, b in zip(linhas_tex, linhas_tex[1:]) if e_texto(a) and e_texto(b)]
    assert not emendas, f"{len(emendas)} pares de paragrafos sem linha em branco entre eles, p.ex. {emendas[0]}"
    # identidade de contagem: os paragrafos do markdown que o conversor emite (corpo + apendices, sem tabelas,
    # titulos e marcadores) tem de ser exatamente os blocos de texto do TeX (fora dos ambientes e da bibliografia)
    n_par_md = sum(1 for p in re.split(r"\n\s*\n", corpo + "\n\n" + apendices)
                   if p.strip() and not p.lstrip().startswith(("|", "#", "<!--", "**Table ", "- ")) and not re.match(r"\d+\. ", p.lstrip()))   # captions vem do v1; listas viram \item
    n_par_tex = sum(1 for l in linhas_tex if e_texto(l))
    assert n_par_md == n_par_tex, f"{n_par_md} paragrafos no markdown, {n_par_tex} blocos de texto no TeX"
    print(f"paragrafos: {n_par_tex} blocos de texto no TeX = {n_par_md} paragrafos no markdown; 0 pares emendados")


if __name__ == "__main__":
    main()
