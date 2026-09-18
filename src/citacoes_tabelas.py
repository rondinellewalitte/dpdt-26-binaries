# -*- coding: utf-8 -*-
"""Rodada catorze, bloco G: a guarda que faltava - todo numeral do CORPO que cita uma entrada
das Tabelas 3 ou 4 e conferido contra o parquet, nao contra a memoria de quem escreveu.

O QUE ELA PEGA. O estado F regenerou as tabelas e deixou dez citacoes do corpo na versao
anterior: V527 Dra -0,0370 onde a tabela diz -0,0368, 232634196 -0,217 onde diz -0,2181, e por
ai. E a familia dos itens 17-19 do Apendice D com a direcao INVERTIDA: ali a cadeia calculava
certo e o conversor publicava outra coisa; aqui a cadeia calculou certo e o texto ficou velho.
Nenhuma das guardas existentes olhava para isso - `numeros_da_nota`/`confere_nota` conferem as
grandezas REGISTRADAS, e uma citacao de coeficiente no meio de uma frase nao e uma delas.

COMO ELA DECIDE. Frase a frase: quais alvos a frase nomeia (por "TIC <numero>" ou pelo nome
curto da Tabela 3), e entao cada numeral da frase e comparado com as grandezas tabeladas
daqueles alvos (dP/dt, sigma, chi2_red, p_gof, p_curv, p_adv, P, span). Um numeral escrito com d
casas CONFERE se for o arredondamento do valor tabelado com d casas. Se nao for, mas estiver a
menos de 5% dele, e citacao VELHA - o texto ficou para tras. Se estiver longe, e outra grandeza
(o coeficiente do bloco completo, um z, um offset) e a guarda nao opina.

Os 5% sao o vao entre "a mesma grandeza, desatualizada" e "outra grandeza": o maior movimento do
estado F foi 1,11 sigma e ainda assim ficou em 0,4% do valor; o menor par (bloco x desenho) que
o texto poe lado a lado difere por 12%.

EXCECOES. Sao os casos em que um numero LEGITIMO se parece com um desatualizado: o mesmo alvo
medido noutra janela ou noutro catalogo, a poucos por cento do valor tabelado. Sem a razao
escrita, a proxima leitura nao distingue excecao de defeito - por isso cada uma esta aqui:

  230386284, +0,0091: e o coeficiente do BLOCO TESS COMPLETO (5.5), nao o do desenho (+0,0081).
      Os dois convivem na mesma frase, que existe para compara-los, e diferem por 12% - a
      distancia mais curta entre bloco e desenho em toda a amostra, e o motivo de a janela da
      guarda ser de 5% e nao de 15%.
  329246824, -0,0066: e o coeficiente do VarAstro sobre o registro inteiro de 1999-2025 (5.5),
      medido por terceiros e citado ao lado do nosso (-0,0204) para mostrar o desacordo. Nao e a
      nossa medida em versao velha; e outra medida, de outra janela e de outro autor.
  232634196, 4,1: e o chi2_red do ajuste a TUDO (bloco completo mais arquivo, 5.5), nao o
      chi2_red do desenho da Tabela 3 (0,72). O texto rotula "reduced chi2" nos dois casos, e e
      por isso que o rotulo sozinho nao resolve este.
  230386284, +0,0083: e o coeficiente do BLOCO completo com A = 2,0 min na varredura de
      amplitude da 5.3.2 (a frase mostra que o coeficiente quase nao se move enquanto o sigma
      cresce), nao o do desenho (+0,0081); a 2,5% um do outro, dentro da janela da guarda.

Uso:
    python src/citacoes_tabelas.py             # lista as divergencias
    python src/citacoes_tabelas.py --corrigir  # troca cada numeral velho pelo tabelado
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="backslashreplace")
except (AttributeError, ValueError):
    pass

BASE = config.DATA / "orquestra" / "oc_lote"
NOTA = config.ROOT / "reports" / "dpdt_note.md"
JANELA_REL = 0.05          # dentro disto e a mesma grandeza, desatualizada; fora, e outra coisa
JANELA_CHARS = 260         # quantos caracteres em volta do numeral definem "a frase" (ver confere)
ANTES_CHARS = 1200         # ate onde procurar, para tras, a ultima mencao de alvo

# (TIC, numeral como esta escrito): por que este numeral proximo do tabelado NAO e a tabela
EXCECOES = {
    (230386284, "0.0091"): "coeficiente do BLOCO completo (5.5), nao o do desenho",
    (230386284, "+0.0091"): "coeficiente do BLOCO completo (5.5), nao o do desenho",
    (329246824, "0.0066"): "coeficiente do VarAstro no registro inteiro (5.5)",
    (329246824, "−0.0066"): "coeficiente do VarAstro no registro inteiro (5.5)",
    (232634196, "4.1"): "chi2_red do ajuste a TUDO (5.5), nao o do desenho da Tabela 3",
    (230386284, "+0.0083"): "coeficiente do BLOCO completo em A = 2,0 min na varredura de amplitude (5.3.2)",
}

# O corpo cita, em prosa, o COEFICIENTE e a BARRA de um alvo ("+0.0141 +- 0.0016"); as outras
# colunas da Tabela 3 (chi2_red, p, span, P) aparecem so dentro da tabela, e inclui-las aqui
# produzia casamento cruzado entre grandezas de escalas parecidas. Os agregados da Tabela 4 vao
# na segunda parte da guarda, por rotulo.
QUANT = {"dPdt": "dP/dt", "s": "sigma(dP/dt)"}

# Segunda via da guarda: quando o proprio texto ROTULA a grandeza antes do numeral, nao ha
# ambiguidade e a exigencia e mais dura - o numeral tem de ser o arredondamento do tabelado,
# sem janela de tolerancia. (chi2_red de V527 Dra saiu 2,03 no texto e 2,04 na tabela e na
# figura; nenhuma guarda via, porque nao e uma grandeza registrada.)
# "chi2/nu" NAO entra: neste artigo esse simbolo e o do modelo de ruido no bloco TESS
# estendido (5.3.2), nao o chi2_red do ajuste do desenho da Tabela 3
# O vao entre o rotulo e o numeral nao pode ter operador de comparacao (< 0,05 e o LIMIAR do
# teste, nao a medida do alvo) nem a palavra "median" (ai o numero e o resumo de uma
# distribuicao simulada, nao a entrada da tabela).
VAO = r"[^0-9−+<>≤≥]{0,24}$"
ROTULOS = [(r"(?:χ²_red|reduced χ²)(?! median)" + VAO, "chi2r", "chi2_red"),
           (r"p_gof" + VAO, "p_gof", "p_gof"),
           (r"p_curv" + VAO, "p", "p_curv"),
           (r"p_adv" + VAO, "padv", "p_adv")]
NUM = re.compile(r"[−+-]?\d+\.\d+")


def tabela():
    t = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    return t


def corpo(nota):
    """O corpo em PROSA: do inicio ate o Apendice D/E (que guardam o historico e podem citar
    valores velhos de proposito), e sem as tabelas geradas - o que esta dentro dos marcadores e
    das linhas "|" vem do parquet por construcao, e compara-lo consigo mesmo so gera ruido. As
    posicoes sao preservadas (o que sai vira espaco) para que --corrigir edite o lugar certo."""
    fim = min(nota.index("## Appendix C."), nota.index("## Appendix D."))   # C e D guardam o historico
    texto = nota[:fim]
    linhas = []
    dentro = False
    for ln in texto.split(chr(10)):
        if re.match(r"<!-- tabela\d+:inicio", ln):
            dentro = True
        fora = ln.startswith("|") or dentro or ln.startswith("<!--")
        linhas.append(" " * len(ln) if fora else ln)
        if re.match(r"<!-- tabela\d+:fim", ln):
            dentro = False
    return chr(10).join(linhas)


def alvos_da_frase(frase, t):
    achados = set()
    for m in re.finditer(r"TIC (\d{6,10})", frase):
        tic = int(m.group(1))
        if tic in t.index:
            achados.add(tic)
    for tic, nome in t.nome.items():
        if isinstance(nome, str) and len(nome) > 3 and nome in frase:
            achados.add(int(tic))
    return achados


def confere(nota=None, t=None):
    """Devolve a lista de citacoes velhas: (posicao, tic, grandeza, escrito, correto, frase)."""
    nota = nota if nota is not None else NOTA.read_text(encoding="utf-8")
    t = t if t is not None else tabela()
    texto = corpo(nota)
    achados = []
    # JANELA, nao frase: um separador de frases que quebre em "." parte os proprios numerais
    # ("-0.0370" vira "-0" e ".0370") e a guarda nao ve nada - foi o que aconteceu na primeira
    # versao desta funcao, que devolveu zero divergencia num texto com dez
    for mn in NUM.finditer(texto):
        # ASSOCIACAO (rodada dezessete): a janela simetrica de 260 caracteres deixou passar duas
        # citacoes de V527 Dra em que o nome estava 284 e 428 caracteres ANTES do numeral (a
        # abertura da 5.4 e a legenda da Fig. 5). Vale agora a ultima mencao ate ANTES_CHARS
        # atras, mais a janela simetrica curta - a mesma regra do ramo rotulado.
        a0, a1 = max(0, mn.start() - JANELA_CHARS), min(len(texto), mn.end() + JANELA_CHARS)
        frase = texto[a0:a1]
        tics = alvos_da_frase(frase, t)
        if not tics:
            tras = texto[max(0, mn.start() - ANTES_CHARS):mn.start()]
            menc = [(m.start(), tic) for tic in t.index
                    for m in re.finditer(r"TIC " + str(tic) + "|" + re.escape(str(t.loc[tic, "nome"])), tras)
                    if isinstance(t.loc[tic, "nome"], str) and len(str(t.loc[tic, "nome"])) > 3]
            if menc:
                tics = {sorted(menc)[-1][1]}
        if True:
            escrito = mn.group(0)
            if (None, escrito) in EXCECOES:
                continue
            v = float(escrito.replace("−", "-").replace("+", ""))
            casas = len(escrito.split(".")[1])
            antes = texto[max(0, mn.start() - 60):mn.start()]
            rotulado = next(((col, rot) for rx, col, rot in ROTULOS if re.search(rx, antes)), None)
            # "p_curv crosses 0.05", "chi2_red below 1": o numeral e o LIMIAR ou um patamar de
            # comparacao, nao a entrada da tabela daquele alvo
            if rotulado and re.search(r"cross\w*|below|above|under|over|threshold|level|at least", antes[-34:]):
                rotulado = None
            if rotulado:
                col, rot = rotulado
                # com a grandeza rotulada a associacao pode ser mais larga: o PARAGRAFO inteiro
                # (a nota escreve um paragrafo por linha), porque o alvo costuma ser nomeado no
                # inicio dele e o numero aparece frases depois
                # o alvo pode ter sido nomeado paragrafos antes (a 5.4 discute V527 Dra por
                # varias frases antes de citar o chi2_red dele): vale a ULTIMA mencao ate 1500
                # caracteres atras, e so se nenhum outro alvo tiver sido nomeado depois dela
                p0 = max(0, mn.start() - 1500)
                par = texto[p0:mn.start()]
                mencoes = [(m.start(), tic) for tic in t.index
                           for m in re.finditer(r"TIC " + str(tic) + "|" + re.escape(str(t.loc[tic, "nome"])), par)
                           if isinstance(t.loc[tic, "nome"], str) and len(str(t.loc[tic, "nome"])) > 3]
                tics_par = {sorted(mencoes)[-1][1]} if mencoes else set()   # a ULTIMA por posicao
                casa = [x for x in tics_par if abs(abs(v) - abs(round(float(t.loc[x, col]), casas))) < 10 ** (-casas) / 2]
                tic = None if casa or not tics_par else min(
                    tics_par, key=lambda x: min(abs(m.start() + p0 - mn.start()) for m in re.finditer(str(x) + "|" + re.escape(str(t.loc[x, "nome"])), par)))
                if tic is not None and (tic, escrito) not in EXCECOES:
                    alvo = float(t.loc[tic, col])
                    certo = round(alvo, casas)
                    if True:
                        achados.append({"pos": mn.start(), "tic": tic, "grandeza": rot + " (rotulado)",
                                        "escrito": escrito, "correto": f"{certo:.{casas}f}",
                                        "frase": texto[max(0, mn.start() - 70):mn.end() + 70].replace(chr(10), " ")})
                        continue
            # sem rotulo, o corpo escreve coeficiente e barra com 3 ou 4 casas; com 1 ou 2 o
            # numeral e outra coisa (uma fracao, um sigma, um minuto) e a guarda nao opina
            if casas < 3 or not tics:
                continue
            # ordem por PROXIMIDADE, nao por numero: numa frase que nomeia dois alvos, o
            # numeral pertence ao que esta mais perto. Ordenar por TIC fazia a guarda testar o
            # alvo errado primeiro e acusar falso positivo (bullets vizinhos da 5.5).
            rel = mn.start() - max(0, mn.start() - JANELA_CHARS)      # posicao do numeral na janela

            def dist(x):
                # o alvo que vem ANTES do numeral tem prioridade: "TIC A - valor ... - TIC B" e
                # uma frase sobre A, e a mencao de B, logo depois, nao pode roubar o numeral
                mm = [(m.start() - rel) for m in re.finditer(r"TIC " + str(x) + "|" + re.escape(str(t.loc[x, "nome"])), frase)]
                antes = [-d for d in mm if d <= 0]
                if antes:
                    return min(antes)
                return 10 ** 5 + (min(mm) if mm else 10 ** 6)   # so foi achado pela busca para tras
            for tic in sorted(tics, key=dist):
                if (tic, escrito) in EXCECOES:
                    continue
                for col, rotulo in QUANT.items():
                    alvo = float(t.loc[tic, col])
                    if alvo == 0 or abs(v) < 1e-12:
                        continue
                    certo = round(alvo, casas)
                    if abs(abs(v) - abs(certo)) < 10 ** (-casas) / 2:
                        break                      # confere: e o arredondamento do tabelado
                    # "mesma grandeza, desatualizada" = perto em termos RELATIVOS ou a poucos
                    # ulps da precisao escrita (CV Dra -0,0014 -> -0,0013 e 1 ulp, mas 7,7%)
                    perto = max(JANELA_REL * abs(alvo), 1.5 * 10 ** (-casas) if casas >= 4 else 0.0)
                    if abs(abs(v) - abs(alvo)) <= perto:
                        sinal = "−" if (alvo < 0 and escrito.startswith("−")) else ("+" if escrito.startswith("+") else "")
                        achados.append({"pos": mn.start(), "tic": tic, "grandeza": rotulo,
                                        "escrito": escrito, "correto": sinal + f"{abs(certo):.{casas}f}",
                                        "frase": texto[max(0, mn.start() - 70):mn.end() + 70].replace(chr(10), " ")})
                        break
                else:
                    continue
                break
    # uma citacao velha pode casar com mais de um alvo da mesma frase: fica a primeira
    vistos, saida = set(), []
    for a in achados:
        if a["pos"] in vistos:
            continue
        vistos.add(a["pos"])
        saida.append(a)
    return saida


def alvo_mais_proximo(texto, pos, t):
    """A ultima mencao de alvo antes de `pos` (ate ANTES_CHARS atras)."""
    tras = texto[max(0, pos - ANTES_CHARS):pos]
    menc = [(m.start(), tic) for tic in t.index
            for m in re.finditer(r"TIC " + str(tic) + "|" + re.escape(str(t.loc[tic, "nome"])), tras)
            if isinstance(t.loc[tic, "nome"], str) and len(str(t.loc[tic, "nome"])) > 3]
    return sorted(menc)[-1][1] if menc else None


def confere_janela(nota=None, t=None):
    """Terceira via: o z do teste de janela citado no corpo contra ruido_admitido_d9 (95%).

    A nota escreve o z de duas formas: "z = +1.0" quando os modelos admitidos dao o mesmo valor,
    e "z from +3.4 to +6.1" quando dao faixa - as mesmas duas formas da Tabela 5. Aqui as duas
    sao lidas e comparadas com (z_min, z_max) do alvo, a uma casa decimal."""
    nota = nota if nota is not None else NOTA.read_text(encoding="utf-8")
    t = t if t is not None else tabela()
    ra = pd.read_parquet(BASE / "ruido_admitido_d9.parquet")
    ra = ra[ra.nivel == 0.95].set_index("tic")
    # A regiao e a dos paragrafos de classificacao da 5.5, onde o z do TESTE DE JANELA e citado
    # por alvo. Fora dela a nota usa "z" para outras comparacoes (o VarAstro no registro inteiro,
    # o rms da sobreposicao), que nao sao esta grandeza.
    texto = corpo(nota)
    try:
        i0 = texto.index("four of the seven are contradicted")
        i1 = texto.index("**What the block shows")
    except ValueError:
        return []
    regiao = texto[i0:i1]
    achados = []
    for m in re.finditer(r"z (?:=|from) ([−+]?\d+\.\d)(?: to ([−+]?\d+\.\d))?", regiao):
        tic = alvo_mais_proximo(texto, i0 + m.start(), t)
        if tic is None or tic not in ra.index:
            continue
        lo, hi = float(ra.loc[tic, "z_min"]), float(ra.loc[tic, "z_max"])
        esperado = (f"{lo:+.1f}" if abs(hi - lo) < 0.05 else f"{lo:+.1f} to {hi:+.1f}")
        escrito = m.group(0)[2:].replace("= ", "").replace("from ", "")

        def num(x):     # ordenado e sem depender do sinal tipografico: a nota as vezes escreve
            return sorted(float(v.replace("−", "-")) for v in re.findall(r"[−+-]?\d+\.\d", x))   # o hifen ASCII do f-string tambem   # o menor |z| primeiro

        if num(escrito) != num(esperado):
            achados.append({"pos": i0 + m.start(), "tic": int(tic), "grandeza": "z do teste de janela",
                            "escrito": escrito, "correto": esperado,
                            "frase": regiao[max(0, m.start() - 70):m.end() + 70].replace(chr(10), " ")})
    return achados


def main():
    nota = NOTA.read_text(encoding="utf-8")
    t = tabela()
    achados = confere(nota, t)
    janela = confere_janela(nota, t)
    print(f"== z do teste de janela citado no corpo x Tabela 5: {len(janela)} divergencia(s)")
    for j in janela:
        print(f"  TIC {j['tic']} escrito {j['escrito']!r}, tabela {j['correto']!r}   ...{j['frase'][:110]}...")
    print(f"== citacoes do corpo a entradas das Tabelas 3 e 4: {len(achados)} divergencia(s)")
    for a in achados:
        print(f"  TIC {a['tic']} {a['grandeza']:<12s} escrito {a['escrito']:>9s}  tabela {a['correto']:>9s}   ...{a['frase'][:110]}...")
    if not achados:
        print("  nenhuma: todo numeral do corpo que cita a Tabela 3 bate com o parquet")
        return 0
    # o que a guarda encontrou fica gravado, para o Apendice D citar a contagem sem digita-la
    (config.ROOT / "reports" / "citacoes_tabelas.json").write_text(
        json.dumps({"divergencias_iniciais": achados}, indent=1, ensure_ascii=False), encoding="utf-8")
    if "--corrigir" in sys.argv:
        for a in sorted(achados, key=lambda x: -x["pos"]):        # de tras para frente: as posicoes nao andam
            i, j = a["pos"], a["pos"] + len(a["escrito"])
            assert nota[i:j] == a["escrito"], (nota[i:j], a["escrito"])
            nota = nota[:i] + a["correto"] + nota[j:]
        NOTA.write_text(nota, encoding="utf-8")
        print(f"\n  {len(achados)} numeral(is) trocado(s) pelo valor da tabela")
        return 0
    print("\n  (rode com --corrigir para trocar cada um pelo valor da tabela)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
