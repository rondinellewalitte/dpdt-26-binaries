# -*- coding: utf-8 -*-
"""Qual barra esta pequena: TESS ou SuperWASP? Quarta rodada, item 1.

O chi2_red mediano de 1,41 (contra 0,73 esperado) diz que a variancia TOTAL
esta ~1,4x subestimada, nao QUAL barra. A nota atribuiu tudo ao TESS com
base na deriva do secundario em 5 de 24 alvos, enquanto a Secao 3.4
declara que a barra SuperWASP carrega um termo de tamanho desconhecido.
Inflar o SuperWASP tambem levaria a mediana a 0,73 e NAO mataria as
deteccoes, que sao internas ao TESS. A manchete depende de qual barra foi
inflada, e a escolha precisa ser justificada pelos dados.

Teste: residuos do ajuste quadratico de cada alvo (gravados no JSON, na
ordem de `pontos`) normalizados pela barra da propria epoca, separados em
epocas TESS e SuperWASP, sobre os 26. Com 3 parametros sobre N = 4-7
pontos os residuos sao encolhidos pela alavanca: E[res_i^2/sigma_i^2] =
1 - h_i, com h_i a diagonal da matriz-chapeu do ajuste ponderado. A
grandeza reportada e o rms de z_i/sqrt(1 - h_i), que vale 1 com barras
corretas, em cada conjunto, com n e dois intervalos: o de chi2 (n r^2 /
chi2_{n; 0,975/0,025}) e o de bootstrap por ALVO (as epocas de um alvo
nao sao independentes entre si). Tambem o rms cru e a decomposicao do
Sigma chi2 = 128,6 (67 gl) entre os dois conjuntos.

Limite do teste, dito antes: a barra SuperWASP e construida DA dispersao
entre temporadas (hypot(max(formal, dispersao), 0,78)), entao os residuos
SuperWASP dentro do grupo sao <= 1 quase por construcao. O teste tem
poder para localizar excesso no TESS; para o SuperWASP ele so ve o que
sobra alem da dispersao que ja entrou na barra. Um deslocamento de GRUPO
do SuperWASP (vies) nao aparece como residuo - a parabola o absorve - e
so um registro externo o veria. Isso fica no relato.

EXPECTATIVA (escrita e commitada antes de rodar): rms corrigido do TESS
entre 1,3 e 1,8 (n ~ 90 epocas); rms corrigido do SuperWASP entre 0,8 e
1,1 (n = 61), em parte por construcao; a decomposicao do Sigma chi2 poe
a maior parte do excesso sobre o esperado no TESS. Se o excesso sair
repartido (SuperWASP tambem > 1,2 com intervalo fora de 1), 1,0 min/ano
e teto, nao estimativa, e a Secao 5 e reescrita - decisao do revisor.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"


def alavancas(E, sig_d):
    """diagonal da matriz-chapeu do ajuste quadratico ponderado (mesma
    construcao de cadeia_oc.ajustar: vander(E, 3) * (1/sigma))."""
    w = 1.0 / sig_d
    A = np.vander(E, 3) * w[:, None]
    H = A @ np.linalg.inv(A.T @ A) @ A.T
    return np.diag(H)


def intervalo_chi2(z2, n):
    r2 = np.mean(z2)
    return np.sqrt(n * r2 / stats.chi2.ppf(0.975, n)), np.sqrt(n * r2 / stats.chi2.ppf(0.025, n))


if __name__ == "__main__":
    t26 = pd.read_parquet(BASE / "tabela_26.parquet")
    linhas = []
    for tic in sorted(t26.TIC):
        j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
        p = pd.DataFrame(j["pontos"])
        res = np.array(j["parabola"]["residuos_min"], float)
        assert len(res) == len(p), tic
        h = alavancas(p.E.values.astype(float), p.sig_min.values / 1440.0)
        for k in range(len(p)):
            linhas.append({"tic": int(tic), "fonte": "TESS" if p.fonte[k].startswith("TESS") else "SuperWASP",
                           "res_min": res[k], "sig_min": float(p.sig_min[k]), "h": float(h[k]),
                           "z": res[k] / float(p.sig_min[k]), "z_corr": res[k] / float(p.sig_min[k]) / np.sqrt(1 - h[k])})
    d = pd.DataFrame(linhas)
    d.to_parquet(BASE / "qual_barra_26.parquet", index=False)
    chi2_total = float((d.z ** 2).sum())
    print(f"Sigma chi2 = {chi2_total:.1f} sobre {len(d)} epocas, {len(d) - 3 * 26} gl (chi2_red global {chi2_total / (len(d) - 3 * 26):.2f}, "
          f"p = {stats.chi2.sf(chi2_total, len(d) - 3 * 26):.1e})")
    rng = np.random.default_rng(4)
    tics = d.tic.unique()
    print("\n== rms dos residuos normalizados por conjunto (1 = barras corretas)")
    print(f"{'conjunto':10s} {'n':>4s} {'rms cru':>8s} {'rms corr.':>10s} {'IC chi2':>14s} {'IC bootstrap/alvo':>18s} {'Sigma z2':>9s} {'esperado':>9s} {'excesso':>8s}")
    for fonte, g in d.groupby("fonte"):
        n = len(g)
        rms_cru = float(np.sqrt(np.mean(g.z ** 2)))
        rms_c = float(np.sqrt(np.mean(g.z_corr ** 2)))
        lo, hi = intervalo_chi2(g.z_corr ** 2, n)
        boot = []
        for _ in range(5000):
            amostra = rng.choice(tics, len(tics), replace=True)
            gg = pd.concat([g[g.tic == t] for t in amostra])
            boot.append(np.sqrt(np.mean(gg.z_corr ** 2)))
        blo, bhi = np.percentile(boot, [2.5, 97.5])
        esperado = float((1 - g.h).sum())
        # estimador de razao de somas: Sigma z2 / Sigma (1 - h) - e a decomposicao
        # do chi2 total, e nao sofre dos denominadores (1 - h) pequenos das
        # epocas de alavanca alta (TESS: 4 pontos, 3 parametros). CI por
        # bootstrap sobre alvos e por chi2 com gl efetivo = esperado.
        razao = float((g.z ** 2).sum()) / esperado
        rlo, rhi = esperado * razao / stats.chi2.ppf(0.975, esperado), esperado * razao / stats.chi2.ppf(0.025, esperado)
        bootr = []
        for _ in range(5000):
            amostra = rng.choice(tics, len(tics), replace=True)
            gg = pd.concat([g[g.tic == t] for t in amostra])
            bootr.append((gg.z ** 2).sum() / (1 - gg.h).sum())
        rblo, rbhi = np.percentile(bootr, [2.5, 97.5])
        print(f"{fonte:10s} {n:4d} {rms_cru:8.2f} {rms_c:10.2f} {lo:6.2f}–{hi:<6.2f} {blo:8.2f}–{bhi:<8.2f} {float((g.z ** 2).sum()):9.1f} {esperado:9.1f} {float((g.z ** 2).sum()) - esperado:+8.1f}"
              f"   | fator de variancia Sigma z2/Sigma(1-h) = {razao:.2f} (chi2 {rlo:.2f}–{rhi:.2f}; bootstrap/alvo {rblo:.2f}–{rbhi:.2f}); fator de barra sqrt = {np.sqrt(razao):.2f}")
    print("\n  (Sigma z2 = soma dos residuos normalizados ao quadrado; esperado = Sigma (1 - h) com barras corretas; excesso = diferenca)")
    print(f"  soma dos esperados = {float((1 - d.h).sum()):.1f} = gl; soma dos excessos = {chi2_total - float((1 - d.h).sum()):+.1f}")
    print(f"  -> {BASE / 'qual_barra_26.parquet'}")
