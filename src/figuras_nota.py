# -*- coding: utf-8 -*-
"""As tres figuras da nota, todas lidas dos artefatos que os scripts da
cadeia gravaram - nada aqui recalcula uma medida.

  fig1  dP/dt dos 26 com barra, ordenados (tabela_26.parquet)
  fig2  leque das 360 fases da injecao do LTTE no V527 Dra (v527_ltte_injetado.parquet)
  fig3  CV Dra: minimos CCD de Brno, parabola do registro inteiro e da nossa
        janela, e as nossas epocas na mesma efemeride (minimos_brno_26.parquet,
        oc__229687624.json, comparar_brno_26.ajustar_parabola)
  fig_completude  completude e supressao por alvo contra |dP/dt| secular
        injetado com barras re-estimadas (completude_secular_26.parquet)
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.ticker  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import comparar_brno_26 as C  # noqa: E402
import config  # noqa: E402

BASE = config.DATA / "orquestra" / "oc_lote"
FIG = config.ROOT / "reports" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
TEXTO_FIG = {}          # nome do arquivo -> todo texto desenhado dentro da figura


def registrar_texto(fig, arquivo):
    """Guarda o texto que a figura desenha (eixos, legendas, tiques, anotacoes).

    O PNG rasteriza esse texto: nenhuma leitura do PDF o alcanca, e foi por ai que rotulos do
    estado A/B ("D10", "1 of 10", "observed -0,0359") sobreviveram quatro rodadas dentro dos
    graficos. Gravado aqui, a lista negra de confere_nota o varre como varre a nota."""
    import matplotlib.text
    fig.canvas.draw()                      # sem isto os rotulos de tique ainda estao vazios
    vistos = []
    for t in fig.findobj(matplotlib.text.Text):
        x = t.get_text().strip()
        if x and x not in vistos:
            vistos.append(x)
    TEXTO_FIG[arquivo] = vistos
    return vistos


def fig1():
    t = pd.read_parquet(BASE / "tabela_26.parquet").sort_values("dPdt").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(4.0, 5.4))   # uma coluna do aa.cls (88 mm)
    y = np.arange(len(t))
    ok = ~t.inadequada
    ax.errorbar(t.dPdt[ok], y[ok], xerr=t.s[ok], fmt="o", color="k", ms=4, capsize=2, lw=0.9, label="quadratic fit adequate (p_gof ≥ 0.05)")
    ax.errorbar(t.dPdt[~ok], y[~ok], xerr=t.s[~ok], fmt="o", mfc="white", mec="k", ms=5, capsize=2, lw=0.9, label="parabola inadequate (p_gof < 0.05)")
    v = t[t.TIC == 424461577]
    ax.plot(v.dPdt, y[v.index], marker="*", ms=14, color="C3", ls="none", label="V527 Dra (published LTTE, no secular term)")
    ax.axvline(0, color="0.6", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{n}" if n else str(tic) for n, tic in zip(t.nome, t.TIC)], fontsize=7)
    ax.set_xlabel("quadratic coefficient as dP/dt (s yr⁻¹)", fontsize=9)
    ax.set_xscale("symlog", linthresh=0.01); ax.tick_params(axis="x", labelsize=8)
    # legenda FORA dos eixos: numa coluna ela cobria os rotulos dos tres ultimos alvos
    # a legenda abaixo do rotulo do eixo x (em -0,09 ela cobria o rotulo; conferido no render)
    ax.legend(fontsize=6.8, loc="upper center", bbox_to_anchor=(0.5, -0.155), frameon=False, ncol=1, handletextpad=0.5)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    registrar_texto(fig, "fig4_dpdt.png")
    fig.savefig(FIG / "fig4_dpdt.png", dpi=180)
    plt.close(fig)


def fig2():
    r = pd.read_parquet(BASE / "v527_ltte_injetado.parquet")
    j = json.loads((BASE / "oc__424461577.json").read_text(encoding="utf-8"))
    obs = j["parabola"]
    fig, axs = plt.subplots(3, 1, figsize=(7.5, 7.5), sharex=True)
    axs[0].plot(r.fase, r.dPdt, color="k", lw=1.2)
    axs[0].fill_between(r.fase, r.dPdt - r.sdPdt, r.dPdt + r.sdPdt, color="0.8")
    axs[0].axhline(obs["dPdt_s_por_ano"], color="C3", ls="--", label=f"observed {obs['dPdt_s_por_ano']:+.4f} ± {obs['sdPdt_s_por_ano']:.4f}")
    axs[0].axhline(0, color="0.6", lw=0.8)
    axs[0].set_ylabel("dP/dt (s yr⁻¹)")
    axs[0].legend(fontsize=8)
    axs[1].semilogy(r.fase, r.p_curv, color="k", lw=1.2, label="p_curv")
    axs[1].semilogy(r.fase, r.p_adv, color="C0", lw=1.2, label="p_adv")
    axs[1].axhline(0.05, color="0.6", ls=":", label="0.05")
    axs[1].set_ylim(1e-56, 1.0)          # p > 1 e impossivel: o topo do eixo e 10^0
    axs[1].set_ylabel("p")
    axs[1].legend(fontsize=8, ncol=3)
    axs[2].plot(r.fase, r.chi2_red, color="k", lw=1.2)
    axs[2].axhline(obs["chi2"] / obs["dof"], color="C3", ls="--", label=f"observed {obs['chi2'] / obs['dof']:.2f}")
    axs[2].set_ylabel("χ²_red")
    axs[2].set_xlabel("orbital phase of the injected LTTE (P₃ = 2.733 yr, A = 6.6 min)")
    axs[2].legend(fontsize=8)
    fig.tight_layout()
    registrar_texto(fig, "fig2_v527_injecao.png")
    fig.savefig(FIG / "fig2_v527_injecao.png", dpi=180)
    plt.close(fig)


def fig3():
    mins = pd.read_parquet(BASE / "minimos_brno_26.parquet")
    m = mins[(mins.tic == 229687624) & (mins.tipo == "P") & mins.metodo.str.contains("ccd", na=False)].sort_values("jd")
    P_va, m0 = float(m.P_va.iloc[0]), float(m.m0_va.iloc[0])
    ano = C._ano
    E = np.round((m.jd.values - m0) / P_va)
    coef_t, *_ = C.ajustar_parabola(E, m.oc_d.values)
    a0, a1 = C._janela_nossa(229687624)
    mj = m[(ano(m.jd) >= a0 - 0.5) & (ano(m.jd) <= a1 + 0.5)]
    Ej = np.round((mj.jd.values - m0) / P_va)
    coef_j, *_ = C.ajustar_parabola(Ej, mj.oc_d.values)
    j = json.loads((BASE / "oc__229687624.json").read_text(encoding="utf-8"))
    p = pd.DataFrame(j["pontos"])
    jd_n = p.t0.values + 2457000.0 - C.TT_UTC_D
    En = np.round((jd_n - m0) / P_va)
    oc_n = (jd_n - (m0 + P_va * En)) * 1440
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.plot(ano(m.jd), m.oc_d * 1440, "o", ms=3.5, color="0.45", label="VarAstro CCD primary minima (83)")
    Eg = np.linspace(E.min(), E.max(), 400)
    ax.plot(ano(m0 + P_va * Eg), np.polyval(coef_t[::-1], Eg) * 1440, "-", color="k", lw=1.3, label="quadratic, whole record 1988–2022 (+0.0088 ± 0.0005 s yr⁻¹)")
    jd_a0, jd_a1 = 2451545.0 + (a0 - 2000.0) * 365.25, 2451545.0 + (a1 - 2000.0) * 365.25
    Ejg = np.linspace((jd_a0 - m0) / P_va, (jd_a1 - m0) / P_va, 200)   # a janela inteira, 2004-2024
    ax.plot(ano(m0 + P_va * Ejg), np.polyval(coef_j[::-1], Ejg) * 1440, "--", color="C0", lw=1.6, label="quadratic, our window 2004–2024 (−0.0003 ± 0.0017 s yr⁻¹)")
    ax.errorbar(ano(jd_n), oc_n, yerr=p.sig_min, fmt="s", color="C3", ms=5, capsize=2, label="this work: SuperWASP seasons and TESS sectors on the same ephemeris")
    ax.axvspan(a0, a1, color="C0", alpha=0.07)
    ax.set_xlabel("year")
    ax.set_ylabel("O − C (min), VarAstro ephemeris")
    ax.legend(fontsize=7.5, loc="upper center")
    fig.tight_layout()
    registrar_texto(fig, "fig3_cvdra.png")
    fig.savefig(FIG / "fig3_cvdra.png", dpi=180)
    plt.close(fig)


def fig_barras():
    """A figura que a tese pede: epsilon e |dP/dt|/sigma contra a barra
    mediana do alvo, com D11 marcado e, se existir, quem sobrevive a
    reinflacao das barras do TESS em cada cenario."""
    t = pd.read_parquet(BASE / "tabela_26.parquet")
    e = pd.read_parquet(BASE / "epsilon_26.parquet").set_index("tic")
    t["eps"] = t.TIC.map(e.epsilon); t["barra"] = t.TIC.map(e.barra_mediana_min); t["z"] = t.dPdt.abs() / t.s
    rf = BASE / "reinflar_tess_26.parquet"
    R = pd.read_parquet(rf) if rf.exists() else None
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.4))
    for ax, col, lab in ((axs[0], "eps", "detection efficiency ε (Sect. 5.4)"), (axs[1], "z", "|dP/dt| / σ")):
        d11 = t.curvatura
        ax.plot(t.barra[~d11], t[col][~d11], "o", mfc="white", mec="0.4", ms=6, label="no surviving curvature")
        ax.plot(t.barra[d11], t[col][d11], "o", color="k", ms=6, label=f"D{int(d11.sum())}: curvature survives adversarial")
        v = t[t.TIC == 424461577]
        ax.plot(v.barra, v[col], "*", color="C3", ms=14, ls="none", label="V527 Dra (published LTTE)")
        # (o marcador "lost with SuperWASP bars x1.65" saiu em 2026-09-13: nao ha mais correcao SuperWASP)
        n11 = int(d11.sum())
        if R is not None:
            for d, mk, cor in ((1.0, "s", "C0"), (2.6, "^", "C1")):
                g = R[(R.cenario_min_por_ano == d) & R.em_D11 & R.detecta]
                sel = t[t.TIC.isin(g.tic)]
                ax.plot(sel.barra, sel[col], mk, mfc="none", mec=cor, ms=13, mew=1.4, ls="none",
                        label=f"still detected with TESS bars +{d} min yr⁻¹, sensitivity ceiling ({len(sel) if len(sel) else 'none'} of {n11})")
        ax.set_xscale("log"); ax.set_xlabel("median epoch bar of the target (min)")
        ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5, 10]); ax.set_xticklabels(["0.1", "0.2", "0.5", "1", "2", "5", "10"])
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_ylabel(lab)
        if col == "z":
            ax.set_yscale("log"); ax.set_ylim(0.1, 60)
            ax.axhline(3, color="0.6", ls=":", lw=0.8, label="|dP/dt|/σ = 3")
    handles, labels = axs[0].get_legend_handles_labels()
    h2, l2 = axs[1].get_legend_handles_labels()
    for h, l in zip(h2, l2):
        if l not in labels:
            handles.append(h); labels.append(l)
    fig.legend(handles, labels, fontsize=7, loc="lower center", ncol=2, frameon=False)   # fora dos paineis: nao cobre ponto nenhum
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    registrar_texto(fig, "fig1_barras.png")
    fig.savefig(FIG / "fig1_barras.png", dpi=180)
    plt.close(fig)


def fig_completude():
    """Fig. 4 da nota: por alvo, completude (painel de cima) e supressao sigma_re/sigma_verdadeiro (de baixo)
    contra o |dP/dt| secular injetado; o D11 em preto, os outros em cinza, a mediana dos 26 em vermelho.
    Uma coluna, paineis empilhados (rodada dez-b)."""
    R = pd.read_parquet(BASE / "completude_secular_26.parquet")
    t = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    fig, axs = plt.subplots(2, 1, figsize=(4.0, 5.8))   # uma coluna, paineis empilhados
    grade = sorted(R.dpdt.unique()); x = [g if g > 0 else 0.002 for g in grade]      # o 0 (falso alarme) entra como 0,002 no eixo log
    for tic, g in R.groupby("tic"):
        g = g.sort_values("dpdt")
        d10 = bool(t.loc[tic, "curvatura"])
        kw = dict(color="k", lw=1.2, alpha=0.9, zorder=3) if d10 else dict(color="0.65", lw=0.8, alpha=0.8, zorder=2)
        axs[0].plot(x, g.completude_re, "-", **kw); axs[1].plot(x, g.supressao_med, "-", **kw)
    med = R.groupby("dpdt").agg(c=("completude_re", "median"), s=("supressao_med", "median")).sort_index()
    axs[0].plot(x, med.c, "-", color="C3", lw=2.2, zorder=4, label="median of the 26"); axs[1].plot(x, med.s, "-", color="C3", lw=2.2, zorder=4)
    axs[0].plot([], [], "-", color="k", lw=1.2, label=f"D{int(t.curvatura.sum())}"); axs[0].plot([], [], "-", color="0.65", lw=0.8, label="other targets")
    for ax in axs:
        ax.set_xscale("log"); ax.set_xticks([0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]); ax.set_xticklabels(["0", "0.005", "0.01", "0.02", "0.05", "0.1", "0.2"])
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_xlabel("injected |dP/dt| (s yr⁻¹)", fontsize=9)
        ax.tick_params(labelsize=8)
    axs[0].set_ylabel("completeness", fontsize=9); axs[0].set_ylim(-0.02, 1.02); axs[0].axhline(0.5, color="0.6", ls=":", lw=0.8)
    axs[1].set_yscale("log"); axs[1].set_ylabel("σ re-estimated / σ true", fontsize=9); axs[1].axhline(1, color="0.6", ls=":", lw=0.8)
    axs[1].set_yticks([0.5, 1, 2, 5, 10, 20]); axs[1].set_yticklabels(["0.5", "1", "2", "5", "10", "20"]); axs[1].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    axs[0].legend(fontsize=8, loc="lower right", frameon=False)
    fig.tight_layout()
    registrar_texto(fig, "fig5_completude.png")
    fig.savefig(FIG / "fig5_completude.png", dpi=180)
    plt.close(fig)


def gravar_texto_das_figuras():
    destino = config.ROOT / "reports" / "figuras_texto.json"
    destino.write_text(json.dumps(TEXTO_FIG, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"  texto das figuras -> {destino} ({sum(len(v) for v in TEXTO_FIG.values())} cadeias)")


def espelhar_no_paper():
    """As figuras do PDF vem de paper/figures. Ate a rodada dez elas eram copias de mao, e o
    PDF ficou com as figuras do estado A/B enquanto o texto ja estava no E. Agora toda geracao
    copia, e o conversor confere byte a byte (md_para_tex)."""
    import shutil
    destino = config.ROOT / "paper" / "figures"
    destino.mkdir(parents=True, exist_ok=True)
    for p in sorted(FIG.glob("*.png")):
        shutil.copy2(p, destino / p.name)
        assert (destino / p.name).read_bytes() == p.read_bytes(), p.name
    print(f"  espelhadas em {destino}: {len(list(FIG.glob('*.png')))} figuras")


if __name__ == "__main__":
    fig_barras(); fig1(); fig2(); fig3(); fig_completude()
    print("figuras em", FIG)
    espelhar_no_paper()
    gravar_texto_das_figuras()
