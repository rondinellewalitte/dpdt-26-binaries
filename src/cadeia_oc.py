# -*- coding: utf-8 -*-
"""A cadeia de O-C inteira, para UM alvo, como uma funcao - e o controle que
decide se ela vale.

Tudo aqui ja existia e esta validado: `epocas_142874476.carregar/medir` (epoca
por setor), `superwasp.carregar` (BJD_TDB), `superwasp.epoca_por_subconjunto`
(barra do ruido vermelho com gate de cobertura por bloco), o vies da cadeia
(-1,26 +- 0,75 min no estado E; -2,14 +- 0,78 ate o estado D; 4 Jupiteres quentes). O que NAO existia era uma funcao que
encadeasse isso e devolvesse os componentes - o 11,2 sigma vivia num
comentario. Antes de rodar em 88 alvos, ela roda no TYC 7024-1046-1 e cada
componente e comparado com `expectativa_oc_142874476.json`, escrito antes.

Duas regras do desenho que o corpo implementa:
  - dP/dt sem grau de liberdade nao e testado: com k aglomerados e 3
    parametros, dof = k - 3; abaixo de 1 o alvo sai como "medido, nao testado".
  - O-C grande nao e achado sem contagem de ciclos: a escada resolve E antes
    de qualquer barra, e o E do SuperWASP e reportado.
"""
import json
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))

import config
import superwasp as sw

# ESTADO E (2026-09-17): vies_hj_indep - media ponderada do O-C das mesmas temporadas SuperWASP dos 4 Jupiteres
# quentes contra efemerides publicadas com dados de solo (Ivshina & Winn 2022 + Bouma et al. 2020 no WASP-4 b):
# -1,26 +- 0,75 min. Ate o estado D era -2,14 +- 0,78 (run_hj_calib: propagacao TESS-only, que levava o
# decaimento orbital do WASP-4 b para dentro do vies).
VIES_CADEIA_MIN = -1.26
VIES_CADEIA_SIG = 0.75

# ESTADO F (2026-09-17, rodada treze, bloco C): alem do vies COMUM acima, cada alvo com teste de
# forma (vies_forma_swasp) ganha a correcao de forma DO PROPRIO ALVO nas epocas SuperWASP, e o
# piso da barra arquival deixa de ser o 0,75 min sozinho e passa a ser hypot(0,75; sigma da media
# do vies de forma). A tabela por alvo e escrita por `vies_forma_por_alvo.py` pela regra declarada
# na rodada treze; alvo sem medida - os outros 62 dos 88 e o CONTROLE 142874476 - nao muda em nada,
# e e por isso que o controle tem de reproduzir depois da troca.
VIES_FORMA_TABELA = config.DATA / "orquestra" / "vies_forma_por_alvo.json"


@lru_cache(maxsize=1)
def _tabela_vies_forma():
    if not VIES_FORMA_TABELA.exists():
        raise FileNotFoundError(f"estado F: falta {VIES_FORMA_TABELA} (rode src/vies_forma_por_alvo.py)")
    return json.loads(VIES_FORMA_TABELA.read_text(encoding="utf-8"))["alvos"]


def vies_forma(tic):
    """(correcao em min, piso da barra arquival em min) do alvo. Sem medida: (0, VIES_CADEIA_SIG)."""
    d = _tabela_vies_forma().get(str(int(tic)))
    if d is None:
        return 0.0, VIES_CADEIA_SIG
    return float(d["vies_min"]), float(np.hypot(VIES_CADEIA_SIG, float(d["sig_vies_min"])))


def ajustar(E, t, sig, grau):
    """Minimos quadrados ponderados de t(E). Devolve (coef, residuos, chi2, cov)."""
    E, t, sig = map(np.asarray, (E, t, sig))
    w = 1.0 / sig
    A = np.vander(E, grau + 1) * w[:, None]
    coef, *_ = np.linalg.lstsq(A, t * w, rcond=None)
    modelo = np.polyval(coef, E)
    res = t - modelo
    chi2 = float(np.sum((res / sig) ** 2))
    cov = np.linalg.inv(A.T @ A)
    return coef, res, chi2, cov


def epocas_tess(tic):
    """Reusa carregar()+medir() do modulo validado, sem tocar no parquet dele."""
    import epocas_142874476 as ep
    if tic != ep.TIC:
        raise RuntimeError(f"FALHA [cadeia]: extracao TESS so existe para {ep.TIC}; "
                           "generalizar e a Tarefa 2, e ela vem depois da reproducao")
    fontes, vistos = [], set()
    for s in ep.carregar():
        k = (s["setor"], s["reducao"])
        if k not in vistos:
            vistos.add(k)
            fontes.append(s)
    anc = [s for s in fontes if s["setor"] == 105][0]
    t0_anc, _, _ = ep.medir(anc["t"], anc["f"])
    linhas = []
    for s in fontes:
        t0, sig, prof = ep.medir(s["t"], s["f"], semente=t0_anc)
        linhas.append({"setor": s["setor"], "reducao": s["reducao"], "t0_btjd": t0,
                       "sigma_min": sig, "n": s["n"], "cadencia_min": s["cadencia_min"]})
    return pd.DataFrame(linhas).sort_values(["setor", "reducao"]), ep.P_REF


def aglomerar(ep_df):
    """Um ponto por setor. Com varias reducoes: media, e a barra e o ESPALHAMENTO."""
    out = []
    for s, g in ep_df.groupby("setor"):
        if len(g) == 1:
            # as duas componentes da barra (formal; dispersao entre metades) seguem
            # junto quando existem: o nulo com barras re-estimadas precisa delas
            out.append({"setor": s, "t0": float(g.t0_btjd.iloc[0]), "sig_min": float(g.sigma_min.iloc[0]),
                        "n_reducoes": 1,
                        "sig_formal_min": float(g.sigma_formal_min.iloc[0]) if "sigma_formal_min" in g else float("nan"),
                        "sig_metades_min": float(g.sigma_metades_min.iloc[0]) if "sigma_metades_min" in g else float("nan")})
        else:
            t0s = g.t0_btjd.values
            # A barra com varias reducoes E o espalhamento entre elas - no
            # controle, 0,51 min entre QLP/TGLC/eleanor, MENOR que as formais
            # (1,15-1,60), e e a escolha validada contra a simulacao (0,47).
            # Por um piso na formal aqui moveu o controle inteiro (ESCALAR em
            # 18 campos) e foi revertido. Espalhamento ZERO nao e barra: e o
            # mesmo arquivo lido duas vezes, e isso levanta em vez de virar
            # peso infinito na escada.
            esp = float(np.std(t0s, ddof=1) * 1440.0)
            if not esp > 0:
                raise RuntimeError(f"setor {s}: {len(g)} reducoes com epocas IDENTICAS "
                                   f"(espalhamento 0) - o mesmo arquivo entrou mais de uma vez")
            out.append({"setor": s, "t0": float(np.mean(t0s)), "sig_min": esp, "n_reducoes": len(g)})
    return pd.DataFrame(out).sort_values("setor")


def medir_alvo(tic, sourceid, ra, dec, t14_h=None):
    """O alvo de controle, pela extracao especifica dele (FFI em cache)."""
    ep_df, P0 = epocas_tess(tic)
    return medir_alvo_de(tic, ep_df, P0, sourceid, ra, dec, t14_h)


def medir_alvo_de(tic, ep_df, P0, sourceid, ra, dec, t14_h=None):
    """A cadeia a partir de uma tabela de epocas TESS ja medida.

    Separar a extracao (que difere entre FFI em cache e SPOC 2-min dos
    manifestos) do resto da cadeia e o que permite distinguir "divergiu porque
    o codigo mudou" de "divergiu porque o dado de entrada nao e o mesmo".
    """
    ag = aglomerar(ep_df)
    # ancora no setor MAIS RECENTE - o controle tinha o s105, os 88 nao
    # necessariamente; um `== 105` fixo derrubou os 3 primeiros do piloto
    t0_ref = float(ag.loc[ag.setor.idxmax(), "t0"])
    ag["E"] = np.round((ag.t0 - t0_ref) / P0).astype(int)
    sig_d = ag.sig_min.values / 1440.0
    cl, rl, chi2l, _ = ajustar(ag.E.values, ag.t0.values, sig_d, 1)
    cp, rp, chi2p, _ = ajustar(ag.E.values, ag.t0.values, sig_d, 2)
    k = len(ag)
    res = {
        "tic": tic, "P_ref_d": P0,
        "tess": {"n_aglomerados": k,
                 "aglomerados": ag.assign(res_linear_min=rl * 1440, res_parabola_min=rp * 1440)
                 .to_dict("records"),
                 "reducoes_s4": ep_df[ep_df.setor == 4][["reducao", "t0_btjd", "sigma_min"]].to_dict("records"),
                 "linear": {"chi2": chi2l, "dof": k - 2, "P_d": float(cl[0])},
                 "parabola": {"chi2": chi2p, "dof": k - 3, "Q_d_por_ciclo2": float(cp[0])},
                 "delta_chi2": chi2l - chi2p,
                 "p_parabola_vs_linear": float(stats.chi2.sf(chi2l - chi2p, 1)) if k >= 4 else float("nan"),
                 "testavel": k - 3 >= 1},
    }

    # ---- SuperWASP: a epoca de 2006-2008, com a barra que o ruido vermelho exige
    t, f, sig = sw.carregar(sourceid, ra, dec)
    P_lin = float(cl[0])
    t_med = float(np.median(t))
    # semente: a efemeride TESS propagada; t em BJD, TESS em BTJD (= BJD - 2457000)
    E_sw = int(np.round((t_med - (t0_ref + 2457000.0)) / P_lin))
    t0_prev = t0_ref + 2457000.0 + E_sw * P_lin
    if t14_h is None:
        t14_h = sw.largura_do_minimo(t, f, P_lin)
    t0g, sigf, sig_sub, subs, status = sw.epoca_por_subconjunto(
        t, f, P_lin, t0_prev, t14_h=t14_h, n_sub=2, modo="temporada", exigir_cobertura=True)
    barra_interna = max(float(sigf), float(sig_sub)) if np.isfinite(sig_sub) else float("nan")
    barra_total = float(np.hypot(barra_interna, VIES_CADEIA_SIG)) if np.isfinite(barra_interna) else float("nan")
    E_real = int(np.round((t0g - (t0_ref + 2457000.0)) / P_lin))
    # O-C contra as duas efemerides dos pontos TESS, com o vies da cadeia tirado
    oc_lin = (t0g - (np.polyval(cl, E_real) + 2457000.0)) * 1440.0 - VIES_CADEIA_MIN
    oc_par = (t0g - (np.polyval(cp, E_real) + 2457000.0)) * 1440.0 - VIES_CADEIA_MIN
    prev_par_2007 = (np.polyval(cp, E_real) - np.polyval(cl, E_real)) * 1440.0
    res["superwasp"] = {
        "sourceid": sourceid, "npts": int(len(t)), "mad_mmag": sig * 1000,
        "t14_h": float(t14_h), "n_temporadas": len(subs), "status_cobertura_por_bloco": status,
        "epoca_bjd_tdb": float(t0g), "E": E_real, "barra_formal_min": float(sigf),
        "barra_entre_temporadas_min": float(sig_sub), "barra_interna_min": barra_interna,
        "vies_cadeia_min": VIES_CADEIA_MIN, "barra_total_min": barra_total,
        "temporadas": [{k2: (round(v, 4) if isinstance(v, float) else v) for k2, v in s.items()} for s in subs],
    }
    res["teste"] = {
        "residuo_no_linear_min": float(oc_lin), "sigmas_linear": float(oc_lin / barra_total),
        "residuo_na_parabola_min": float(oc_par), "sigmas_parabola": float(abs(oc_par) / barra_total),
        "parabola_prevista_min": float(prev_par_2007),
        "graus_de_liberdade_com_swasp": k + 1 - 3,
        "classe": "medido, nao testado" if k + 1 - 3 < 1 else "testado",
    }
    return res, ag, ep_df


def adversarial(ag, P0, t0_ref):
    """Desloca cada epoca 1 sigma na direcao que mais ENFRAQUECE a curvatura."""
    E = ag.E.values.astype(float)
    t0 = ag.t0.values.copy()
    sig = ag.sig_min.values / 1440.0
    _, _, chi2l, _ = ajustar(E, t0, sig, 1)
    _, _, chi2p, _ = ajustar(E, t0, sig, 2)
    base = chi2l - chi2p
    pior = t0.copy()
    for i in range(len(t0)):
        melhor_d, melhor = None, base
        for sgn in (-1, 1):
            tt = pior.copy()
            tt[i] += sgn * sig[i]
            _, _, a, _ = ajustar(E, tt, sig, 1)
            _, _, b, _ = ajustar(E, tt, sig, 2)
            if a - b < melhor:
                melhor, melhor_d = a - b, sgn
        if melhor_d is not None:
            pior[i] += melhor_d * sig[i]
    _, _, a, _ = ajustar(E, pior, sig, 1)
    _, _, b, _ = ajustar(E, pior, sig, 2)
    return {"delta_chi2_base": float(base), "p_base": float(stats.chi2.sf(base, 1)),
            "delta_chi2_adversarial": float(a - b), "p_adversarial": float(stats.chi2.sf(a - b, 1)),
            "deslocamentos_min": ((pior - t0) * 1440.0).round(3).tolist()}


if __name__ == "__main__":
    exp = json.load(open(config.DATA / "orquestra" / "expectativa_oc_142874476.json", encoding="utf-8"))
    # posicao: a que o identificador 1SWASP codifica (a 0,15" do alvo)
    import superwasp_cobertura as sc
    ra, dec = sc._coord_do_id(exp["superwasp"]["sourceid"])
    res, ag, ep_df = medir_alvo(142874476, exp["superwasp"]["sourceid"], ra, dec)
    res["adversarial"] = adversarial(ag, res["P_ref_d"], float(ag[ag.setor == 105].t0.iloc[0]))
    out = config.DATA / "orquestra" / "cadeia_oc__142874476.json"
    out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
    print(json.dumps(res, indent=2, default=float))
    print("\ngravado:", out)
