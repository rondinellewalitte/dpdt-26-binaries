# -*- coding: utf-8 -*-
"""A cadeia de O-C em lote: SPOC 2-min dos manifestos + SuperWASP, por alvo.

O que muda em relacao a `cadeia_oc.medir_alvo` (o alvo de controle) e SO a
extracao das epocas TESS: la, FFI de tres reducoes num cache proprio; aqui,
os arquivos SPOC de 2 min que os 22 manifestos listam, baixados por URL. A
partir da tabela de epocas, a cadeia e a mesma funcao (`medir_alvo_de`).

DOIS CONTROLES, que distinguem defeito de generalizacao:
  A. `medir_alvo_de` alimentada com as epocas ESPECIFICAS do TIC 142874476
     tem que reproduzir a expectativa v2 - testa o trecho compartilhado.
  B. o leitor SPOC 2-min + `medir` aplicados ao MESMO arquivo s105 do alvo
     de controle tem que devolver a epoca do s105 (4219,602785) - testa a
     extracao generica onde ha resposta conhecida. O alvo de controle so
     aparece em UM dos 22 manifestos (s105), entao a cadeia generica
     inteira nele nao teria aglomerados para comparar - e por isso o
     controle e do elo, nao da cadeia.

Por alvo, o produto e distribuicao, na ordem de solidez do desenho:
  1. dP/dt com barra e dof (a distribuicao na populacao e resultado por si);
  2. curvatura que sobrevive ao adversarial (o subconjunto forte);
  3. "medido, nao testado" declarado como tal (dof < 1).
A convencao da parabola e a de `cadeia_oc.ajustar` - a unica reconstruivel -
e cada linha de saida a carrega por nome.

DEFEITO NA COSTURA ENTRE DOIS CHAMADORES (achado em 2026-09-13, revisao
externa, item 1): `superwasp.epoca_por_subconjunto` devolve `sig_sub` =
std(desvios das temporadas, ddof=1) / sqrt(n) - o ERRO DA MEDIA das
temporadas, correto para o desenho de `cadeia_oc.medir_alvo_de` em que o
bloco SuperWASP entra como UMA epoca combinada (o caminho do controle). Este
lote poe CADA temporada como epoca propria e dava a cada uma esse mesmo
`sig_sub`: a dispersao de uma temporada e std(desvios), nao std/sqrt(n).
Com n = 2 (17 alvos) ou 3 (9 alvos), barra sqrt(2)-sqrt(3) pequena demais
em 49 das 61 epocas SuperWASP dos 26. So de leitura dos JSONs: a barra
certa supera a usada por mediana 1,41 (quartis 1,28-1,63, max 1,73); em
variancia, ponderado pelos pesos do ajuste, 1,74 - contra o fator 2,73
(1,79-4,64) que `qual_barra_26` mediu e a nota chamou de "barras SuperWASP
1,65x pequenas demais". O comentario acima de `disp_sw` dizia "a dispersao";
o codigo passava a dispersao/sqrt(n).

EXPECTATIVA (escrita e commitada ANTES de corrigir e rodar), para a cadeia
nos 88 com os manifestos em cache (mesmos desenhos, barras corrigidas):
  - controle TIC 142874476 (A e B) NAO se move: usa o erro da media para
    uma epoca unica, que e o uso certo. Se mover, a mudanca foi longe demais.
  - continuam 26 medidos de 88 (a escada e so TESS; a cobertura por
    temporada nao muda).
  - barras SuperWASP por temporada crescem por 1,0-1,73 (mediana 1,41);
    dP/dt muda pouco; sigma(dP/dt) cresce nos alvos cuja alavanca e SuperWASP.
  - curvatura sobrevivendo ao adversarial: 11 -> 10 (TIC 359552377 cai; as
    rodadas SuperWASP x1,34 e x1,65 concordam nisso); 9 ou 8 em D9.
  - mediana do chi2_red de 1,41 para 0,8-0,9 (x1,34 deu 0,86; x1,65 deu
    0,63); agregado 1,92 -> ~1,2; marcados p_gof < 0,05: 4 -> 2 ou 3.
  - `qual_barra_26`: fator SuperWASP 2,73 -> ~1,6 (intervalo ~1,0-2,7); TESS
    1,03 inalterado. O ~1,6 restante e do tamanho que s estimado com
    nu = 1-2 produz sozinho (E[z^2] = nu/(nu-2)); isso se testa DEPOIS, com
    o nulo re-estimando as barras das temporadas simuladas.
  - Sigma epsilon: 6,7 -> ~6,0.
Desvio em qualquer direcao vai para o humano antes de qualquer outra coisa.
"""
import glob
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cadeia_oc as co
import config
import download as dl
import epocas_142874476 as ep
import superwasp_cobertura as sc

CACHE = Path(config.CACHE) / "oc_2min"
SAIDA = config.DATA / "orquestra" / "oc_lote"
CONVENCAO = "parabola ponderada nos aglomerados TESS, avaliada em E do SuperWASP (cadeia_oc.ajustar)"


def _norm(t, f, q):
    ok = np.isfinite(t) & np.isfinite(f) & (np.asarray(q) == 0)
    t, f = np.asarray(t)[ok], np.asarray(f)[ok]
    return t, f / np.percentile(f, 95)


def ler_spoc_2min(path):
    """PDCSAP com QUALITY == 0, normalizado pelo percentil 95 - como `carregar`."""
    from astropy.io import fits
    with fits.open(path) as h:
        d = h[1].data
        setor = int(h[0].header["SECTOR"])
        cad = float(h[1].header.get("TIMEDEL", np.nan)) * 1440
        t, f = _norm(d["TIME"], d["PDCSAP_FLUX"], d["QUALITY"])
    if len(t) < 200:
        raise RuntimeError(f"FALHA [oc_lote]: {Path(path).name} com so {len(t)} pontos validos")
    return setor, cad, t, f


def manifestos():
    """Os 22 manifestos por setor mais o suplementar (setores intermediarios
    que o MAST lista e os 22 nao tinham - ver recuperar_14.py)."""
    # `manifest_s*` casava tambem com `manifest_suplementar...`: o suplementar
    # entrava DUAS vezes, cada setor dele virava duas "reducoes" identicas,
    # o espalhamento entre elas dava sigma = 0, o peso 1/0 e a escada devolvia
    # NaN - o "cannot convert float NaN to integer" dos 4 CVZ. Glob por digito
    # e deduplicacao por (tic, arquivo): as duas, porque uma so protegeria
    # contra este caso e nao contra o proximo.
    out = []
    for m in sorted(glob.glob(str(config.CATALOGS / "manifest_s[0-9]*_2min.parquet"))) +              sorted(glob.glob(str(config.CATALOGS / "manifest_suplementar*_2min.parquet"))):
        d = pd.read_parquet(m, columns=["tic", "filename", "url", "sector"])
        out.append(d)
    man = pd.concat(out, ignore_index=True)
    n0 = len(man)
    man = man.drop_duplicates(["tic", "filename"]).reset_index(drop=True)
    if len(man) != n0:
        print(f"  [manifestos] {n0 - len(man)} linhas duplicadas (tic, arquivo) removidas", flush=True)
    return man


def baixar_alvo(tic, man):
    CACHE.mkdir(parents=True, exist_ok=True)
    arqs = []
    for r in man[man.tic == tic].itertuples():
        dest = CACHE / r.filename
        if not dest.exists():
            ok, nb, err = dl.download_one(r.url, dest)
            if not ok:
                raise RuntimeError(f"FALHA [oc_lote]: download de {r.filename}: {err}")
        arqs.append((int(r.sector), dest))
    return sorted(arqs)


def epocas_2min(tic, P0, t14_h, man, t0_catalogo_btjd, span_min=None):
    """Uma epoca por setor SPOC 2-min. Ancora no setor MAIS RECENTE.

    A SEMENTE DA ANCORA E A EFEMERIDE DO CATALOGO, nao o minimo da curva: o
    docstring de `medir` registra o caso em que o minimo global caiu no
    SECUNDARIO (TGLC do alvo de controle) e a epoca saiu 1.008 min fora. Em
    binaria de contato, com eclipses quase iguais, isso e a regra e nao a
    excecao. O catalogo (Prsa BJD0 / Kostov T0-pri) diz qual e o primario.
    """
    arqs = baixar_alvo(tic, man)
    if not arqs:
        raise RuntimeError(f"FALHA [oc_lote]: TIC {tic} sem arquivo nos manifestos")
    curvas = [ler_spoc_2min(p) for _, p in arqs]
    s_anc, cad_anc, t_anc, f_anc = max(curvas, key=lambda c: c[2].max())
    if not np.isfinite(t0_catalogo_btjd):
        raise RuntimeError(f"FALHA [oc_lote]: TIC {tic} sem T0 de catalogo para semear o primario")
    t0_anc, sig_anc, prof_anc = ep.medir(t_anc, f_anc, P=P0, t14_h=t14_h, semente=float(t0_catalogo_btjd),
                                          span_min=span_min or max(120.0, 1.5 * t14_h * 60))
    linhas, descartados = [], []
    sp = span_min or max(120.0, 1.5 * t14_h * 60)
    for setor, cad, t, f in curvas:
        t0, sig, prof = ep.medir(t, f, P=P0, t14_h=t14_h, semente=t0_anc, span_min=sp)
        # A BARRA SAI DO PROPRIO DADO, nao so do delta-chi2 = 1. A formal e o
        # ruido de foton; manchas e O'Connell deslocam o minimo de ciclo para
        # ciclo e nao estao nela. O lote inteiro com barras formais deu 24 de
        # 31 alvos com |z| > 3 e chi2 reduzido da PARABOLA mediano 6,4 - a
        # assinatura de barra 3x pequena, nao de populacao curva. Mesma regra
        # de `epoca_por_subconjunto`: a maior entre formal e dispersao vale.
        # Aqui o bloco e a METADE do setor (o gap de meio de setor nao e um
        # corte de 30 d), e a dispersao e |t0_a - t0_b| / sqrt(2).
        meio = int(np.searchsorted(t, np.median(t)))
        disp = np.nan
        if meio > 100 and len(t) - meio > 100:
            ta, _, _ = ep.medir(t[:meio], f[:meio], P=P0, t14_h=t14_h, semente=t0, span_min=sp)
            tb, _, _ = ep.medir(t[meio:], f[meio:], P=P0, t14_h=t14_h, semente=t0, span_min=sp)
            ka, kb = np.round((ta - t0) / P0), np.round((tb - t0) / P0)
            disp = abs((ta - ka * P0) - (tb - kb * P0)) * 1440.0 / np.sqrt(2.0)
        sig_uso = max(float(sig), float(disp)) if np.isfinite(disp) else float(sig)
        if not (np.isfinite(t0) and np.isfinite(sig_uso) and sig_uso > 0):
            # a grade nao fechou delta-chi2 = 1 (barra indefinida) ou a epoca
            # nao e finita: o setor sai, e o motivo fica - antes ele entrava
            # como NaN e a escada morria com "cannot convert float NaN"
            descartados.append({"setor": setor, "motivo": "epoca ou barra nao finita",
                                "t0": float(t0), "sigma_formal_min": float(sig)})
            continue
        linhas.append({"setor": setor, "reducao": "SPOC", "t0_btjd": t0, "sigma_min": sig_uso,
                       "sigma_formal_min": float(sig), "sigma_metades_min": float(disp),
                       "prof_ajustada_ppm": prof * 1e6, "n": len(t), "cadencia_min": cad})
    if len(linhas) < 2:
        raise RuntimeError(f"so {len(linhas)} setor(es) com epoca finita "
                           f"(descartados: {[(d['setor'], d['motivo']) for d in descartados]})")
    out = pd.DataFrame(linhas).sort_values("setor")
    out.attrs["descartados"] = descartados
    return out



def medir_alvo_lote(tic, ep_df, P0, sourceid, ra, dec, t14_h):
    """A cadeia para um alvo do lote: AJUSTE CONJUNTO de todas as epocas.

    Difere do controle num ponto de DESENHO, nao de maquinaria: la, com 4
    aglomerados TESS, a parabola se ajusta no TESS e o SuperWASP a TESTA. Aqui
    32 dos 88 tem so 2 setores TESS - a parabola nem se ajusta neles. Entao o
    ajuste e conjunto, com CADA TEMPORADA do SuperWASP como epoca propria
    (barra = hypot(formal da temporada, 0,78 do vies), vies tirado), e o grau
    de liberdade e n_pontos - 3. E o que o dimensionamento contou quando
    somou temporadas a aglomerados; uma epoca SuperWASP unica daria 3 pontos
    para 3 parametros em 32 alvos, e "medido" em todos.

    Onde ha 3+ setores TESS, o teste predizer-e-conferir do controle tambem
    sai, por comparabilidade.
    """
    import superwasp as sw
    from scipy import stats
    ag = co.aglomerar(ep_df)
    k = len(ag)
    # A tolerancia de 3 min da escada nasceu com barras formais de 0,3-0,5 min
    # (6-10 sigma). Com a barra do proprio dado - mediana 3 min nos alvos de
    # recuperacao - 3 min e 1 sigma, e a escada rejeitava degraus corretos
    # (128790976 fecha a 10 min com residuo maximo 7,7 = 2,5 sigma). Escala
    # com a barra, sem baixar do piso: o controle (sigma ~0,4) fica em 3.
    tol = max(3.0, 3.0 * float(ag.sig_min.median()))
    esc = sw.escada_de_periodo(ag.rename(columns={"sig_min": "sig"})[["t0", "sig"]], P0, max_res_min=tol)
    if esc is None:
        raise RuntimeError(f"escada nao fechou no TESS (tolerancia {tol:.1f} min = 3 x barra mediana): "
                           "contagem de ciclos ambigua")
    P_t, T0_t, sP_t, E_t, _ = esc
    ag["E"] = E_t.astype(int)

    t, f, sig_mad = sw.carregar(sourceid, ra, dec)
    t_med = float(np.median(t))
    E_sw = int(np.round((t_med - (T0_t + 2457000.0)) / P_t))
    t0_prev = T0_t + 2457000.0 + E_sw * P_t
    t0g, sigf, sig_sub, subs, status = sw.epoca_por_subconjunto(
        t, f, P_t, t0_prev, t14_h=t14_h, n_sub=2, modo="temporada", exigir_cobertura=True)
    if not str(status).startswith("testado"):
        raise RuntimeError(f"cobertura por bloco: {status}")
    pontos = [{"fonte": "TESS s%d" % int(r.setor), "t0": float(r.t0), "sig_min": float(r.sig_min),
               "sig_formal_min": float(getattr(r, "sig_formal_min", float("nan"))),
               "sig_metades_min": float(getattr(r, "sig_metades_min", float("nan")))}
              for r in ag.itertuples()]
    # a dispersao ENTRE temporadas e a medida do ruido vermelho; cada temporada
    # entra com a maior entre a formal dela e essa dispersao. A dispersao de
    # UMA temporada e std(desvios) - `dispersao_entre_blocos_min`, posta em
    # cada bloco por `epoca_por_subconjunto` -, NAO `sig_sub` (std/sqrt n, o
    # erro da media, que e a barra da epoca combinada do controle). Ver o
    # docstring: era `sig_sub` aqui, barra sqrt(n) pequena em 49 de 61 epocas.
    disp_sw = float(subs[0]["dispersao_entre_blocos_min"]) if np.isfinite(sig_sub) else 0.0
    assert all(abs(sb["dispersao_entre_blocos_min"] - disp_sw) < 1e-9 for sb in subs) if np.isfinite(sig_sub) else True
    assert abs(disp_sw - float(sig_sub) * np.sqrt(len(subs))) < 1e-6 if np.isfinite(sig_sub) else True, (disp_sw, sig_sub, len(subs))
    for sb in subs:
        pontos.append({"fonte": "SuperWASP %.0f" % sb["t_medio"],
                       "t0": float(sb["t0"]) - 2457000.0 - co.VIES_CADEIA_MIN / 1440.0,
                       "sig_min": float(np.hypot(max(sb["sigma_formal_min"], disp_sw), co.VIES_CADEIA_SIG)),
                       "sig_formal_min": float(sb["sigma_formal_min"])})
    pt = pd.DataFrame(pontos).sort_values("t0").reset_index(drop=True)
    # A ESCADA NAO ATRAVESSA O VAO, E NAO DEVE. Ela resolve alias entre epocas
    # proximas exigindo residuo < 3 min numa efemeride LINEAR - e a 19 anos um
    # dP/dt real produz dezenas de minutos de O-C (o controle tem +22,7). Usa-
    # la no vao rejeitava como "alias nao resolvido" justamente os alvos com
    # sinal: 16 dos primeiros 28 do lote, inclusive com e_Per propagando a
    # 0,0005 ciclo. O controle nunca fez isso: la o E do SuperWASP e o
    # arredondamento com o P refinado no TESS, e a AMBIGUIDADE e a incerteza
    # do periodo vezes o numero de ciclos - nao o residuo, que e o sinal.
    if sP_t <= 0 or not np.isfinite(sP_t):
        raise RuntimeError(f"escada no TESS devolveu sP = {sP_t!r}: sem incerteza nao ha como julgar o vao")
    E_j = np.round((pt.t0.values - T0_t) / P_t)
    N = float(np.abs(E_j).max())
    amb = N * sP_t / P_t                 # fracao de ciclo que sP_t propaga ate a epoca mais distante
    if amb > 0.25:
        raise RuntimeError(f"contagem de ciclo ambigua no vao: sP {sP_t:.2e} d x {N:.0f} ciclos = "
                           f"{amb:.2f} ciclo (limite 0,25)")
    pt["E"] = E_j.astype(int)
    P_j, sP_j = P_t, sP_t
    sig_d = pt.sig_min.values / 1440.0
    cl, rl, chi2l, _ = co.ajustar(pt.E.values, pt.t0.values, sig_d, 1)
    n = len(pt)
    dof = n - 3
    res = {"tic": tic, "P_ref_d": P0, "P_escada_d": P_j, "sP_escada_d": sP_j,
           "n_pontos": n, "n_aglomerados_tess": k, "n_temporadas_swasp": len(subs), "dof": dof,
           "pontos": pt.assign(res_linear_min=rl * 1440).to_dict("records"),
           "linear": {"chi2": chi2l, "dof": n - 2},
           "superwasp": {"sourceid": sourceid, "npts": int(len(t)), "status_cobertura_por_bloco": status,
                         "barra_formal_min": float(sigf), "barra_entre_temporadas_min": float(sig_sub), "dispersao_entre_temporadas_min": disp_sw,
                         "temporadas": [{a: (round(v, 4) if isinstance(v, float) else v) for a, v in sb.items()} for sb in subs]}}
    if dof >= 1:
        cp, rp, chi2p, cov = co.ajustar(pt.E.values, pt.t0.values, sig_d, 2)
        Q = float(cp[0]); sQ = float(np.sqrt(cov[0, 0]))
        res["parabola"] = {"chi2": chi2p, "dof": dof, "Q_d_por_ciclo2": Q, "sQ": sQ,
                           "dPdt_s_por_ano": 2 * Q / P_j * 365.25 * 86400,
                           "sdPdt_s_por_ano": 2 * sQ / P_j * 365.25 * 86400,
                           "residuos_min": (rp * 1440).round(3).tolist()}
        res["delta_chi2"] = chi2l - chi2p
        res["p_curvatura"] = float(stats.chi2.sf(chi2l - chi2p, 1))
        res["classe"] = "testado"
        ag2 = pt.copy(); ag2["setor"] = range(n)
        res["adversarial"] = co.adversarial(ag2, P_j, float(pt.t0.iloc[-1]))
    else:
        res["classe"] = "medido, nao testado"
    if k >= 3:
        sg = ag.sig_min.values / 1440
        cl_t, _, _, _ = co.ajustar(ag.E.values, ag.t0.values, sg, 1)
        cp_t, _, _, _ = co.ajustar(ag.E.values, ag.t0.values, sg, 2)
        E_g = int(np.round((t0g - 2457000.0 - T0_t) / P_t))
        bt = float(np.hypot(max(sigf, sig_sub if np.isfinite(sig_sub) else 0.0), co.VIES_CADEIA_SIG))
        oc_l = (t0g - 2457000.0 - np.polyval(cl_t, E_g)) * 1440 - co.VIES_CADEIA_MIN
        oc_p = (t0g - 2457000.0 - np.polyval(cp_t, E_g)) * 1440 - co.VIES_CADEIA_MIN
        res["teste_predizer_conferir"] = {"E_swasp": E_g, "barra_total_min": bt,
                                          "residuo_no_linear_min": float(oc_l), "sigmas_linear": float(oc_l / bt),
                                          "residuo_na_parabola_min": float(oc_p), "sigmas_parabola": float(abs(oc_p) / bt)}
    return res

def rodar_um(r, man):
    t_ini = time.time()
    tic = int(r.tic)
    try:
        # T14 do catalogo maior que 30% do periodo nao e duracao de eclipse, e
        # largura de curva de contato mal definida; o trapezio degenera. Sai
        # como inaplicavel em vez de medir com forma errada (36,6 h no piloto)
        if not np.isfinite(r.t14_h) or r.t14_h / 24.0 > 0.30 * r.period_d:
            raise RuntimeError(f"T14 do catalogo inaplicavel: {r.t14_h:.2f} h para P = {r.period_d:.3f} d")
        ep_df = epocas_2min(tic, float(r.period_d), float(r.t14_h), man, float(r.t0_btjd))
        # profundidade que troca de sinal ou some entre setores = eclipse
        # errado (primario/secundario) ou semente fora: nao entra calado
        if (ep_df.prof_ajustada_ppm <= 0).any():
            raise RuntimeError(f"profundidade ajustada <= 0 em algum setor: "
                               f"{ep_df.prof_ajustada_ppm.round(0).tolist()}")
        res = medir_alvo_lote(tic, ep_df, float(r.period_d), r.sourceid,
                              float(r.ra), float(r.dec), float(r.t14_h))
        res["convencao_parabola"] = CONVENCAO
        res["status"] = "ok"
    except Exception as e:  # noqa: BLE001 - contado, nunca lido como resultado
        res = {"tic": tic, "status": "falha", "erro": str(e)[:200],
               "traceback": traceback.format_exc()[-600:]}
    res["segundos"] = time.time() - t_ini
    return res


def linha_resumo(res):
    """Uma linha por alvo, com as grandezas do desenho e os n ao lado."""
    if res.get("status") != "ok":
        return {"tic": res["tic"], "status": "falha", "erro": res.get("erro")}
    d = {"tic": res["tic"], "status": "ok", "P_d": res["P_ref_d"], "P_escada_d": res["P_escada_d"],
         "n_pontos": res["n_pontos"], "n_aglomerados_tess": res["n_aglomerados_tess"],
         "n_temporadas_swasp": res["n_temporadas_swasp"], "dof": res["dof"], "classe": res["classe"],
         "npts_swasp": res["superwasp"]["npts"],
         "cobertura_por_bloco": res["superwasp"]["status_cobertura_por_bloco"],
         "convencao": CONVENCAO}
    if "parabola" in res:
        d.update({"dPdt_s_por_ano": res["parabola"]["dPdt_s_por_ano"],
                  "sdPdt_s_por_ano": res["parabola"]["sdPdt_s_por_ano"],
                  "delta_chi2": res["delta_chi2"], "p_curvatura": res["p_curvatura"],
                  "p_adversarial": res["adversarial"]["p_adversarial"]})
    if "teste_predizer_conferir" in res:
        d.update({"oc_linear_min": res["teste_predizer_conferir"]["residuo_no_linear_min"],
                  "sigmas_linear": res["teste_predizer_conferir"]["sigmas_linear"],
                  "sigmas_parabola": res["teste_predizer_conferir"]["sigmas_parabola"]})
    return d


def controles():
    """A e B, e os dois tem que passar antes de qualquer alvo novo."""
    import comparar_oc as C
    exp = json.load(open(config.DATA / "orquestra" / "expectativa_oc_142874476_v2.json", encoding="utf-8"))
    ra, dec = sc._coord_do_id(exp["superwasp"]["sourceid"])

    # A: a cadeia compartilhada com as epocas ESPECIFICAS do controle
    ep_df, P0 = co.epocas_tess(142874476)
    res, ag, _ = co.medir_alvo_de(142874476, ep_df, P0, exp["superwasp"]["sourceid"], ra, dec)
    res["adversarial"] = co.adversarial(ag, P0, float(ag[ag.setor == 105].t0.iloc[0]))
    dA = C.comparar(exp, res)
    print(f"CONTROLE A (cadeia compartilhada, epocas especificas): "
          f"{'REPRODUZIU' if not dA else 'ESCALAR ' + str(dA)}")

    # B: o extrator generico no MESMO arquivo s105 do controle
    man = manifestos()
    arqs = baixar_alvo(142874476, man)
    setor, cad, t, f = ler_spoc_2min(arqs[0][1])
    esperado = exp["tess"]["s105"]["epoca_btjd"]
    t0, sig, prof = ep.medir(t, f, P=P0, t14_h=1.967, semente=esperado, span_min=120.0)
    dif = (t0 - esperado) * 1440
    okB = abs(dif) <= exp["tess"]["s105"]["tol_min"] and abs(sig - exp["tess"]["s105"]["sigma_min"]) <= 0.05
    print(f"CONTROLE B (leitor 2-min + medir no s105): epoca {t0:.6f} (esperado {esperado:.5f}, "
          f"dif {dif:+.3f} min), sigma {sig:.2f} (esperado 0,30) -> {'REPRODUZIU' if okB else 'ESCALAR'}")
    return (not dA) and okB


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--controles", action="store_true")
    ap.add_argument("--alvos", default=None, help="parquet com tic, ra, dec, period_d, t14_h, t0_btjd, sourceid")
    ap.add_argument("--limite", type=int, default=None)
    a = ap.parse_args()
    if a.controles or a.alvos is None:
        ok = controles()
        if not ok or a.alvos is None:
            sys.exit(0 if ok else 1)
    alvos = pd.read_parquet(a.alvos)
    if a.limite:
        alvos = alvos.head(a.limite)
    man = manifestos()
    SAIDA.mkdir(parents=True, exist_ok=True)
    resumos = []
    for i, r in enumerate(alvos.itertuples(), 1):
        res = rodar_um(r, man)
        (SAIDA / f"oc__{int(r.tic)}.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
        resumos.append(linha_resumo(res))
        lr = resumos[-1]
        print(f"[{i}/{len(alvos)}] TIC {int(r.tic)}: {lr['status']}"
              + ((f"  {lr['classe']}  n={lr['n_pontos']} dof={lr['dof']}"
                  + (f"  dP/dt {lr['dPdt_s_por_ano']:+.3f}+-{lr['sdPdt_s_por_ano']:.3f} s/ano"
                     f"  p {lr['p_curvatura']:.3f}  p_adv {lr['p_adversarial']:.3f}"
                     if "dPdt_s_por_ano" in lr else ""))
                 if lr["status"] == "ok" else f"  {lr['erro']}"), flush=True)
        pd.DataFrame(resumos).to_parquet(SAIDA / "resumo_oc_lote.parquet", index=False)
    print("gravado:", SAIDA / "resumo_oc_lote.parquet")


# Binaria de periodo curto varia tipicamente em 1e-2 a 1 s/ano. Acima disto e
# contagem de ciclo errada, nao astrofisica - a mesma logica do teste de Roche:
# "impossivel" e afirmacao forte, "cabe" e fraca. E a guarda que faltou quando
# a escada devolveu E embaralhado: nenhuma das guardas internas pegou, quem
# pegou foi a implausibilidade de -24.122 s/ano. Um embaralhamento entre duas
# epocas proximas daria -0,5 s/ano, plausivel para Algol com transferencia de
# massa, e teria entrado na distribuicao como achado.
DPDT_SANIDADE_S_ANO = 10.0


def consolidar(alvos_parquet=None):
    """Reconstroi o resumo a partir dos JSON por alvo, com a guarda de sanidade
    e a morfologia ao lado - re-executavel sem refazer nenhuma medida."""
    linhas = []
    for p in sorted(SAIDA.glob("oc__*.json")):
        res = json.loads(p.read_text(encoding="utf-8"))
        lr = linha_resumo(res)
        if lr.get("status") == "ok" and "dPdt_s_por_ano" in lr:
            if abs(lr["dPdt_s_por_ano"]) > DPDT_SANIDADE_S_ANO:
                lr["classe"] = "contagem de ciclo suspeita (|dP/dt| > %g s/ano)" % DPDT_SANIDADE_S_ANO
        linhas.append(lr)
    d = pd.DataFrame(linhas)
    if alvos_parquet:
        a = pd.read_parquet(alvos_parquet)[["tic", "morph", "period_d", "npts", "mag_swasp", "n_setores_2min", "t14_h"]]
        d = d.merge(a, on="tic", how="left", suffixes=("", "_alvo"))
    d.to_parquet(SAIDA / "resumo_oc_lote.parquet", index=False)
    return d


def relatar(d):
    from orquestra.contrato import Taxa
    out = []
    n = len(d)
    ok = d[d.status == "ok"]
    out.append(f"alvos: {n} | mediram: {Taxa.de(len(ok), n).texto()}")
    falhas = d[d.status != "ok"].get("erro", pd.Series(dtype=str)).fillna("?").str.slice(0, 48).value_counts()
    out.append("  nao mediram, por motivo declarado:")
    for k, v in falhas.items():
        out.append(f"    {v:3d}  {k}")
    if "morph" in d.columns:
        t14 = d[d.get("erro", pd.Series("", index=d.index)).fillna("").str.contains("T14")]
        if len(t14):
            out.append(f"  os que cairam por T14: morph mediana {t14.morph.median():.2f} (n={int(t14.morph.notna().sum())}) "
                       f"- {int((t14.morph >= 0.7).sum())} de contato, {int((t14.morph < 0.5).sum())} destacadas")
    for cl, g in ok.groupby("classe"):
        out.append(f"\n  {cl}: {len(g)}")
        if "dPdt_s_por_ano" in g.columns and g.dPdt_s_por_ano.notna().any():
            x = g.dPdt_s_por_ano.dropna()
            out.append(f"    dP/dt s/ano: mediana {x.median():+.4f} (n={len(x)}), p5 {x.quantile(.05):+.4f}, p95 {x.quantile(.95):+.4f}")
            out.append(f"    |dP/dt| mediana {x.abs().median():.4f} | barra mediana {g.sdPdt_s_por_ano.median():.4f}")
            if "p_curvatura" in g.columns:
                forte = g[(g.p_curvatura < 0.05) & (g.p_adversarial < 0.05)]
                out.append(f"    curvatura p < 0,05: {int((g.p_curvatura < 0.05).sum())} | "
                           f"que SOBREVIVE ao adversarial (p_adv < 0,05): {Taxa.de(len(forte), len(g)).texto()}")
                if len(forte):
                    out.append("    " + forte[["tic", "n_pontos", "dof", "dPdt_s_por_ano", "sdPdt_s_por_ano", "p_curvatura", "p_adversarial"]]
                               .to_string(index=False, float_format=lambda v: f"{v:.4f}").replace("\n", "\n    "))
    return "\n".join(out)


def relatar_distribuicoes(d):
    """As duas distribuicoes, juntas: dP/dt diz o que a populacao faz;
    dP/dt / sigma diz quantos sao significativos. Uma no lugar da outra
    mistura sinal com precisao, porque as barras variam muito entre alvos
    (n de 4 a 8 epocas, arquivos de qualidade diferente)."""
    ok = d[(d.status == "ok") & d.get("dPdt_s_por_ano", pd.Series(dtype=float)).notna()].copy()
    if not len(ok):
        return "sem alvos com dP/dt medido"
    ok["z"] = ok.dPdt_s_por_ano / ok.sdPdt_s_por_ano
    out = [f"\nDISTRIBUICOES (n={len(ok)} com dP/dt medido, todas as classes)"]
    out.append("  dP/dt (s/ano) - o que a populacao faz:")
    for q in (0.05, 0.25, 0.5, 0.75, 0.95):
        out.append(f"    p{int(q*100):<3d} {ok.dPdt_s_por_ano.quantile(q):+10.4f}")
    out.append(f"    |dP/dt| mediana {ok.dPdt_s_por_ano.abs().median():.4f} | barra mediana {ok.sdPdt_s_por_ano.median():.4f}")
    out.append("  dP/dt / sigma - quantos sao significativos:")
    for lim in (2, 3, 5):
        out.append(f"    |z| > {lim}: {int((ok.z.abs() > lim).sum()):3d} de {len(ok)}")
    out.append(f"  GUARDA DE SANIDADE: |dP/dt| > {DPDT_SANIDADE_S_ANO:g} s/ano em {int((ok.dPdt_s_por_ano.abs() > DPDT_SANIDADE_S_ANO).sum())} de {len(ok)}")
    sus = ok[ok.dPdt_s_por_ano.abs() > DPDT_SANIDADE_S_ANO]
    if len(sus):
        out.append("    " + sus[["tic", "n_pontos", "dof", "dPdt_s_por_ano", "sdPdt_s_por_ano", "P_d"]]
                   .to_string(index=False, float_format=lambda v: f"{v:.3f}").replace("\n", "\n    "))
    return "\n".join(out)
