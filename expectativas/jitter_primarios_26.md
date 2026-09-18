# Expectativa: `src/jitter_primarios_26.py` (Passo 1)

Escrita e commitada ANTES da primeira execucao. Nada rodado. Nada em
`cadeia_oc.py`, `oc_lote.py`, `superwasp.py` muda: o script IMPORTA
`oc_lote.epocas_2min` (o estimador do estado D: trapezio, semente do
catalogo, ancora no setor mais recente, barra max(formal, |a - b|/2)) e
`reinflar_tess_26.refazer` (refit a partir dos `pontos` gravados).

## O que o script faz

1. Manifesto `data/catalogs/pendente/completo26_2min.parquet` (799 curvas dos
   26). Desenho = setores de `oc_lote.manifestos()` por alvo (os que a cadeia
   usou); extras = os demais (715; 666 a baixar, ~1,3 GB, cache oc_2min).
2. Por alvo (paralelo, 6-8 processos; um JSON por alvo em
   `data/orquestra/oc_lote/jitter_26/`): `epocas_2min` no manifesto desenho +
   extras, primario (semente t0 catalogo) e secundario (semente + P/2). Controle
   embutido: os setores do desenho reproduzem os `pontos` do registro a
   < 0,05 min (como secundarios_26). Ancora = setor mais recente de todos.
3. Contagem de ciclos com P e T0 da ESCADA do desenho (T0 = t0 - E P do ultimo
   TESS do registro): E = round((t - T0)/P); falha se sP x |Delta E| >= 0,25,
   se a fase |t - T0 - E P| > 0,1 P (semente caiu no outro minimo), se
   profundidade <= 0 ou epoca/barra nao finita. Falhas registradas por setor.
4. sigma_j por alvo, primario e secundario: MV contra modelo linear local ao
   bloco TESS estendido (desenho + extras), mesmo perfil do Passo 0 (grade
   0-10 min, IC 68/95 por perfil). n setores extras e cobertura de lag (min,
   max, em dias) sobre todos os pares TESS. Informativo = >= 4 extras.
5. Secundario: sigma_j e a deriva d = t_sec - t_prim - P/2 (inclinacao
   ponderada em min/ano sobre todos os setores) comparada com a de 5.3.1
   (primeiro-ultimo setor do desenho).
6. Variograma agregado dos residuos primarios (informativos): pares (i,k) do
   mesmo alvo, tau = |t_k - t_i|; gama_norm(tau) = <(r_k - r_i)^2/(sig_i^2 + sig_k^2)>
   (esperanca 1 se branco com as barras do estado D) e gama_min2(tau) =
   1/2 <(r_k - r_i)^2 - (sig_i^2 + sig_k^2)> em min^2 (ruido de medida
   subtraido); bins 0-30 d, 30-100 d, 100 d-0,5 a, 0,5-1, 1-2, 2-3, 3-4, 4-5,
   5-7 a; bootstrap por alvo (2 000). Ajuste gama_min2 = c + q tau (tau em
   anos, pesos pelo bootstrap); q > 0 a 2 sigma (percentil 2,5 do bootstrap
   > 0) dispara o Passo 2. Convencao: passeio aleatorio com C = q_RW min(t_i,t_j)
   da 1/2 <Delta r^2> = c + q_RW tau / 2, logo q_RW = 2 q (declarado aqui e
   usado no Passo 2). Residuo de ajuste linear ao bloco absorve parte da
   deriva de baixa frequencia: q e um piso, dito assim.
7. Teste que substitui o teto de 5.3.2: barras TESS do DESENHO (os `pontos`,
   nada re-extraido) -> hypot(sigma, sigma_j) com sigma_j = MV do proprio alvo
   (informativos) e limite superior do IC 95% (nao informativos; 10 min se o
   perfil nao fecha). Cenario B declarado: limite superior do IC 95% para
   todos. `refazer`: p_curv, p_adv, p_gof, vao; quantos de D11/D9 sobrevivem
   e quais.
8. Saidas: `jitter_primarios_26.parquet` (uma linha por alvo e setor, primario
   e secundario), `jitter_primarios_26_resumo.parquet` (por alvo),
   `variograma_26.parquet` (bins + ajuste + bootstrap),
   `jitter_inflacao_26.parquet` (os 26 refeitos) e um JSON de resumo.

## Numeros esperados

- Downloads: 666 arquivos, 0 falhas (falha alta senao). Setores descartados
  (ciclo/fase/profundidade/nao finito): **< 5% dos 715**; alvos com ciclo
  ambiguo: **0-1** (sP x N do desenho <= 0,05 a 4 anos; os extras chegam a
  2025,0 -> < 0,1).
- Controle embutido: 84 de 84 epocas do desenho reproduzidas a < 0,05 min.
- Informativos: **21** (os do norte, 21-39 extras); nao informativos: **5**
  (115244268, 126763885, 126945917, 139256217, 144194304; 1-3 extras).
- sigma_j do PRIMARIO nos 21: **mediana 0,6-2,0 min** (abaixo da regra de
  parada de 3 min), faixa por alvo 0,2-4 min; 232634196 dentro do IC do
  Passo 0. Limite inferior do IC 95% > 0 em **>= 12 dos 21** (jitter exigido
  na maioria: barras de meia-setor de 0,3-2 min contra deriva de manchas).
- sigma_j do SECUNDARIO: **mediana 1-3 min**, razao sec/prim mediana
  **1,0-2,0** (minimo mais raso; e 5.3.1 achou 6 derivas diferenciais). Os
  alvos com deriva significativa em 5.3.1 (229476285, 229914020, 232634196,
  329289186, 390021728) ficam na metade superior de sigma_j secundario; a
  inclinacao de d sobre todos os setores tem o mesmo sinal da de 5.3.1 em
  >= 4 dos 5 e modulo compativel a 2 sigma.
- Variograma: gama_norm ~1,0-1,5 em tau < 30 d (dentro do setor a barra e da
  meia-setor); sobe para **1,5-3** em tau > 2 a (232634196 mostra +3,4 min/ano
  em 10 meses e residuos de 2020 vs 2024 fora da parabola). q > 0 a 2 sigma:
  **plausivel, ~50%**; se sim, Passo 2 roda. Se gama_norm ficar ~1 em todos os
  bins, o ruido entre setores e branco e o teto de 5.3.2 era so teto.
- Inflacao (cenario A, MV nos informativos, IC95 sup nos 5 do sul): barras TESS
  passam de 0,3-2,3 para ~1,0-3 min; sobrevivem **3-7 de D11** (entre o 1 de
  11 a 4,2 min de reinflar_tess_26 e os 11 sem inflacao); **232634196
  sobrevive**; **144194304 (sul, nao informativo, IC sup grande) cai**;
  D9 sobreviventes 2-6. Cenario B (IC95 sup em todos): 1-4 de 11.
- sP x N apos inflacao <= 0,1 em todos (vao nunca ambiguo).

Regra de parada: se a mediana de sigma_j dos primarios passar de 3 min, paro e
reporto ao fim do Passo 1. Desvio em qualquer direcao: relatorio, nao texto.
