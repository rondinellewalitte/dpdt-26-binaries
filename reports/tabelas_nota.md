# Tabelas da nota (geradas por src/tabela_nota.py)

## Tabela 2

| Reason | N | Nature |
|---|---|---|
| SuperWASP curve with < 500 valid points after cleaning | 22 | archive depth; structural |
| Phase coverage of the eclipse not testable in ≥ 1 season block (Section 3.3; a season whose fitted depth is ≤ 0 counts as not testable, Section 3.1) | 22 | archive sampling; recoverable with a second archive |
| SuperWASP light-curve download returned an empty file in four attempts (a different operation from the cone search of Table 1, step 4) | 4 | archive-side failure; recoverable if the archive serves the file |
| Cycle count between TESS sectors ambiguous (period ladder did not close) | 8 | local sector coverage; intermediate sectors were added from MAST and none was recovered (per-target outcome in the deposit) |
| No applicable eclipse duration in the catalogue (4 NaN; 1 contact-binary width; 1 W UMa) | 6 | catalogue gap; T14 was measured on the TESS curve for five and none was recovered (per-target outcome in the deposit) |

## Tabela 3

| TIC | Name | RA, Dec (J2000, deg) | Tmag | P (d) | N (TESS+SW) | d.o.f. | span (yr) | quadratic coefficient as dP/dt (s yr⁻¹) | χ²_red | p_gof | p_curv | p_adv | class | archival offset that defeats it (min) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 115244268 | TYC 7970-378-1 | 317.6428, -40.5239 | 10.9 | 2.0089 | 5 (3+2) | 2 | 19.9 | +0.0174 ± 0.0055 | 0.70 | 0.499 | 0.001 | 0.312 | no significant curvature (fails adversarial) | — |
| 126763885 | UCAC4 219-184178 | 316.2145, -46.3262 | 12.4 | 4.3574 | 4 (2+2) | 1 | 19.9 | -0.0501 ± 0.0422 | 1.20 | 0.273 | 0.234 | 0.623 | no significant curvature | — |
| 126945917 | HD 202042 | 318.7324, -43.3954 | 9.1 | 1.3411 | 6 (4+2) | 3 | 19.9 | +0.0099 ± 0.0035 | 0.53 | 0.664 | 0.004 | 0.665 | no significant curvature (fails adversarial) | — |
| 139256217 | TYC 8448-418-1 | 347.1903, -46.1102 | 11.0 | 2.1992 | 4 (2+2) | 1 | 19.9 | -0.0274 ± 0.0408 | 1.07 | 0.300 | 0.501 | 0.861 | no significant curvature | — |
| 144194304 | V Gru | 327.9727, -42.3732 | 9.2 | 0.4835 | 5 (2+3) | 2 | 19.9 | +0.0173 ± 0.0014 | 0.53 | 0.589 | < 10⁻³ | < 10⁻³ | curvature (parabola adequate) | 7 |
| 198388252 | [GGM2006] 2908282 | 258.1953, +57.0529 | 11.6 | 0.3952 | 7 (4+3) | 4 | 19.6 | +0.0141 ± 0.0016 | 7.66 | < 10⁻³ | < 10⁻³ | < 10⁻³ | curvature (parabola inadequate) | 7 |
| 198408416 | NSVS 2910034 | 258.9133, +58.8645 | 12.1 | 0.3647 | 7 (4+3) | 4 | 19.5 | +0.0099 ± 0.0011 | 0.64 | 0.631 | < 10⁻³ | < 10⁻³ | curvature (parabola adequate) | 5 |
| 199688409 | BPS BS 16080-0095 | 253.0515, +57.7255 | 12.0 | 0.8683 | 7 (4+3) | 4 | 17.6 | -0.0056 ± 0.0026 | 0.22 | 0.925 | 0.032 | 0.677 | no significant curvature (fails adversarial) | — |
| 207496824 | V353 Dra | 246.9546, +58.8398 | 10.8 | 0.8027 | 4 (2+2) | 1 | 16.7 | +0.0190 ± 0.0072 | 2.36 | 0.124 | 0.008 | 0.386 | no significant curvature (fails adversarial) | — |
| 224605072 | [DCO2008] T-Dra0-00919 | 255.5992, +58.8700 | 11.4 | 0.6649 | 7 (4+3) | 4 | 17.6 | -0.0510 ± 0.0893 | 0.68 | 0.606 | 0.568 | 0.978 | no significant curvature | — |
| 229461186 | BL Dra | 295.1026, +60.9203 | 10.9 | 0.4027 | 5 (3+2) | 2 | 15.7 | +0.0087 ± 0.0036 | 0.96 | 0.381 | 0.014 | 0.332 | no significant curvature (fails adversarial) | — |
| 229476285 | UCAC4 753-059681 | 295.5345, +60.4830 | 12.9 | 1.0965 | 5 (3+2) | 2 | 15.7 | +0.2189 ± 0.3975 | 0.79 | 0.455 | 0.582 | 0.762 | no significant curvature | — |
| 229687624 | CV Dra | 262.9885, +57.1784 | 9.0 | 0.6176 | 7 (4+3) | 4 | 19.6 | -0.0013 ± 0.0021 | 0.56 | 0.691 | 0.542 | 0.686 | no significant curvature | — |
| 229914020 | 1RXS J191035.4+635611 | 287.6436, +63.9373 | 11.9 | 0.6760 | 5 (3+2) | 2 | 15.7 | +0.0133 ± 0.0023 | 5.42 | 0.004 | < 10⁻³ | < 10⁻³ | curvature (parabola inadequate) | 2.5 |
| 230386284 | StKM 1-1676 | 285.8220, +63.9928 | 9.0 | 0.3415 | 5 (3+2) | 2 | 15.7 | +0.0081 ± 0.0006 | 0.07 | 0.928 | < 10⁻³ | < 10⁻³ | curvature (parabola adequate) | 5 |
| 232634196 | TYC 3906-474-1 | 269.4385, +56.1960 | 9.9 | 1.2511 | 7 (4+3) | 4 | 19.6 | -0.2181 ± 0.0207 | 0.72 | 0.578 | < 10⁻³ | < 10⁻³ | curvature (parabola adequate) | > 10 |
| 237116051 | V391 Dra | 269.7931, +58.7165 | 10.2 | 2.4629 | 4 (2+2) | 1 | 19.6 | -0.0099 ± 0.0087 | 13.16 | < 10⁻³ | 0.255 | 0.582 | no significant curvature (parabola inadequate) | — |
| 243352373 | V584 Dra | 290.2686, +56.3283 | 10.2 | 0.3343 | 4 (2+2) | 1 | 16.7 | -0.0042 ± 0.0026 | 0.74 | 0.389 | 0.113 | 0.839 | no significant curvature | — |
| 329246824 | V564 Dra | 263.8847, +57.8025 | 11.0 | 0.5877 | 7 (4+3) | 4 | 19.6 | -0.0204 ± 0.0026 | 0.18 | 0.949 | < 10⁻³ | < 10⁻³ | curvature (parabola adequate) | 7 |
| 329248002 | NP Dra | 263.8179, +55.0034 | 8.8 | 3.1089 | 6 (4+2) | 3 | 16.7 | +0.0275 ± 0.0251 | 0.81 | 0.487 | 0.273 | 0.778 | no significant curvature | — |
| 329289186 | ATO J269.1883+52.8656 | 269.1883, +52.8657 | 13.8 | 0.5525 | 5 (2+3) | 2 | 19.6 | -0.0063 ± 0.0030 | 0.55 | 0.576 | 0.039 | 0.966 | no significant curvature (fails adversarial) | — |
| 359552377 | TYC 3913-792-1 | 279.3983, +57.7443 | 10.9 | 3.4391 | 5 (3+2) | 2 | 16.7 | -0.2498 ± 0.0775 | 0.67 | 0.512 | 0.001 | 0.226 | no significant curvature (fails adversarial) | — |
| 377253090 | NSVS 3056525 | 286.1923, +58.9513 | 11.8 | 0.4245 | 6 (4+2) | 3 | 16.7 | -0.0061 ± 0.0010 | 0.52 | 0.670 | < 10⁻³ | < 10⁻³ | curvature (parabola adequate) | 2 |
| 390021728 | [GGM2006] 3034154 | 277.7934, +56.4175 | 11.5 | 0.3020 | 6 (4+2) | 3 | 16.6 | -0.0052 ± 0.0010 | 1.15 | 0.328 | < 10⁻³ | 0.004 | curvature (parabola adequate) | 2 |
| 392536812 | NSVS 2953494 | 280.4674, +58.0565 | 11.2 | 0.5910 | 6 (4+2) | 3 | 16.7 | +0.0536 ± 0.0029 | 0.24 | 0.866 | < 10⁻³ | < 10⁻³ | curvature (parabola adequate) | > 10 |
| 424461577 | V527 Dra | 279.7343, +60.3982 | 10.2 | 0.7448 | 6 (4+2) | 3 | 15.7 | -0.0368 ± 0.0033 | 2.04 | 0.106 | < 10⁻³ | < 10⁻³ | curvature (parabola adequate) | 10 |

## Tabela 4

| Quantity | n = 26 (all) | n = 23 (quadratic fit adequate, p_gof ≥ 0.05) | n = 16 (adequate and σ(dP/dt) < 0.01 s yr⁻¹) |
|---|---|---|---|
| dP/dt, median | -0.0027 s yr⁻¹ | -0.0042 s yr⁻¹ | +0.0034 s yr⁻¹ |
| \|dP/dt\|, median | 0.0157 s yr⁻¹ | 0.0174 s yr⁻¹ | 0.0093 s yr⁻¹ |
| dP/dt, interquartile range | [-0.0178, +0.0139] s yr⁻¹ | [-0.0239, +0.0136] s yr⁻¹ | [-0.0057, +0.0117] s yr⁻¹ |
| dP/dt, minimum and maximum | [-0.250, +0.219] s yr⁻¹ | [-0.250, +0.219] s yr⁻¹ | [-0.037, +0.054] s yr⁻¹ |
| Sign of dP/dt: positive / negative (two-sided sign test p) | 12 / 14 (p = 0.85) | 10 / 13 (p = 0.68) | 8 / 8 (p = 1.00) |
| σ(dP/dt), median | 0.0032 s yr⁻¹ | 0.0033 s yr⁻¹ | 0.0026 s yr⁻¹ |
| χ²_red of quadratic fit, median (25th–75th) | 0.71 (0.53–1.13) | 0.68 (0.53–0.89) | 0.56 (0.45–0.80) |
| \|dP/dt\|/σ > 2 / > 3 / > 5 | 18 / 13 / 11 | 16 / 11 / 9 | 14 / 9 / 8 |
| Curvature with p < 0.05 | 18 | 16 | 14 |
| Curvature surviving the adversarial test | 11 of 26 — 6 > 0, 5 < 0 | 9 of 23 — 4 > 0, 5 < 0 | 8 of 16 — 4 > 0, 4 < 0 |
| Same, over targets with detection efficiency ε > 0 (Section 5.4) | 11 of 24 | 9 of 21 | — |
| Same, with TESS bars inflated by 1.0 / 2.6 min yr⁻¹ (Section 5.3, sensitivity ceiling) | 1 / 0 of 26 | 1 / 0 of 23 | 0 / 0 of 16 |
| Measured but not tested (d.o.f. < 1) | 0 | 0 | 0 |
| \|dP/dt\| > 10 s yr⁻¹ (cycle-count guard) | 0 | 0 | 0 |

## Tabela 6

| TIC | Name | block, admitted (s yr⁻¹) | N | z (adm.) | class | MDD | pooled GP (s yr⁻¹, χ²/ν) | Δ per 2.14 min |
|---|---|---|---|---|---|---|---|---|
| 115244268 | TYC 7970-378-1 | — | — | — | not in D9 | — | — | 0.0053 |
| 126763885 | UCAC4 219-184178 | — | — | — | not in D9 | — | — | 0.0117 |
| 126945917 | HD 202042 | — | — | — | not in D9 | — | — | 0.0036 |
| 139256217 | TYC 8448-418-1 | — | — | — | not in D9 | — | — | 0.0059 |
| 144194304 | V Gru | — | 4 | — | not testable | — | — | 0.0013 |
| 198388252 | [GGM2006] 2908282 | — | — | — | not in D9 | — | — | 0.0011 |
| 198408416 | NSVS 2910034 | +0.0468 ± 0.0107 | 39 | +3.4 to +6.1 | contradicted | 0.0321 | +0.0452 ± 0.0120 (0.33) | 0.0011 |
| 199688409 | BPS BS 16080-0095 | — | — | — | not in D9 | — | — | 0.0030 |
| 207496824 | V353 Dra | — | — | — | not in D9 | — | — | 0.0028 |
| 224605072 | [DCO2008] T-Dra0-00919 | — | — | — | not in D9 | — | — | 0.0020 |
| 229461186 | BL Dra | — | — | — | not in D9 | — | — | 0.0016 |
| 229476285 | UCAC4 753-059681 | — | — | — | not in D9 | — | — | 0.0026 |
| 229687624 | CV Dra | — | — | — | not in D9 | — | — | 0.0017 |
| 229914020 | 1RXS J191035.4+635611 | — | — | — | not in D9 | — | — | 0.0026 |
| 230386284 | StKM 1-1676 | +0.0091 ± 0.0008 | 38 | +1.0 | not contradicted | 0.0024 | +0.0088 ± 0.0110 (0.05) | 0.0013 |
| 232634196 | TYC 3906-474-1 | -0.6485 ± 0.1328 | 43 | -9.5 to -3.2 | contradicted | 0.3983 | -0.7022 ± 0.0473 (2.16) | 0.0036 |
| 237116051 | V391 Dra | — | — | — | not in D9 | — | — | 0.0057 |
| 243352373 | V584 Dra | — | — | — | not in D9 | — | — | 0.0012 |
| 329246824 | V564 Dra | -0.0001 ± 0.0066 | 42 | +2.9 | ambiguous | 0.0197 | -0.0036 ± 0.0196 (0.30) | 0.0017 |
| 329248002 | NP Dra | — | — | — | not in D9 | — | — | 0.0111 |
| 329289186 | ATO J269.1883+52.8656 | — | — | — | not in D9 | — | — | 0.0016 |
| 359552377 | TYC 3913-792-1 | — | — | — | not in D9 | — | — | 0.0123 |
| 377253090 | NSVS 3056525 | -0.0050 ± 0.0049 | 40 | -0.0 to +0.2 | not contradicted | 0.0146 | -0.0022 ± 0.0139 (0.56) | 0.0015 |
| 390021728 | [GGM2006] 3034154 | +0.0286 ± 0.0088 | 41 | +3.8 to +5.5 | contradicted | 0.0263 | +0.0272 ± 0.0103 (1.98) | 0.0011 |
| 392536812 | NSVS 2953494 | +0.0209 ± 0.0074 | 42 | -4.1 | contradicted | 0.0223 | +0.0184 ± 0.0216 (0.45) | 0.0021 |
| 424461577 | V527 Dra | -0.0090 ± 0.0564 | 40 | — | published LTTE | — | +0.0513 ± 0.0264 (6.37) | 0.0029 |

## Tabela 7

| TIC | class (window test) | 0.75 | 1 | 1.5 | 2 | 2.5 | 3 | 4 | 5 | 7 | 10 | offset that defeats it |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 144194304 | — | — | — | — | — | — | — | — | — | —* | —* | 7 |
| 198388252 | — | — | — | — | — | — | — | — | — | —* | —* | 7 |
| 198408416 | contradicted | 3.4 | 3.4 | 3.5 | 3.5 | 3.5 | 3.5 | 3.5 | 3.5* | 3.6* | 3.6* | 5 |
| 229914020 | — | — | — | — | — | —* | —* | —* | —* | —* | —* | 2.5 |
| 230386284 | not contradicted | 1.0 | 0.9 | 0.8 | 0.6 | 0.5 | 0.5 | 0.4 | 0.3* | 0.2* | 0.1* | 5 |
| 232634196 | contradicted | 3.2 | 3.2 | 3.2 | 3.2 | 3.2 | 3.2 | 3.2 | 3.2 | 3.2 | 3.2 | > 10 |
| 329246824 | ambiguous | 2.9 | 2.9 | 2.8 | 2.8 | 2.7 | 2.7 | 2.6 | 2.4 | 2.2* | 1.9* | 7 |
| 377253090 | not contradicted | 0.0 | 0.0 | 0.0 | 0.0* | 0.0* | 0.0* | 0.0* | 0.1* | 0.1* | 0.2* | 2 |
| 390021728 | contradicted | 3.8 | 3.8 | 3.8 | 3.8* | 3.8* | 3.7* | 3.7* | 3.6* | 3.4* | 3.1* | 2 |
| 392536812 | contradicted | 4.1 | 4.1 | 4.0 | 4.0 | 3.9 | 3.8 | 3.6 | 3.4 | 3.0 | 2.5 | > 10 |
| 424461577 | published LTTE | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5* | 10 |

## Resumo

- Flag p_gof < 0.05: 3 of 26 flagged: 198388252 (χ²_red 7.66, 4 dof, p 0.000), 229914020 (χ²_red 5.42, 2 dof, p 0.004), 237116051 (χ²_red 13.16, 1 dof, p 0.000)
-   of which with significant surviving curvature: 2; without: 1
- Old rule χ²_red ≥ 2: 5 of 26 flagged: 198388252 (4 dof, p_gof 0.000), 207496824 (1 dof, p_gof 0.124), 229914020 (2 dof, p_gof 0.004), 237116051 (1 dof, p_gof 0.000), 424461577 (3 dof, p_gof 0.106)
-   implied false-flag probability of χ²_red ≥ 2 by dof: 1 dof: 0.157, 2 dof: 0.135, 3 dof: 0.112, 4 dof: 0.092
-   chance that a flag of 5/26 hits one given target by coincidence: 19%
- Marginal curvature (0.01 < p_curv < 0.05): 199688409 (p 0.032, p_adv 0.677), 229461186 (p 0.014, p_adv 0.332), 329289186 (p 0.039, p_adv 0.966) — surviving adversarial: 0
- Constellations: Dra 21, Mic 2, Gru 2, Ind 1
- D11 (curvature surviving adversarial, all 26): 144194304, 198388252, 198408416, 229914020, 230386284, 232634196, 329246824, 377253090, 390021728, 392536812, 424461577
- D9 (D11 with adequate fit): 144194304, 198408416, 230386284, 232634196, 329246824, 377253090, 390021728, 392536812, 424461577
- D11 minus D9 (surviving but flagged): 198388252, 229914020
- Flagged without significant curvature: 237116051
- Column 3 of Table 4 (adequate and sigma < 0.01): n=16, surviving 8: 144194304, 198408416, 230386284, 329246824, 377253090, 390021728, 392536812, 424461577
- Sign test n=23: 10+/13- p=0.68; n=26: 12+/14- p=0.85
- Surviving curvature: 9 of 23; 11 of 26
- Surviving curvature over targets with epsilon > 0: all 11 of 24 | adequate 9 of 21 | col3 8 of 16
- Median |dP/dt| by partition: all 0.0157 (n=26) | adequate 0.0174 (n=23) | sigma<0.01 0.0093 (n=16); 0.0157 s/yr = 1.82e-07 d/yr
- Sign test sigma<0.01: 8+/8- p=1.00; surviving 8/16
