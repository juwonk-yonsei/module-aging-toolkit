# 모듈–유전자 대응 감사 보고서

결정 규칙에 따라 **단계 3에서 중단**했다. 공식 Supplementary Table 5B의 모듈별 유전자 목록이 집중도 검정에서 무작위 수준이었으므로 저장소 코드는 수정하지 않았고, 재학습과 영향 범위 목록도 만들지 않았다.

## 1. 요약

오류는 확인됐다. 다만 위치는 이 저장소의 병합 코드가 아니다. 파이프라인이 쓰는 23개 모듈 목록은 공식 Supplementary Table 5B의 비영(非零) 계수 유전자와 같거나 그 상위집합이고, 그 공식 목록 자체가 유전자 계열 집중도 검정에서 무작위 수준이다(모든 계열 p > 0.4). 수정하지 않았다. 공개 계수 재현과 재학습은 수행하지 않았다.

## 2. 환경

- 시작 커밋과 끝 커밋은 같다. `a3dd3098b78f8eee8ea8a0a2d363cf87f32d0e96` (v0.2.0 이후, 메시지 `Add Zenodo DOI to README and CITATION.cff`). 브랜치 `audit/module-mapping`을 이 커밋에서 만들었고, 수정 커밋은 만들지 않았다. push는 하지 않았다.
- 이 브랜치의 작업 트리에는 감사와 무관한 기존 수정(원고 리비전 산출물)이 이미 있었다. 그 파일들은 읽기만 했고 더 고치지 않았다.
- Python 3.11.7, numpy 1.26.4, pandas 2.1.4, scikit-learn 1.2.2, scipy 1.11.4, openpyxl 3.0.10. 시드 20260808은 계수 표본 추출에만 썼다. 집중도 검정 스크립트의 순열 시드는 부록 A 기본값 0이다.
- 실행한 명령, 모두 `module-aging-toolkit` 루트에서:
  - `git checkout -b audit/module-mapping`
  - `curl`로 MOESM6, MOESM7, MOESM8 다운로드 (`results/audit_module_mapping/source/`)
  - `python scripts/audit/build_steps_1_to_3.py`
  - `python scripts/audit/module_coherence_audit.py results/audit_module_mapping/pipeline_modules.csv module_colour mouse_symbol annotation`
  - 같은 스크립트를 `s3a_mortality_modules.csv`, `original_modules.csv`에 각각 실행

지시서와 달랐던 경로: 툴킷 `data/`는 비어 있다. `train_module_clocks.py`가 읽는 `data/supp/module_membership_rodent.csv`와 학습 행렬은 실제로는 `/S4T/jkang/Research/recent/B_module_aging/data/`에 있다. 공식 xlsx의 md5는 그 로컬 복사본과 같다.

## 3. 입력 파일

| 역할 | 경로 또는 URL | md5 | 시트 / 헤더 |
|---|---|---|---|
| (c) Table 4, MOESM6 | https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-026-10542-3/MediaObjects/41586_2026_10542_MOESM6_ESM.xlsx | `8ec1cc76c1ac3d77450386d7fdad1f6f` | `(B) Module rodent clocks`. 헤더 1행. `Outcome`은 아래로 채워져 있다. |
| (d)(b) Table 5, MOESM7 | 같은 접두 + `MOESM7_ESM.xlsx` | `990a4d6cc6a9aaa892606f47fa2e1731` | `(A) Composite clocks` 헤더 1행. `(B) Module rodent clocks`는 4행 헤더(outcome, module, annotation, Entrez/symbol) 뒤 절편 행, 그 다음 2,141유전자. |
| (a)에 해당하는 Table 6, MOESM8 | 같은 접두 + `MOESM8_ESM.xlsx` | `e44f1a3dfbfa6f510adc31d8806cd777` | `(B) Rodent module dictionary`. 헤더 1행. 열은 Module ID, Functional annotation, Module size뿐이고 **유전자 열이 없다.** |
| 파이프라인 membership | `B_module_aging/data/supp/module_membership_rodent.csv` | `af977d9f8cb509cd4d7be92437ee59ab` | 열 `entrez, symbol, module, annotation`. 1,995행. |
| 학습 행렬 | Zenodo 18763485 `Expression_data_relative_rodents_Scaled.csv`, 로컬 복사본 동일 | `c283869cf69e35afa200c6a6fc281519` | 유전자 18,286 × 샘플 4,539. |
| 주석 | 같은 Zenodo 기록의 `Data_annotation_relative_rodents.xlsx` | `bb56d389a63c76f02c4e1b98319d6bd8` | 이번 단계에서 내용을 다시 읽지는 않았다. |
| S3a | `supplementary_revision/Supplementary_Tables.xlsx` 시트 `S3a_Module_gene_weights` | 워크북 전체가 작업 트리에서 수정돼 있으나, 유전자–모듈 쌍은 HEAD와 동일하다(3,990쌍, 불일치 0). | 헤더 1행. |

28개 모듈의 유전자 membership 표는 보충자료에 없다. 유전자 단위 원 목록으로 쓴 것은 Table 5B에서 계수가 0이 아닌 유전자이다. 이 표는 칸이 모두 숫자라서, 계수가 0인 멤버와 비멤버를 구분할 수 없다.

## 4. 단계 1 결과

덤프는 `results/audit_module_mapping/pipeline_modules.csv`이다. `train_module_clocks.py`가 읽는 CSV를 `module_colour, annotation, mouse_entrez, mouse_symbol`로 바꿨다. 1,995유전자, 23모듈.

S3a(mortality와 chronological 각각)와 유전자–모듈 쌍이 완전히 같다. `results/audit_module_mapping/s3a_vs_pipeline.csv`: 각 outcome에서 pipeline 1,995, S3a 1,995, 교집합 1,995, 한쪽에만 있는 쌍 0, `identical_pairs=True`, annotation 불일치 0, symbol 불일치 0. S3a를 내보내는 단계에서 생긴 어긋남은 없다. S3a 유전자는 학습된 pkl의 `genes`에서 오므로, 저장된 clock의 유전자 목록도 이 CSV와 같다.

집중도 검정(`pipeline_modules_coherence.csv`, `s3a_mortality_modules_coherence.csv`):

| family | n | top module | top share | null median | p | modules spanned |
|---|---:|---|---:|---:|---:|---:|
| ribosomal (Rpl/Rps) | 66 | turquoise (Innate immunity) | 0.182 | 0.197 | 0.8022 | 19 |
| OXPHOS | 80 | brown (Chromatin modification 1) | 0.212 | 0.200 | 0.4156 | 17 |
| interferon-stimulated | 45 | blue (Respiration) | 0.156 | 0.200 | 0.9650 | 18 |
| cell cycle | 39 | blue (Respiration) | 0.128 | 0.205 | 0.9990 | 17 |
| collagen | 12 | blue (Respiration) | 0.250 | 0.250 | 0.7862 | 7 |
| heat shock (Hsp/Dnaj) | 19 | brown (Chromatin modification 1) | 0.316 | 0.211 | 0.2048 | 7 |

S3a mortality 표는 위와 같은 n, top module, top share, null median, modules spanned이다. p만 순열 난수로 0.01 이내에서 다르다(`s3a_mortality_modules_coherence.csv`).

지시서 1절 기준선과 비교하면 n, top module, top share, null median, modules spanned는 같다. p는 소수 둘째 자리에서 다르다(가장 큰 차이는 heat shock 0.2048 대 기준선 0.18). 점추정은 같고 결론도 같다. p < 0.001인 계열은 없다.

## 5. 단계 2 결과

공식 Table 5B 유전자 목록의 집중도 검정은 `original_modules_coherence.csv`이다. 역시 모두 무작위 수준이다.

| family | n | top module | top share | null median | p | modules spanned |
|---|---:|---|---:|---:|---:|---:|
| ribosomal | 64 | turquoise (Innate immunity) | 0.188 | 0.203 | 0.7942 | 19 |
| OXPHOS | 77 | brown (Chromatin modification 1) | 0.182 | 0.208 | 0.7922 | 17 |
| interferon-stimulated | 43 | blue (Respiration) | 0.163 | 0.209 | 0.9441 | 18 |
| cell cycle | 37 | blue (Respiration) | 0.135 | 0.216 | 0.9970 | 17 |
| collagen | 12 | blue (Respiration) | 0.250 | 0.250 | 0.7782 | 7 |
| heat shock | 18 | brown (Chromatin modification 1) | 0.278 | 0.222 | 0.4416 | 7 |

원 모듈 수와 크기(`original_module_dictionary_28.csv`, Table 6B): 28개, 크기 합 2,141. 범위는 12(royalblue, salmon)에서 357(blue)이다. 본문의 “30–630”과는 숫자가 다르다. Table 6B에 적힌 크기가 이 파일의 모듈 크기이다.

23개 retained 모듈(`module_size_reconciliation.csv`): Table 4B에 clock이 없는 purple(67, 주석 `-`)와 royalblue(12, 주석 `-`), 그리고 두 outcome 평균 Pearson r이 0.3 미만인 white(0.283), lightcyan(0.239), darkturquoise(0.253)를 뺀 23개이다. 이 23개의 colour–주석 대응은 Table 6B와 파이프라인 CSV가 같다. Table 5B 헤더의 colour–주석도 같다. 저장소에만 라벨이 뒤바뀐 것이 아니다.

학습 행렬 18,286유전자 안에 Table 5B의 2,141유전자가 모두 있다(`training_matrix_overlap.csv`: `n_query_in_matrix_index=2141`).

1,995 대 2,141 (`module_size_reconciliation.csv`):

- 2,141은 28개 모듈 크기의 합이고, global model의 비영 계수 수와 같다(Table 4B, All module genes).
- 위에 적은 5개 제외 모듈의 크기 합은 134이다. 2,141 − 134 = 2,007. 이것이 Table 6 기준 23개 모듈의 유전자 수이다.
- 파이프라인은 1,995개로, 2,007보다 12개 적다. 부족한 모듈은 brown 9개(255 대 246), darkgreen 2개(36 대 34), darkorange 1개(17 대 16)이다.
- Table 5B에서 23개 모듈의 비영 계수는 1,922개이다. 나머지 219개(2,141 − 1,922)는 global 열에만 계수가 있고 개별 모듈 열은 0이다. 파이프라인 1,995는 이 1,922의 상위집합이다(단계 3).

Table 6에는 유전자 목록이 없으므로, 2,007이라는 수는 크기 합이지 유전자 ID로 확인한 합집합이 아니다.

## 6. 단계 3 결과

Jaccard는 `jaccard_matrix.csv`, 모듈별 요약은 `jaccard_summary.csv`이다. 행은 파이프라인 모듈, 열은 Table 5B 비영 유전자 모듈이다.

| repo module | n | best original | best Jaccard | same-name Jaccard | 분류 |
|---|---:|---|---:|---:|---|
| black, blue, darkred, darkgrey, greenyellow, green, lightgreen, lightyellow, magenta, orange, pink, red, salmon, skyblue, tan, turquoise, yellow (17개) | 표 참조 | 같은 이름 | 1 | 1 | (i) 정확 일치 |
| brown | 246 | brown | 0.7764 | 0.7764 | 같은 이름의 상위집합. (ii)도 (iii)도 아님 |
| darkgreen | 34 | darkgreen | 0.8235 | 0.8235 | 같음 |
| grey60 | 52 | grey60 | 0.8846 | 0.8846 | 같음 |
| midnightblue | 52 | midnightblue | 0.9231 | 0.9231 | 같음 |
| darkorange | 16 | darkorange | 0.9375 | 0.9375 | 같음 |
| cyan | 20 | cyan | 0.9500 | 0.9500 | 같음 |

(ii) 라벨 뒤바뀜은 없다. 최댓값이 항상 같은 이름에서 나온다. (iii)에 해당하는, 어떤 원 모듈과도 Jaccard가 낮은 경우는 없다. 17개는 Table 5B 지지집합과 정확히 같고, 6개는 그 지지집합을 포함하고 계수가 0인 유전자를 더 가지고 있다. 추가분은 73유전자이다(246−191, 34−28, 52−46, 52−48, 16−15, 20−19).

계수 대조(`coefficient_spotcheck_summary.csv`): S3a의 `coef_published_SuppTable5B`와 Table 5B 값을 outcome × module × Entrez로 맞추면 3,827쌍이 모두 일치한다(불일치 0, max abs diff 0). S3a에만 있는 쌍은 163개(두 outcome 합)이고, Table 5B에만 있는 쌍은 0개이다. 시드 20260808로 뽑은 50쌍은 `coefficient_spotcheck_50.csv`에 있다. 유전자와 계수 값의 연결은 공식 표와 같다. 어긋난 것은 그 유전자가 어느 모듈에 속한다고 적혀 있는지이고, 그 소속은 공식 표에서 이미 생물학적으로 맞지 않는다.

예: skyblue(Heat stress response)는 파이프라인과 Table 5B가 같은 20유전자이다. 지시서에 적힌 Higd1a, Mmp2, Acads, …, Zfp609가 그 목록이고 heat shock 유전자가 아니다.

## 7. 단계 4 결과

수행하지 않았다. 공식 Table 5B가 집중도 검정에서 무작위 수준이므로 결정 규칙에 따라 수정 전에 중단했다.

## 8. 단계 5 결과

수행하지 않았다. 수정 diff, 수정 후 대조, 공개 계수 재현 표, composite 대조는 없다.

## 9. 단계 6 결과

수행하지 않았다. 재학습하지 않았다.

## 10. 단계 7 결과

수행하지 않았다. 영향 범위 표와 계산 시간 추정은 없다.

## 11. 미해결, 불확실, 예상 밖 사항

- 28개 모듈의 유전자 membership 표는 공개 보충자료에 없다. Table 6B는 이름, 주석, 크기만 있다. 따라서 “원 membership은 정상인데 저장소만 뒤섞였다”를 유전자 ID로 검증할 수 없다. 검증 가능한 원 유전자 목록은 Table 5B이고, 그 목록이 무작위 수준이다.
- 파이프라인이 Table 5B 비영 집합에 더하는 73유전자가 어떻게 붙었는지는 저장소에 생성 스크립트가 없어 확인하지 못했다. 모듈 크기는 Table 6에 가깝지만 brown, darkgreen, darkorange는 각각 9, 2, 1개 부족하다.
- 지시서의 재현 실패(공개 계수 in-sample r 중앙값 0.04, composite r = 0.998)는 이번 단계에서 다시 계산하지 않았다. 계수 값이 S3a와 공식 표에서 일치하므로, 그 실패가 이 저장소의 계수 붙여넣기 오류 때문이라는 가설은 이번 대조와 맞지 않는다. 추정: 공식 표의 유전자–모듈 소속 자체가 공동발현 모듈과 맞지 않거나, 적용 규약(스케일, 결측, 절편)이 다르다. 어느 쪽인지는 단계 5를 하지 않아 구분하지 못했다.
- 본문 “30–630 유전자”와 Table 6B의 12–357은 서로 다르다. 어느 쪽이 필터 전 크기인지는 유전자 목록이 없어 확인하지 못했다.
- 집중도 p값이 지시서 기준선과 소수 둘째 자리에서 다르다. top share는 같다.

## 12. 첨부 파일 목록

`results/audit_module_mapping/` 아래:

- `pipeline_modules.csv`
- `original_modules.csv` (Table 5B chronological, 비영 계수 유전자, 23모듈)
- `pipeline_modules_coherence.csv`, `s3a_mortality_modules_coherence.csv`, `original_modules_coherence.csv`
- `jaccard_matrix.csv`, `jaccard_summary.csv`
- `s3a_vs_pipeline.csv`, `s3a_mortality_modules.csv`
- `published_5b_nonzero_coefficients.csv`
- `original_module_dictionary_28.csv`
- `original_table4b_performance.csv`
- `module_size_reconciliation.csv`
- `training_matrix_overlap.csv`
- `coefficient_spotcheck_50.csv`, `coefficient_spotcheck_summary.csv`
- `source/41586_2026_10542_MOESM6_ESM.xlsx`, `MOESM7_ESM.xlsx`, `MOESM8_ESM.xlsx`

만들지 않은 파일: `published_reproduction.csv`, `retrain_cv.csv`, `impact_inventory.csv`, `fix.patch`. 단계 5–7을 하지 않았기 때문이다.

스크립트: `scripts/audit/module_coherence_audit.py`, `scripts/audit/build_steps_1_to_3.py`.
