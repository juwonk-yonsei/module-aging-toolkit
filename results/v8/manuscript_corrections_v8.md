# 원고 수정 목록

#7 목록(`results/v7/manuscript_corrections_v7.md`)에 이번 계산을 더했다. 기존 결과 파일은 덮어쓰지 않았다.

## 대조군 수

GSE136831 통계 격자의 대조군은 27명이다. control donor 244C는 통계 세포유형마다 세포가 50개 미만이라 격자 밖에 있다. 원고에 이 격자의 대조군을 28명으로 적은 곳이 있으면 27명으로 고친다.

## N과 MI

원고 2.8절의 정의는 성분 24개(모듈 23개와 모듈 밖 1개)에 대해 S = Σ|Δ_k|, N = |Σ Δ_k|, MI = (S − N) / S 이다.

r5 전역 표 `results/v3/human/r5/masking_global_r5.csv`의 N은 모듈 밖 성분의 절댓값 0.035669이고, MI는 0.807이다. 그 계산은 `scripts/v3/r5/r5_masking.py`에서 23개 모듈만 절댓값으로 더하고 마지막 성분의 절댓값을 N으로 두었다.

같은 점수와 같은 donor 순열로 정의를 다시 적용한 값(`results/v8/masking/masking_fixed_v8.csv`)은 다음과 같다.

- S = 0.220563 (p = 0.0010, 부트스트랩 구간 0.156–0.418)
- N = 0.023309 (p = 0.658, 구간 0.0028–0.178)
- MI = 0.894321 (p = 0.141, 구간 0.306–0.989)
- composite 계수 = −0.023309 (구간 −0.159–0.155)

v0.2.0 전역 표 `results/revision/rev07_masking_global.csv`의 N은 0.023309로, 수정된 N과 같다. 세포유형별 점추정은 r5 세포유형 표가 이미 |composite|를 N으로 두고 있어 값이 같다(`results/v8/masking/masking_celltype_fixed_v8.csv`).

## 원고에 넣을 수치의 출처

`results/v8/headline_v8.csv`의 각 행에 값, 단위, 출처 파일이 있다. 원고와 그림에 쓰는 숫자는 그 출처 파일에서 읽는다.

## v0.2.0 원고와 저장 결과가 다른 항목

- "Five effects passed both nulls." `results/revision/rev06_hits_competitive.csv`에서 두 경쟁 귀무의 FDR이 모두 0.10 미만인 효과는 10개이다. 그 문장은 이 파일에서 복원되지 않는다.
- S = 0.226, p = 0.0005. 저장된 masking 표에 S = 0.226은 없다. v0.2.0 전역 S는 0.216570이고 p_S = 0.0010이다(`results/revision/rev07_masking_global.csv`). r5가 저장한 S 0.184894는 23개 모듈만 더한 값이다. 24개 정의로 다시 계산한 S는 0.220563이다.
