"""Expressions that recompute each computed number of the manuscript from the result tables.

EXPR maps a claim id (claims_r2.csv) to an expression over the helpers of verify_claims_r2.py.
NOTE carries a short explanation, and for numbers that the released tables do not contain, the
reason (prefix 'NOT IN TABLES:'); such rows have no expression.
"""
EXPR, NOTE, KIND = {}, {}, {}   # KIND: final kind where it differs from the extraction rule

def E(cid, expr, note=''):
    EXPR[cid] = expr
    if note:
        NOTE[cid] = note
    else:
        NOTE.pop(cid, None)          # a later mapping replaces the note of an earlier one

def X(cid, note):
    NOTE[cid] = 'NOT IN TABLES: ' + note

AD, GRID = 'GSE136831 (Adams)', 'GSE136831 statistical grid'

# ----------------------------------------------------------------------------- Abstract (M009)
E('M0001', "v('S16_screen','k',gate='program',item='Q1')")
E('M0002', "v('S16_screen','n',gate='program',item='Q1')")
E('M0008', "100*v('S16_screen','rate',gate='clock',item='Q2')")
E('M0009', "100*v('S16_screen','rate',gate='program',item='Q2')")
E('M0010', "100*v('S16_input','set_failure_rate',tag='corrupted',f=0.1)")
E('M0013', f"v('S01_cohorts','n',cohort='{AD}',group='IPF')")
E('M0014', f"v('S01_cohorts','n',cohort='{GRID}',group='Control')")
E('M0015', f"v('S01_cohorts','n',cohort='{GRID}',group='COPD')")
E('M0016', "v('S09_global','beta',module='Composite')")
E('M0018', "v('S09_global','ci_lo',module='Composite')")
E('M0019', "v('S09_global','ci_hi',module='Composite')")
E('M0020', "v('S11_global','p_S')")
E('M0021', "v('S09_funnel','amended_r5',quantity='tests')")
E('M0022', "v('S09_funnel','amended_r5',quantity='hits_fdr10')")
E('M0023', "v('S09_funnel','amended_r5',quantity='ipf_specific')")
E('M0024', "100*v('S18_recovery_ci','recovery',unit='whole')")
E('M0025', "100*mn('S18_recovery_ci','recovery',query='unit.str.startswith(\"celltype\")')")
E('M0026', "100*mx('S18_recovery_ci','recovery',query='unit.str.startswith(\"celltype\")')")

# ----------------------------------------------------------------------------- Methods
E('M0067', "secs(v('S21_plans','first_result_mtime',plan='results/v3/human/analysis_plan_v3.md'),"
           " v('S21_plans','write_time_utc',plan='results/v3/human/analysis_plan_v3.md'))",
  'seconds from the first result file to the writing time of the v3 plan')
E('M0085', "col('S03_modules','n_genes_serialized').sum()")
E('M0087', "v('S02_detected_all','composite_genes_serialized')")
# cohorts (Section 2.5)
for cid, grp, field in (('M0124', 'IPF', 'n'), ('M0129', 'IPF', 'n_male'), ('M0130', 'Control', 'n'),
                        ('M0135', 'Control', 'n_male'), ('M0136', 'COPD', 'n'), ('M0139', 'COPD', 'n_male')):
    E(cid, f"v('S01_cohorts','{field}',cohort='{AD}',group='{grp}')")
for cid, grp, i in (('M0125', 'IPF', 0), ('M0126', 'IPF', 1), ('M0127', 'IPF', 2), ('M0128', 'IPF', 3),
                    ('M0131', 'Control', 0), ('M0132', 'Control', 1), ('M0133', 'Control', 2), ('M0134', 'Control', 3),
                    ('M0137', 'COPD', 0), ('M0138', 'COPD', 1)):
    E(cid, f"num(v('S01_cohorts','age',cohort='{AD}',group='{grp}'),{i})")
E('M0143', "n('S02_celltypes',atlas='GSE136831',in_statistics_set=True)")
E('M0145', "v('S09_funnel','planned_v3',quantity='tests')")
E('M0148', f"v('S01_cohorts','n',cohort='{AD}',group='IPF')")
E('M0149', f"v('S01_cohorts','n',cohort='{GRID}',group='Control')")
E('M0150', f"num(v('S01_cohorts','age',cohort='{GRID}',group='Control'),0)")
E('M0151', f"num(v('S01_cohorts','age',cohort='{GRID}',group='Control'),1)")
E('M0152', f"v('S01_cohorts','n_male',cohort='{GRID}',group='Control')")
E('M0153', f"v('S01_cohorts','n',cohort='{GRID}',group='COPD')")
E('M0154', f"num(v('S01_cohorts','age',cohort='{GRID}',group='COPD'),0)")
E('M0155', f"num(v('S01_cohorts','age',cohort='{GRID}',group='COPD'),1)")
E('M0156', f"v('S01_cohorts','n_male',cohort='{GRID}',group='COPD')")
E('M0161', "v('S01_cohorts','n',cohort='GSE135893 (Habermann)',group='IPF')")
E('M0162', "v('S01_cohorts','n',cohort='GSE135893 (Habermann)',group='Control')")
E('M0164', "v('S01_cohorts','n',cohort='GSE134692 (bulk)',group='IPF')")
E('M0165', "v('S01_cohorts','n',cohort='GSE134692 (bulk)',group='Control')")
# detection (Section 2.6)
E('M0170', "v('S02_detected_all','module_genes_detected')")
E('M0171', "v('S02_detected_all','module_genes_serialized')")
E('M0172', "100*v('S02_detected_all','module_genes_detected')/v('S02_detected_all','module_genes_serialized')")
E('M0173', "v('S02_detected_all','composite_genes_detected')")
E('M0174', "v('S02_detected_all','composite_genes_serialized')")
E('M0175', "100*v('S02_detected_all','composite_genes_detected')/v('S02_detected_all','composite_genes_serialized')")
E('M0176', "n('S02_detected')")
E('M0179', "v('S09_funnel','amended_r5',quantity='tests')")
E('M0180', "v('S12_sensitivity','beta_main',analysis='main',celltype='cMonocyte',module='black')", 'constant black-module score in classical monocytes (a hit before the estimability rule)')
E('M0181', "v('S09_funnel','planned_v3',quantity='tests')")

# ----------------------------------------------------------------------------- Methods, continued
E('M0188', "v('S09_funnel','amended_r5',quantity='tests')")
E('M0204', "v('S09_funnel','amended_r5',quantity='tests')")
E('M0210', "v('S02_detected_all','composite_genes_serialized')")
E('M0217', f"num(v('S01_cohorts','age',cohort='{AD}',group='IPF'),2)", 'youngest IPF donor')
E('M0232', "v('S02_detected_all','module_genes_serialized')")
E('M0237', "num(v('S01_cohorts','note',cohort='HLCA analyzed'),2)", 'eligible donors, from the note of the HLCA row')
E('M0239', "v('S01_cohorts','n',cohort='HLCA analyzed')")
E('M0240', "v('S01_cohorts','n',cohort='GTEx v8 lung')")
E('M0243', "v('S18_positive','n_young_mice')")
E('M0250', "v('S18_resolution','expected_shift',species='mouse')/21", 'slope = expected 21-month shift / 21')
E('M0251', "v('S18_resolution','expected_shift',species='mouse')")

# ----------------------------------------------------------------------------- Results 3.1 (input integrity)
SPR = "source='supp_table',convention='as_released_na0'"
MODS = "query=\"module not in ['All module genes','composite_5A'] and retained == True\""
E('M0282', f"med('S05_insample','r_insample',{MODS},{SPR},outcome='Mortality')")
E('M0284', f"v('S05_insample','r_insample',{SPR},outcome='Mortality',module='All module genes')")
E('M0285', f"v('S05_insample','r_insample',{SPR},outcome='Chronological',module='All module genes')")
E('M0288', "100*v('S05_sign','simple',source='5B',outcome='Mortality',module='All module genes')")
E('M0289', "v('S05_sign','p_simple',source='5B',outcome='Mortality',module='All module genes')")
E('M0290', "n('S05_sign',source='5B',outcome='Mortality',query=\"module != 'All module genes' and p_simple >= 0.05\")",
  'modules whose sign agreement is consistent with the shuffled null (p >= 0.05)')
E('M0292', "med('S05_coexpr_table','mean_pairwise_r')")
E('M0294', "mn('S05_coherence','p',source='table5b') > 0.45", 'bound')
E('M0296', f"v('S05_insample','r_insample',{SPR},outcome='Mortality',module='composite_5A')")
SER = "query=\"module != 'allmodulegenes'\""
E('M0298', f"med('S05_serialized','r_native',{SER},outcome='Mortality')")
E('M0299', f"med('S05_serialized','r_native',{SER})")
E('M0301', "n('S05_serialized',query=\"module != 'allmodulegenes' and r_native >= r_reported\")")
E('M0304', "v('S05_serialized','r_native',outcome='Mortality',module='allmodulegenes')")
E('M0305', "v('S05_serialized','r_native',outcome='Chronological',module='allmodulegenes')")
E('M0308', "med('S17_input','rbar',tag='official')")
E('M0309', "mn('S17_input','z',tag='official')")
FAM = {'ribo': 'ribosomal (Rpl/Rps)', 'ox': 'OXPHOS (Nduf/Cox/Uqcr/Sdh/Atp5)', 'ifn': 'interferon-stimulated',
       'cc': 'cell cycle', 'col': 'collagen (Col)', 'hs': 'heat shock (Hsp/Dnaj)'}
for cid, fam, field, scale in (('M0310', 'ribo', 'top_share', 100), ('M0311', 'ribo', 'n_genes', 1),
                               ('M0312', 'ox', 'top_share', 100), ('M0313', 'ox', 'n_genes', 1),
                               ('M0314', 'ifn', 'top_share', 100), ('M0315', 'ifn', 'n_genes', 1),
                               ('M0316', 'cc', 'top_share', 100), ('M0317', 'cc', 'n_genes', 1),
                               ('M0318', 'col', 'n_genes', 1), ('M0320', 'hs', 'top_share', 100),
                               ('M0321', 'hs', 'n_genes', 1), ('M0322', 'hs', 'p', 1)):
    E(cid, f"{scale}*v('S05_coherence','{field}',source='official',family='{FAM[fam]}')")
E('M0319', f"mx('S05_coherence','p',source='official',query=\"family != '{FAM['hs']}'\")", 'largest p of the five families')

# ----------------------------------------------------------------------------- 3.1 corruption and scrambled memberships
E('M0326', "v('S16_input','set_failure_k',tag='corrupted',f=0.1)")
E('M0328', "100*v('S16_input','set_failure_rate',tag='corrupted',f=0.1)")
E('M0332', "mn('S16_input','set_failure_k',tag='corrupted',query='f >= 0.25')", 'all 20 failed at f = 0.25, 0.5 and 1')
CORR = "tag='corrupted'"
E('M0339', f"mean('S17_corrupt','n_program_pass',{CORR},path='V1')")
E('M0340', f"mean('S17_corrupt','n_clock_pass',{CORR},path='V1')")
E('M0341', f"mean('S17_corrupt','n_program_pass',{CORR},path='V2')")
E('M0342', f"mean('S17_corrupt','n_robust',{CORR},path='V1')")
E('M0343', f"mn('S17_corrupt','n_robust',{CORR},path='V1')")
E('M0344', f"mx('S17_corrupt','n_robust',{CORR},path='V1')")
E('M0345', f"mean('S17_corrupt','n_robust',{CORR},path='V2')")
E('M0346', f"mn('S17_corrupt','n_robust',{CORR},path='V2')")
E('M0347', f"mx('S17_corrupt','n_robust',{CORR},path='V2')")
E('M0348', "v('S17_corrupt','n_robust',tag='official')")
E('M0351', "v('S17_scrambled','classical_median')")
E('M0354', "v('S17_scrambled','classical_p2_5')")
E('M0355', "v('S17_scrambled','classical_p97_5')")
E('M0356', "med('S17_scrambled_raw','hits_Macrophage')")
E('M0357', "med('S17_scrambled_raw','hits_Macrophage_Alveolar')")
E('M0358', "v('S09_funnel','amended_r5',quantity='hits_fdr10')")
E('M0359', "v('S17_scrambled','percentile_amended_10')")
E('M0360', "v('S17_scrambled','original_submission_hits')")
E('M0361', "v('S17_scrambled','percentile_original_submission_30')")

# ----------------------------------------------------------------------------- 3.3 rodent domain (S04, S06)
E('M0418', "med('S04_cv','r_oof_all',outcome='Mortality')")
E('M0419', "med('S04_cv','r_oof_all',outcome='Chronological')")
E('M0420', "med('S04_cv','r_oof_all_encv',outcome='Mortality')")
E('M0421', "med('S04_cv','r_oof_all_encv',outcome='Chronological')")
E('M0422', "med('S04_cv','r_oof_mouse',outcome='Mortality')")
E('M0423', "med('S04_cv','r_oof_rat',outcome='Mortality')")
E('M0424', "med('S04_cv','r_mouse_to_rat',outcome='Mortality')")
E('M0425', "med('S04_cv','r_rat_to_mouse',outcome='Mortality')")
E('M0426', "med('S04_cv','r_within_rat_lodo',outcome='Mortality')")
E('M0427', "med('S04_cv','r_random_mean',outcome='Mortality')")
E('M0428', "med('S04_cv','r_random_mean',outcome='Chronological')")
E('M0429', "n('S04_cv',exceeds_all_random=True)")
E('M0431', "v('S04_original','median_heldout_r',outcome='Mortality')")
E('M0432', "med('S04_cv','r_oof_all_encv',outcome='Mortality')")
E('M0433', "float((rows('S04_cv')['r_oof_all'] - rows('S04_cv')['r_oof_equal_sign']).median())",
  'median of held-out r minus the equal-weight signed score')
E('M0439', "v('S06_meff','median_abs_r',outcome='Mortality',method='spearman')")
E('M0440', "v('S06_meff','meff_liji',outcome='Mortality',method='spearman')")
E('M0442', "v('S06_meff','median_abs_r',query='dataset.notna()')")
E('M0443', "v('S06_meff','meff_liji',query='dataset.notna()')")

# ----------------------------------------------------------------------------- 3.2 spike-ins (S16)
E('M0362', "v('S16_screen','k',gate='program',item='Q1')")
E('M0363', "v('S16_screen','n',gate='program',item='Q1')")
E('M0365', "v('S16_program','k',item='P1')")
E('M0366', "v('S16_program','n',item='P1')")
E('M0376', "v('S16_screen','n',gate='program',item='Q3')")
E('M0377', "100*v('S16_screen','rate',gate='program',item='Q3')")
E('M0378', "100*v('S16_screen','rate',gate='clock',item='Q3')")
E('M0379', "100*v('S16_screen','rate',gate='program',item='Q3_planted')")
E('M0380', "v('S16_program','k',item='P3')")
E('M0381', "v('S16_program','n',item='P3')")
E('M0383', "100*v('S16_program','rate',item='P2')")
E('M0385', "100*v('S16_bias','rate_set_specific',table='by_class',scenario='S2a',d='all',class_or_denominator='false_hit')")
E('M0386', "v('S16_bias','n_hits',table='by_class',scenario='S2a',d='all',class_or_denominator='false_hit')")
E('M0387', "v('S16_screen','n',gate='program',item='Q2')")
E('M0390', "100*v('S16_screen','rate',gate='program',item='Q2')")
E('M0391', "100*v('S16_screen','rate',gate='clock',item='Q2')")
E('M0401', "n('S02_detected')")

# ----------------------------------------------------------------------------- 3.3 Klotho (S07)
KC = 'S07_klotho_composite'
E('M0444', f"v('{KC}','delta',tissue='Kidney')")
E('M0446', f"v('{KC}','ci_lo',tissue='Kidney')")
E('M0447', f"v('{KC}','ci_hi',tissue='Kidney')")
E('M0448', f"v('{KC}','p_exact',tissue='Kidney')")
E('M0449', f"v('{KC}','delta',tissue='Skeletal muscle')")
E('M0450', f"v('{KC}','ci_lo',tissue='Skeletal muscle')")
E('M0451', f"v('{KC}','ci_hi',tissue='Skeletal muscle')")
E('M0452', f"v('{KC}','p_exact',tissue='Skeletal muscle')")
E('M0453', "float(col('S07_klotho','delta',tissue='Kidney').abs().median())")
E('M0454', "float(col('S07_klotho','delta',tissue='Skeletal muscle').abs().median())")
E('M0457', "v('S07_klotho_count','n_positive',tissue='Kidney')")
E('M0459', "v('S07_klotho_count','n_positive',tissue='Skeletal muscle')")
E('M0460', "v('S07_klotho_count','n_modules',tissue='Skeletal muscle')")
E('M0462', "v('S07_klotho_spearman','n_positive_in_both')")
E('M0463', "v('S07_klotho_spearman','pairing_null_median')")
E('M0464', "v('S07_klotho_spearman','pairing_one_sided_p')")
E('M0465', "v('S07_klotho_spearman','rho')")
E('M0466', "v('S07_klotho_spearman','two_sided_p')")
E('M0467', "v('S07_klotho_count','null_count_median',tissue='Kidney')")
E('M0468', "v('S07_klotho_count','p_count_exact',tissue='Kidney')")
E('M0469', "v('S07_klotho_count','p_count_exact',tissue='Skeletal muscle')")
E('M0470', "n('S07_klotho',tissue='Kidney',query='fdr_perm < 0.10')")
E('M0471', "n('S07_klotho',tissue='Skeletal muscle',query='fdr_perm < 0.10')")

# ----------------------------------------------------------------------------- 3.3 human age (S08)
HC = "cohort='donors_all',outcome='Mortality'"
E('M0474', "v('S01_cohorts','n',cohort='HLCA analyzed')")
E('M0475', "num(v('S01_cohorts','age',cohort='HLCA analyzed'),0)")
E('M0476', "num(v('S01_cohorts','age',cohort='HLCA analyzed'),1)")
E('M0477', f"v('S08_hlca_compare','n_donors_min',{HC})")
E('M0478', f"v('S08_hlca_compare','n_donors_max',{HC})")
E('M0480', f"v('S08_hlca_compare','median_rho',{HC})")
E('M0481', "mn('S08_hlca77','min_abs_rho_80power',outcome='Mortality')")
E('M0482', "mx('S08_hlca77','min_abs_rho_80power',outcome='Mortality')")
E('M0485', "v('S01_cohorts','n',cohort='HLCA 42-donor cap')")
E('M0486', "v('S08_hlca_compare','fdr10_raw',cohort='donors_42',outcome='Mortality')")
E('M0487', "v('S08_hlca_compare','fdr10_raw_negative',cohort='donors_42',outcome='Mortality')")
AC = "outcome='Mortality',query=\"module != 'Composite'\""
E('M0488', f"mn('S08_adams_controls','n',{AC})")
E('M0489', f"mx('S08_adams_controls','n',{AC})")
E('M0491', f"med('S08_adams_controls','rho',{AC})")
E('M0492', "v('S08_gtex','n',outcome='Mortality',module='Composite')")
E('M0493', "v('S08_gtex','adj_spearman',outcome='Mortality',module='Composite')")
E('M0494', "v('S08_gtex','adj_p',outcome='Mortality',module='Composite')")
E('M0495', "v('S08_gtex','adj_spearman',outcome='Chronological',module='Composite')")
E('M0498', "v('S08_gtex','adj_spearman',outcome='Mortality',module='blue')")
E('M0499', "v('S08_gtex','adj_spearman',outcome='Mortality',module='tan')")
E('M0500', "v('S08_gtex','adj_spearman',outcome='Mortality',module='salmon')")

# ----------------------------------------------------------------------------- 3.3 aging mouse lung (S18)
RC = 'S18_recovery_ci'
CT = "query='unit.str.startswith(\"celltype\")'"
E('M0502', f"v('{RC}','beta',unit='whole')")
E('M0506', f"v('{RC}','ci_lo',unit='whole')")
E('M0507', f"v('{RC}','ci_hi',unit='whole')")
E('M0508', f"v('{RC}','p',unit='whole')")
E('M0509', f"100*v('{RC}','recovery',unit='whole')")
E('M0510', f"v('{RC}','expected_shift',unit='whole')")
E('M0511', f"100*mn('{RC}','recovery',{CT})")
E('M0512', f"100*mx('{RC}','recovery',{CT})")
for cid, lab, field in (('M0514', 'AT2', 'beta'), ('M0515', 'AT2', 'p'), ('M0516', 'Alv. mac', 'beta'),
                        ('M0517', 'Alv. mac', 'p'), ('M0518', 'Ciliated', 'beta'), ('M0519', 'Ciliated', 'p')):
    E(cid, f"v('{RC}','{field}',label='{lab}')")
E('M0520', f"v('{RC}','beta',unit='whole_detection_limited')")
E('M0521', f"100*v('{RC}','recovery',unit='composition_only')")
E('M0522', f"v('{RC}','beta',unit='composition_only')")
E('M0524', f"v('{RC}','ci_lo',unit='composition_only')")
E('M0525', f"v('{RC}','ci_hi',unit='composition_only')")
E('M0526', f"v('{RC}','p',unit='composition_only')")
E('M0527', f"100*v('{RC}','recovery',unit='expression_only')")
E('M0528', f"v('{RC}','beta',unit='expression_only')")
E('M0529', f"v('{RC}','ci_lo',unit='expression_only')")
E('M0530', f"v('{RC}','ci_hi',unit='expression_only')")
E('M0531', f"v('{RC}','p',unit='expression_only')")
E('M0532', "v('S18_resolution','ratio_whole_over_median',species='mouse')")
E('M0533', "100*v('S18_units','recovery',species='mouse',unit='whole',score='composite_Chronological')")
E('M0534', "v('S18_units','p',species='mouse',unit='whole',score='composite_Chronological')")
CC = "species='mouse',outcome='Mortality'"
for cid, ct, field, scale in (('M0536', 'Ciliated_cells', 'contrib_ratio', 100), ('M0537', 'Ciliated_cells', 'beta_c', 1),
                              ('M0539', 'Ciliated_cells', 'ci_lo', 1), ('M0540', 'Ciliated_cells', 'ci_hi', 1),
                              ('M0541', 'Alveolar_macrophage', 'contrib_ratio', 100), ('M0542', 'Alveolar_macrophage', 'beta_c', 1),
                              ('M0543', 'Alveolar_macrophage', 'ci_lo', 1), ('M0544', 'Alveolar_macrophage', 'ci_hi', 1),
                              ('M0546', 'Type_2_pneumocytes', 'contrib_ratio', 100),
                              ('M0550', 'Ciliated_cells', 'standalone_beta', 1), ('M0551', 'Ciliated_cells', 'beta_c', 1),
                              ('M0552', 'Alveolar_macrophage', 'standalone_beta', 1), ('M0553', 'Alveolar_macrophage', 'beta_c', 1),
                              ('M0555', 'Type_2_pneumocytes', 'standalone_beta', 1), ('M0556', 'Type_2_pneumocytes', 'beta_c', 1)):
    E(cid, f"{scale}*v('S18_contrib_cells','{field}',{CC},celltype='{ct}')")
E('M0547', "100*v('S18_contrib','additivity_ratio',species='mouse',outcome='Mortality')")
E('M0557', f"100*mx('{RC}','recovery',{CT})", 'largest cell-type recovery')
E('M0558', "v('S18_resolution','whole_beta',species='human')")
E('M0559', "v('S18_resolution','whole_p',species='human')")
E('M0560', "v('S18_resolution','ratio_whole_over_median',species='human')")
E('M0561', "v('S18_resolution','median_celltype_recovery',species='human')")
E('M0564', "v('S18_contrib_cells','p_c',species='human',outcome='Mortality',celltype='expression_only')")
for cid, ct in (('M0565', 'Macrophage'), ('M0566', 'T'), ('M0567', 'Stromal')):
    E(cid, f"100*v('S18_contrib_cells','contrib_ratio',species='human',outcome='Mortality',celltype='{ct}')")
E('M0569', "v('S16_criteria','k2',criterion='C')")
E('M0570', "v('S16_criteria','n2',criterion='C')")
E('M0571', "100*v('S16_criteria','rate2',criterion='C')")
E('M0572', "v('S18_positive','positive_p')")
E('M0573', "v('S18_tms','beta',unit='whole',score='composite_Mortality')")
E('M0574', "v('S18_tms','p',unit='whole',score='composite_Mortality')")

# ----------------------------------------------------------------------------- 3.4 composite decomposition (S09_global, S11)
E('M0575', "v('S09_global','beta',module='Composite')")
E('M0577', "v('S09_global','ci_lo',module='Composite')")
E('M0578', "v('S09_global','ci_hi',module='Composite')")
E('M0579', "v('S11_global','composite_ci_lo',level='global')")
E('M0580', "v('S11_global','composite_ci_hi',level='global')")
E('M0581', "v('S09_global','p_cr1_t',module='Composite')")
E('M0584', "n('S11_components',query='fdr < 0.10',level='component')")
for cid, field in (('M0586', 'S'), ('M0588', 'S_ci_lo'), ('M0589', 'S_ci_hi'), ('M0591', 'p_S'), ('M0592', 'N'),
                   ('M0593', 'p_N'), ('M0594', 'MI'), ('M0596', 'MI_ci_lo'), ('M0597', 'MI_ci_hi'), ('M0599', 'p_MI')):
    E(cid, f"v('S11_global','{field}',level='global')")
E('M0590', "v('S11_null','null_median',statistic='S')")
E('M0598', "v('S11_null','null_median',statistic='MI')")
E('M0601', "n('S11_celltype')")
E('M0604', "v('S09_global','beta',module='red')")
E('M0605', "100*v('S09_red_library','ribosomal_coef_fraction',scope='global')")
E('M0606', "v('S09_red_library','rho',scope='Macrophage',contrast='spearman score vs lib_size')")
E('M0607', "v('S09_red_library','p',scope='Macrophage',contrast='spearman score vs lib_size')")
E('M0608', "v('S09_red_library','rho',scope='cMonocyte',contrast='spearman score vs lib_size')")
E('M0609', "v('S09_red_library','p',scope='cMonocyte',contrast='spearman score vs lib_size')")
E('M0610', "v('S09_red_library','rho',scope='GSE134692')")
E('M0611', "v('S09_red_library','rho',scope='GSE135893:Macrophages')")
E('M0612', "v('S09_red_library','rho',scope='GSE135893:Monocytes')")

# ----------------------------------------------------------------------------- 3.5 IPF grid (S09_ipf, S09_funnel, S10_screen)
HIT = "query='fdr_IPF < 0.10'"
E('M0613', "n('S09_ipf')", 'estimable tests')
E('M0614', f"n('S09_ipf',{HIT})")
E('M0615', f"n('S09_ipf',{HIT},robust=True)")
E('M0616', f"n('S09_ipf',{HIT},robust=True,IPF_specific=True)")
E('M0618', "mx('S19_mde','abs_standardized',celltype='Macrophage',module='grey60')")
E('M0619', "-mx('S19_mde','abs_standardized',celltype='Macrophage_Alveolar',module='green')", 'sign of beta_IPF (negative)')
E('M0620', "v('S09_ipf','beta_IPF',celltype='cMonocyte',module='red')")
E('M0621', "v('S09_ipf','beta_COPD',celltype='cMonocyte',module='red')")
E('M0622', "n('S10_screen')")
E('M0623', "min(mn('S10_screen','fdr_screen_random'), mn('S10_screen','fdr_screen_perm'))")
E('M0624', "n('S10_screen')")
E('M0625', "(v('S10_screen','fdr_screen_random',celltype='Macrophage',module='grey60'), "
           "v('S10_screen','fdr_screen_random',celltype='Macrophage_Alveolar',module='green'))")
E('M0626', "v('S10_screen','p_screen_random',celltype='Macrophage',module='grey60')")
E('M0627', "v('S10_screen','p_screen_random',celltype='Macrophage_Alveolar',module='green')")
E('M0631', f"sum(v('S03_modules','heldout_r_Mortality',colour=m) < 0.30 for m in rows('S09_ipf',{HIT}).module) == 9 "
           "and all(v('S03_modules','heldout_r_Mortality',colour=m) < 0.30 for m in ('grey60','green'))",
  'nine of the ten hits, including both IPF-specific ones, in modules with dataset-held-out r < 0.30')
E('M0632', "v('S09_funnel','planned_v3',quantity='tests')")
E('M0634', "v('S09_funnel','planned_v3',quantity='hits_fdr10')")
E('M0635', "v('S09_funnel','planned_v3',quantity='robust')")
E('M0636', "v('S09_funnel','planned_v3',quantity='ipf_specific')")
E('M0637', "v('S12_sensitivity','fdr',analysis='main',celltype='NK',module='greenyellow')", 'S12 main analysis holds all 253 planned tests with FDR over 253')
E('M0638', "v('S09_ipf','fdr_IPF',celltype='NK',module='greenyellow')")

# ----------------------------------------------------------------------------- 3.5 sensitivity (S12) and baselines (S10)
def S12H(analysis):
    """S12 rows of one analysis for the ten hits of the 249-test grid."""
    return (f"rows('S12_sensitivity',analysis='{analysis}')"
            ".merge(rows('S09_ipf',query='fdr_IPF < 0.10')[['celltype','module']])")
AGE = S12H('age_overlap') + ".dropna(subset=['ratio'])"
E('M0639', f"100*{AGE}.ratio.min()")
E('M0640', f"100*{AGE}.ratio.max()")
E('M0641', f"len({AGE}) == 6 and {AGE}.p.max() < 0.05", 'six testable hits, all p < 0.05')
TECH = ('cpm', 'min_cells_30', 'min_cells_100', 'depth_group_matched', 'cells_group_matched')
E('M0642', "min(" + ", ".join(f"float({S12H(a)}.ratio.median())" for a in TECH) + ")")
E('M0643', "max(" + ", ".join(f"float({S12H(a)}.ratio.median())" for a in TECH) + ")")
E('M0646', f"float({S12H('tech_covariates')}.ratio.median())")
for cid, ct, mod, an in (('M0647', 'Macrophage', 'turquoise', 'tech_covariates'),
                         ('M0648', 'ncMonocyte', 'green', 'cells_group_matched'),
                         ('M0649', 'ncMonocyte', 'green', 'cells_equalized'),
                         ('M0650', 'ncMonocyte', 'green', 'cells_loss_control'),
                         ('M0651', 'ncMonocyte', 'tan', 'depth_equalized'),
                         ('M0652', 'ncMonocyte', 'tan', 'depth_loss_control'),
                         ('M0653', 'Ciliated', 'lightyellow', 'min_cells_30')):
    E(cid, f"v('S12_sensitivity','ratio',celltype='{ct}',module='{mod}',analysis='{an}')")
E('M0658', "n('S10_baselines')")
E('M0660', "n('S10_baselines',query='fdr_unweighted_249 < 0.10')")
E('M0661', "n('S10_baselines',query='fdr_signed_249 < 0.10')")
E('M0662', "corr('S10_baselines','t_clock','t_geneset','spearman')")
E('M0663', "corr('S10_baselines','t_clock','t_geneset_signed','spearman')")

# ----------------------------------------------------------------------------- 3.5 alveolar-macrophage cell-cycle candidate (S19, S13_specific)
def DIR(measure, field):
    return f"v('S19_direction','{field}',celltype='Macrophage_Alveolar',module='green',measure='{measure}')"
def DIRCI(measure, i):
    return f"tci({DIR(measure, 'coefficient')}, {DIR(measure, 't')}, {DIR(measure, 'dof')})[{i}]"
E('M0665', DIR('official_clock', 'coefficient'))
E('M0666', DIR('unweighted_mean', 'correlation_with_clock'))
E('M0667', DIR('pc1', 'coefficient'))
E('M0669', DIRCI('pc1', 0), '95% CI from coefficient, t and residual df')
E('M0670', DIRCI('pc1', 1))
E('M0671', DIR('pc1', 'correlation_with_clock'))
E('M0672', DIR('unweighted_mean', 'coefficient'))
E('M0673', DIRCI('unweighted_mean', 0))
E('M0674', DIRCI('unweighted_mean', 1))
for cid, mk, field in (('M0675', 'MKI67', 'coefficient'), ('M0676', 'MKI67', 'ci_lo'), ('M0677', 'MKI67', 'ci_hi'),
                       ('M0678', 'TOP2A', 'coefficient'), ('M0679', 'TOP2A', 'ci_lo'), ('M0680', 'TOP2A', 'ci_hi')):
    E(cid, f"v('S19_markers','{field}',marker='{mk}')")
E('M0681', "(v('S19_markers','n',marker='MKI67'), v('S19_markers','n',marker='TOP2A'))")
E('M0683', "v('S13_specific','p_H',celltype='Macrophage_Alveolar',module='green')")
E('M0684', "v('S13_specific','p_H',celltype='Macrophage',module='grey60')")
E('M0685', "v('S13_replication','n_IPF_GSE135893',coarse='Macrophage')")
E('M0686', "v('S13_replication','n_Control_GSE135893',coarse='Macrophage')")
PW = ("power_t(mx('S19_mde','abs_standardized',celltype='{ct}'), v('S13_replication','n_IPF_GSE135893',coarse='Macrophage'), "
      "v('S13_replication','n_Control_GSE135893',coarse='Macrophage'))")
E('M0687', PW.format(ct='Macrophage'), 'power of a two-sided t test, df = n - 4, as in S19_power')
E('M0688', PW.format(ct='Macrophage_Alveolar'))

# ----------------------------------------------------------------------------- 3.6 replication (S13), checklist (S14, S15)
REP = (('M0690', 'NK', 'r_centred'), ('M0691', 'NK', 'p_two_sided'), ('M0693', 'T cell', 'r_centred'),
       ('M0694', 'T cell', 'p_two_sided'), ('M0695', 'Macrophage', 'r_centred'), ('M0696', 'Macrophage', 'p_two_sided'),
       ('M0697', 'Monocyte', 'r_centred'), ('M0698', 'Monocyte', 'p_two_sided'), ('M0699', 'Ciliated', 'r_centred'),
       ('M0700', 'Ciliated', 'p_two_sided'), ('M0701', 'B cell', 'r_centred'), ('M0702', 'B cell', 'p_two_sided'),
       ('M0704', 'AT2', 'r_centred'), ('M0705', 'AT2', 'p_two_sided'), ('M0706', 'NK', 'r_partial_uniform'),
       ('M0707', 'NK', 'nk_lodo_r_min'), ('M0708', 'NK', 'nk_lodo_r_max'), ('M0709', 'NK', 'splithalf_r'),
       ('M0710', 'NK', 'splithalf_tau'))
for cid, ct, field in REP:
    E(cid, f"v('S13_replication','{field}',coarse='{ct}')")
E('M0711', "n('S13_specific',amended_249_hit=True,same_sign=True)")
E('M0712', "n('S13_specific',amended_249_hit=True)")
E('M0713', "mn('S13_specific','p_H',amended_249_hit=True) >= 0.05", 'none of the ten hits reached p < 0.05 in GSE135893')
for cid, ct, field in (('M0716', 'Macrophage', 'r_splithalf'), ('M0717', 'Monocyte', 'r_splithalf'),
                       ('M0718', 'Macrophage', 'tau_splithalf'), ('M0719', 'Monocyte', 'tau_splithalf'),
                       ('M0720', 'Macrophage', 'n_IPF_H'), ('M0721', 'Macrophage', 'n_Ctrl_H'),
                       ('M0722', 'Monocyte', 'n_IPF_H'), ('M0723', 'Monocyte', 'n_Ctrl_H')):
    E(cid, f"v('S15_real','{field}',coarse='{ct}')")
SUB = "v('S14_substate','r',split='{sp}',celltype='{ct}',substate='{ss}')"
E('M0724', SUB.format(sp='marker sub-state', ct='Macrophage', ss='monoHi'))
E('M0725', SUB.format(sp='marker sub-state', ct='Macrophage', ss='resHi'))
E('M0726', "v('S14_substate_p','bonferroni_p',celltype='Macrophage')")
E('M0727', SUB.format(sp='marker sub-state', ct='Monocyte', ss='monoHi'))
E('M0728', SUB.format(sp='marker sub-state', ct='Monocyte', ss='resHi'))
E('M0729', "v('S14_substate_p','bonferroni_p',celltype='Monocyte')")
for cid, sp, ct in (('M0731', 'marker sub-state', 'Macrophage'), ('M0732', 'random split', 'Macrophage'),
                    ('M0733', 'marker sub-state', 'Monocyte'), ('M0734', 'random split', 'Monocyte')):
    E(cid, f"mean('S14_substate','r',split='{sp}',celltype='{ct}')", 'mean of the two sub-states')
for cid, coh, field in (('M0736', 'GSE136831', 'median_IPF'), ('M0737', 'GSE136831', 'median_Control'),
                        ('M0738', 'GSE135893', 'median_IPF'), ('M0739', 'GSE135893', 'median_Control')):
    E(cid, f"v('S13_composition','{field}',cohort='{coh}',cell_type='Macrophage')")
E('M0740', "v('S13_monocyte','mono_frac_IPF',celltype='Macrophage')")
E('M0741', "v('S13_monocyte','mono_frac_Control',celltype='Macrophage')")
E('M0742', "v('S13_monocyte','p_mannwhitney',celltype='Macrophage')")
E('M0744', "n('S14_interaction',celltype='Macrophage',query='fdr_interaction < 0.10')")
E('M0746', "n('S14_interaction',celltype='Monocyte',query='fdr_interaction < 0.10')")
for cid, pair, field in (('M0749', 'GSE136831–GSE135893', 'r'), ('M0750', 'GSE136831–GSE135893', 'p'),
                         ('M0751', 'GSE136831–GSE134692', 'r'), ('M0752', 'GSE136831–GSE134692', 'p'),
                         ('M0753', 'GSE135893–GSE134692', 'r'), ('M0754', 'GSE135893–GSE134692', 'p')):
    E(cid, f"v('S13_perm','{field}',pair='{pair}')")
for cid, coh, comp in (('M0755', 'GSE136831', 'composition_only'), ('M0756', 'GSE135893', 'composition_only'),
                       ('M0757', 'GSE136831', 'expression_only'), ('M0758', 'GSE135893', 'expression_only')):
    E(cid, f"v('S13_decomposition','r_with_observed',cohort='{coh}',component='{comp}')")
for cid, sc, field in (('M0760', 'shared effect', 'replicated'), ('M0761', 'sub-state composition', 'composition'),
                       ('M0762', 'confound in one cohort', 'between-cohort difference'),
                       ('M0763', 'no effect', 'inconclusive'), ('M0764', 'weak shared effect', 'replicated'),
                       ('M0765', 'weak shared effect', 'between-cohort difference'),
                       ('M0766', 'composition shift only', 'replicated')):
    E(cid, f"100*v('S15_confusion','{field}',scenario='{sc}')")

# ----------------------------------------------------------------------------- remaining Results and Discussion
E('M0441', "v('S06_meff','n_modules',outcome='Mortality',method='spearman')")
E('M0775', "all(1.0 <= 100*v('S16_screen','rate',gate=g,item='Q2') < 3.0 for g in ('program','clock'))",
  'about 2%: screen-definition false-hit rates 2.4% (set) and 1.4% (weight)')
E('M0779', "100*v('S18_contrib','max_contrib_ratio',species='mouse',outcome='Mortality')")
E('M0780', "v('S01_cohorts','n',cohort='HLCA analyzed')")
E('M0784', "n('S10_baselines')")
E('M0785', "corr('S10_baselines','t_clock','t_geneset','spearman')")
E('M0786', "corr('S10_baselines','t_clock','t_geneset_signed','spearman')")
E('M0789', "35 <= 100*v('S15_confusion','between-cohort difference',scenario='weak shared effect') <= 45",
  'about 40%: weak shared effects attributed to a between-cohort difference (39%)')

# ----------------------------------------------------------------------------- figure legends
NCT = "len(set(rows('S09_ipf').celltype))"
E('M0811', NCT)
E('M0829', "tuple(col('S05_coherence','p',source='official',query='family != \"heat shock (Hsp/Dnaj)\"'))",
  'five families at the smallest attainable p, 1/1001')
E('M0830', "v('S05_coherence','p',source='official',family='heat shock (Hsp/Dnaj)')")
E('M0831', "mn('S05_coherence','p',source='table5b') > 0.45")
E('M0841', "n('S09_ipf',query='fdr_IPF < 0.10')")
E('M0842', "v('S17_scrambled','original_submission_hits')")
E('M0844', "v('S16_criteria','k',criterion='P3')")
E('M0845', "v('S16_criteria','n',criterion='P3')")
E('M0849', "v('S16_screen','n',gate='program',item='Q2')")
E('M0852', NCT)
E('M0857', "v('S01_cohorts','n',cohort='HLCA analyzed')")
E('M0858', "len(set(T('S08_hlca77').celltype))")
E('M0860', "v('S01_cohorts','n',cohort='GTEx v8 lung')")
E('M0863', "v('S08_gtex','adj_spearman',outcome='Mortality',module='Composite')")
E('M0866', "v('S18_recovery_ci','expected_shift',unit='whole')")
E('M0868', "v('S01_cohorts','n',cohort='GSE124872 Angelidis',group='young')")
E('M0869', "v('S01_cohorts','n',cohort='GSE124872 Angelidis',group='old')")
E('M0871', "v('S18_contrib','additivity_ratio',species='mouse',outcome='Mortality')")
E('M0872', "v('S18_contrib','additivity_ratio',species='human',outcome='Mortality')")
E('M0879', "n('S09_ipf')")
E('M0881', NCT)
E('M0883', "v('S02_celltypes','donors_Control',atlas='GSE136831',celltype='Macrophage_Alveolar')")
E('M0884', "v('S02_celltypes','donors_IPF',atlas='GSE136831',celltype='Macrophage_Alveolar')")
E('M0885', "v('S02_celltypes','donors_COPD',atlas='GSE136831',celltype='Macrophage_Alveolar')")

# ----------------------------------------------------------------------------- Table 1 and Table 2
# Table 2 rows are read from S16_criteria, which holds the 14 criteria exactly as Table 2 lists them
# (k, n, rate; k2, n2, rate2 for the second population of A2 and the second part of C; p for C)
T2 = "v('S16_criteria','{f}',criterion='{c}')"
for cid, c, f in (('M0928', 'A1', 'k'), ('M0929', 'A1', 'n'),
                  ('M0931', 'A2', 'k'), ('M0932', 'A2', 'n'), ('M0933', 'A2', 'rate'),
                  ('M0934', 'A2', 'k2'), ('M0935', 'A2', 'n2'), ('M0936', 'A2', 'rate2'),
                  ('M0938', 'A3', 'k'), ('M0939', 'A3', 'n'),
                  ('M0942', 'B', 'k'), ('M0943', 'B', 'n'), ('M0944', 'B', 'rate'),
                  ('M0948', 'C', 'k'), ('M0949', 'C', 'n'), ('M0950', 'C', 'k2'), ('M0951', 'C', 'n2'),
                  ('M0952', 'C', 'rate2'), ('M0953', 'C', 'p'),
                  ('M0956', 'P1', 'k'), ('M0957', 'P1', 'n'), ('M0959', 'P2', 'k'), ('M0960', 'P2', 'n'),
                  ('M0961', 'P2', 'rate'), ('M0963', 'P3', 'k'), ('M0964', 'P3', 'n'),
                  ('M0967', 'Q1_program', 'k'), ('M0968', 'Q1_program', 'n'),
                  ('M0970', 'Q2_program', 'k'), ('M0971', 'Q2_program', 'n'), ('M0972', 'Q2_program', 'rate'),
                  ('M0974', 'Q3_program', 'k'), ('M0975', 'Q3_program', 'n'), ('M0976', 'Q3_program', 'rate'),
                  ('M0979', 'Q1_clock', 'k'), ('M0980', 'Q1_clock', 'n'),
                  ('M0982', 'Q2_clock', 'k'), ('M0983', 'Q2_clock', 'n'), ('M0984', 'Q2_clock', 'rate'),
                  ('M0986', 'Q3_clock', 'k'), ('M0987', 'Q3_clock', 'n'), ('M0988', 'Q3_clock', 'rate')):
    E(cid, ('100*' if f.startswith('rate') else '') + T2.format(f=f, c=c))
E('M1013', "100*v('S15_confusion','between-cohort difference',scenario='weak shared effect')")
E('M1014', "100*v('S15_confusion','replicated',scenario='composition shift only')")

# ----------------------------------------------------------------------------- Supplementary Information
E('S1054', "tuple(col('S05_insample','n_in_matrix',module='All module genes'))")
E('S1055', "len(set(zip(rows('S05_insample',source='pkl_archive').outcome, rows('S05_insample',source='pkl_archive').module)))",
  '23 modules + all-module clock, two targets')
E('S1058', "sum(col('S03_modules','n_genes_serialized'))")
E('S1061', "tuple(col('S02_detected','composite_genes_serialized'))")
E('S1086', "sum(col('S03_modules','n_genes_serialized'))")
E('S1107', "med('S02_rle','genes_nonzero_all')")
E('S1108', "mn('S02_rle','genes_nonzero_all')")
E('S1109', "mx('S02_rle','genes_nonzero_all')")
E('S1117', "num(v('S01_cohorts','age',cohort='GSE136831 (Adams)',group='IPF'),2)", 'youngest IPF donor')
E('S1153', "num(v('S01_cohorts','note',cohort='HLCA analyzed'),2)", 'eligible donors, from the note of the HLCA row')
E('S1154', "num(v('S01_cohorts','age',cohort='HLCA analyzed'),0)")
E('S1155', "num(v('S01_cohorts','age',cohort='HLCA analyzed'),1)")
E('S1156', "v('S01_cohorts','n',cohort='HLCA analyzed')")
E('S1159', "v('S01_cohorts','n',cohort='GTEx v8 lung')")
E('S1165', "v('S18_positive','n_young_mice')")
E('S1172', "100*v('S18_positive','symbol_to_entrez_rate')")
E('S1173', "v('S18_positive','n_symbols')")
E('S1174', "v('S18_resolution','expected_shift',species='mouse')/21")
E('S1177', "v('S18_recovery_ci','expected_shift',unit='whole')")
E('S1199', "v('S13_replication','n_IPF_GSE136831',coarse='Macrophage')")
E('S1200', "v('S13_replication','n_Control_GSE136831',coarse='Macrophage')")
E('S1201', "v('S13_replication','n_IPF_GSE135893',coarse='Macrophage')")
E('S1202', "v('S13_replication','n_Control_GSE135893',coarse='Macrophage')")
E('S1216', "tuple(col('S15_confusion','tau_splithalf'))")
for cid, sc, field in (('S1217', 'shared effect', 'replicated'), ('S1218', 'sub-state composition', 'composition'),
                       ('S1219', 'confound in one cohort', 'between-cohort difference'),
                       ('S1220', 'no effect', 'inconclusive'), ('S1221', 'weak shared effect', 'replicated'),
                       ('S1222', 'weak shared effect', 'between-cohort difference'),
                       ('S1223', 'composition shift only', 'replicated'),
                       ('S1224', 'composition shift only', 'between-cohort difference'),
                       ('S1225', 'composition shift only', 'inconclusive')):
    E(cid, f"100*v('S15_confusion','{field}',scenario='{sc}')")
E('S1228', "v('S01_cohorts','n',cohort='HLCA analyzed')")
E('S1229', "v('S01_cohorts','n',cohort='HLCA 42-donor cap')")
E('S1230', "n('S12_sensitivity',analysis='main')", 'all planned tests (primary analysis rows of S12)')

# ----------------------------------------------------------------------------- kinds changed after the S14 correction
KIND.update({'M0745': 'computed', 'M0747': 'computed', 'M0748': 'computed'})
E('M0745', "n('S14_interaction',celltype='Macrophage')", 'module scores tested per cell type')
E('M0747', "n('S14_interaction',celltype='Monocyte')")
E('M0748', "n('S14_interaction')", 'size of the BH family (reproduced from p_interaction)')

# ----------------------------------------------------------------------------- reported values of the original study stored in the tables
E('M0283', "med('S05_serialized','r_reported',outcome='Mortality',query=\"module != 'allmodulegenes'\")")
E('M0286', "v('S05_serialized','r_reported',outcome='Mortality',module='allmodulegenes')")
E('M0287', "v('S05_serialized','r_reported',outcome='Chronological',module='allmodulegenes')")
E('M0306', "v('S05_serialized','r_reported',outcome='Mortality',module='allmodulegenes')")
E('M0307', "v('S05_serialized','r_reported',outcome='Chronological',module='allmodulegenes')")
E('M0826', "med('S05_serialized','r_reported',outcome='Mortality',query=\"module != 'allmodulegenes'\")")

# ----------------------------------------------------------------------------- numbers written as words (W ids, words_r2.csv)
HITS = "rows('S09_ipf',query='fdr_IPF < 0.10')"
W_COMPUTED = {
    'W0006': "sum(secs(a, b) > 0 for a, b in zip(rows('S21_plans').first_result_mtime, rows('S21_plans').write_time_utc))",
    'W0007': "n('S21_plans',query='amendment.notna()')",
    'W0008': "n('S16_program',met='not met',item=['P1','P2','P3'])",
    'W0010': "sum(col('S03_modules','n_genes_serialized')) == mx('S02_detected','module_genes_serialized')",
    'W0016': "tuple(num(v('S01_cohorts','note',cohort='Klotho knockout',group=g),0) for g in ('Kidney','Skeletal muscle'))",
    'W0017': "tuple(num(v('S01_cohorts','note',cohort='Klotho knockout',group=g),1) for g in ('Kidney','Skeletal muscle'))",
    'W0018': "n('S02_celltypes',atlas='GSE136831',in_statistics_set=True)",
    'W0020': f"v('S01_cohorts','n',cohort='{AD}',group='Control') - v('S01_cohorts','n',cohort='{GRID}',group='Control')",
    'W0021': f"v('S01_cohorts','n',cohort='{AD}',group='COPD') - v('S01_cohorts','n',cohort='{GRID}',group='COPD')",
    'W0022': "tuple(v('S02_celltypes','donors_IPF',atlas='GSE136831',celltype=c) for c in ('ATII','Lymphatic'))",
    'W0024': "len({x.rsplit(' ',1)[0] for g in ('Control','IPF') for x in "
             "v('S01_cohorts','tissue',cohort='GSE135893 (Habermann)',group=g).split('; ')[1].split(', ')})",
    'W0025': "v('S09_funnel','amended_r5',quantity='unestimable')",
    'W0039': "n('S16_program',met='not met',item=['P1','P2','P3'])",
    'W0050': "n('S05_coherence',source='table5b')",
    'W0051': "n('S05_serialized',query=\"module != 'allmodulegenes' and r_native < 0.9*r_reported\")",
    'W0054': "max(mx('S17_corrupt','n_program_pass',tag='corrupted'), mx('S17_corrupt','n_clock_pass',tag='corrupted'))",
    'W0055': "n('S16_program',met='met',item=['P1','P2','P3'])",
    'W0056': "n('S16_program',item=['P1','P2','P3'])",
    'W0064': "len(set(T('S08_adams_controls').celltype))",
    'W0065': "n('S08_gtex',outcome='Mortality',query=\"adj_fdr < 0.10 and module != 'Composite'\")",
    'W0067': "v('S18_contrib','n_ge_10',species='mouse',outcome='Mortality')",
    'W0068': "v('S18_positive','n_celltypes')",
    'W0069': "v('S16_criteria','k',criterion='C')",
    'W0070': "n('S18_tms',score='composite_Mortality',query=\"unit.str.startswith('celltype')\")",
    'W0071': "n('S11_components',level='component',sign='+',query='fdr < 0.10')",
    'W0072': "n('S11_components',level='component',sign='-',query='fdr < 0.10')",
    'W0073': "n('S11_celltype',query='fdr_S < 0.10')",
    'W0074': "n('S09_ipf',robust=True,IPF_specific=False)",
    'W0075': "n('S09_ipf',robust=True,celltype='ncMonocyte',query='fdr_IPF < 0.10')",
    'W0076': "tuple(sorted(set(rows('S09_ipf',celltype='ncMonocyte').n_Control)))",
    'W0077': "n('S09_ipf',IPF_specific=True)",
    'W0078': f"sum(v('S03_modules','heldout_r_Mortality',colour=m) < 0.30 for m in {HITS}.module)",
    'W0079': "n('S09_ipf',query='fdr_IPF < 0.10')",
    'W0080': "n('S12_sensitivity',analysis='main',query='fdr < 0.10') - n('S09_ipf',query='fdr_IPF < 0.10')",
    'W0081': "v('S09_funnel','amended_r5',quantity='unestimable')",
    'W0082': f"len({S12H('age_overlap')}.dropna(subset=['ratio']))",
    'W0083': "sum(min(v('S12_sensitivity','ratio',celltype=c,module=m,analysis=a) for a in "
             "('tech_covariates','depth_group_matched','cells_group_matched')) >= 0.5 "
             f"for c, m in zip({HITS}.celltype, {HITS}.module))",
    'W0084': "n('S09_ipf',query='fdr_IPF < 0.10')",
    'W0085': "sum(v('S12_sensitivity','ratio',celltype=c,module=m,analysis=k+'_equalized') < 0.5 <= "
             "v('S12_sensitivity','ratio',celltype=c,module=m,analysis=k+'_loss_control') "
             f"for c, m in zip({HITS}.celltype, {HITS}.module) for k in ('depth','cells'))",
    'W0086': f"len(rows('S10_within',query='fdr_hits_random < 0.10').merge({HITS}[['celltype','module']]))",
    'W0087': f"len(rows('S10_within',query='fdr_hits_perm < 0.10').merge({HITS}[['celltype','module']]))",
    'W0088': "n('S15_real',attribution='replicated')",
    'W0090': "v('S13_replication','n_IPF_GSE136831',coarse='AT2')",
    'W0091': "v('S13_replication','nk_lodo_n',coarse='NK')",
    'W0092': "1 if v('S13_replication','vanderbilt_equals_full',coarse='NK') else 0",
    'W0094': "n('S15_real',attribution='replicated')",
    'W0106': "v('S17_scrambled','original_submission_hits') / n('S09_ipf',query='fdr_IPF < 0.10')",
    'W0108': "0.60 <= v('S16_bias','rate_set_specific',table='by_class',scenario='S2a',d='all',"
             "class_or_denominator='false_hit') <= 0.70",
    'W0112': "n('S16_program',met='not met',item=['P1','P2','P3'])",
    'W0114': "v('S18_positive','n_celltypes')",
    'W0116': "v('S09_funnel','amended_r5',quantity='unestimable')",
    'W0117': "n('S13_replication')",
    'W0119': "n('S05_coherence',source='official',query='p < 0.002')",
    'W0120': "n('S05_coherence',source='table5b')",
    'W0122': "tuple(num(v('S01_cohorts','note',cohort='Klotho knockout',group=g),0) for g in ('Kidney','Skeletal muscle'))",
    'W0123': "tuple(num(v('S01_cohorts','note',cohort='Klotho knockout',group=g),1) for g in ('Kidney','Skeletal muscle'))",
    'W0125': "v('S18_positive','n_celltypes')",
    'W0128': "v('S13_replication','n_IPF_GSE136831',coarse='AT2')",
    'W0139': "sum(col('S03_modules','n_genes_serialized')) == mx('S02_detected','module_genes_serialized')",
    'W0146': "tuple(num(v('S01_cohorts','note',cohort='Klotho knockout',group=g),0) for g in ('Kidney','Skeletal muscle'))",
    'W0147': "tuple(num(v('S01_cohorts','note',cohort='Klotho knockout',group=g),1) for g in ('Kidney','Skeletal muscle'))",
    'W0161': "n('S16_program',met='not met',item=['P1','P2','P3'])",
}
for cid, expr in W_COMPUTED.items():
    E(cid, expr)
    KIND[cid] = 'computed'
# number words that are not quantities (pronouns, headings, rhetorical counts, layout of the files)
W_STRUCTURAL = {'W0002', 'W0029', 'W0031', 'W0053', 'W0061', 'W0100', 'W0105', 'W0132', 'W0160', 'W0178',
                'W0180', 'W0181'}

# ----------------------------------------------------------------------------- training data, second spike-in experiment and other checks
# S01_training: the released training matrix and the check of its mortality target (Supplementary Methods S1)
def TR(q):
    return f"v('S01_training','value',quantity='{q}')"
for cid, q in (('M0073', 'training_genes'), ('M0117', 'training_genes'), ('M0436', 'training_genes'),
               ('S1092', 'training_genes'), ('M0074', 'training_samples'), ('M0102', 'training_samples'),
               ('M0893', 'training_samples'), ('S1019', 'training_samples'), ('S1064', 'training_samples'),
               ('M0075', 'training_datasets'), ('M0894', 'training_datasets'), ('S1022', 'training_datasets'),
               ('M0076', 'training_samples_mouse'), ('S1020', 'training_samples_mouse'),
               ('M0077', 'training_samples_rat'), ('S1021', 'training_samples_rat'),
               ('M0088', 'composite_features_in_modules'), ('M0209', 'composite_features_in_modules'),
               ('M0090', 'composite_features_outside_modules'), ('M0211', 'composite_features_outside_modules'),
               ('M0098', 'control_mice_slope_per_month'), ('S1029', 'control_mice_slope_per_month'),
               ('M0099', 'control_mice_r_target_age'), ('S1176', 'control_mice_r_target_age'),
               ('M0100', 'control_mice'), ('S1175', 'control_mice'),
               ('S1018', 'samples_without_reference_group'), ('S1023', 'samples_with_reference_group'),
               ('S1024', 'target_reproduced_exactly'), ('S1027', 'control_animals_r_target_age'),
               ('S1028', 'control_animals'), ('S1031', 'control_animals_slope_per_0.1_lifespan'),
               ('S1032', 'interventions_diet_n'), ('S1033', 'interventions_diet_mean_target'),
               ('S1035', 'interventions_compounds_n'), ('S1036', 'interventions_compounds_mean_target'),
               ('S1038', 'interventions_genotypes_n'), ('S1039', 'interventions_genotypes_mean_target'),
               ('S1040', 'interventions_genotypes_median_target'), ('S1041', 'target_min'), ('S1042', 'target_max'),
               ('S1179', 'expected_chronological_shift_21_months')):
    E(cid, TR(q))
E('S1025', f"100*{TR('target_reproduced_exactly')}/{TR('training_samples')}")
E('S1026', f"100*{TR('target_reproduced_exactly')}/{TR('samples_with_reference_group')}")
E('S1034', f"10**{TR('interventions_diet_mean_target')}", 'hazard ratio = 10^(mean log10 target)')
E('S1037', f"10**{TR('interventions_compounds_mean_target')}", 'hazard ratio = 10^(mean log10 target)')

# donors (S01_cohorts, S02_excluded_donors)
E('M0147', "v('S02_excluded_donors','max_cells_in_grid_celltype',cohort='GSE136831',group='COPD')")
E('S1125', f"v('S01_cohorts','n',cohort='{AD}',group='donors with pseudobulks')")
E('S1160', "num(v('S01_cohorts','age',cohort='GTEx v8 lung'),0)", 'lower bound of the youngest age bracket')
E('S1161', "num(v('S01_cohorts','age',cohort='GTEx v8 lung'),1)", 'upper bound of the oldest age bracket')

# input check: gene-set overlap (S05_overlap), retraining check (S17_retrain_check), within-hit passes (S17_within_corrupt)
E('M0323', "med('S05_overlap','jaccard')")
E('S1089', "v('S17_retrain_check','spearman_rho_coefficients',module='red')")
E('M0335', "mn('S17_within_corrupt','mean_within_hit_passes',path='V1')")
E('M0336', "mx('S17_within_corrupt','mean_within_hit_passes',path='V1')")
E('M0337', "mn('S17_within_corrupt','mean_within_hit_passes',path='V2')")
E('M0338', "mx('S17_within_corrupt','mean_within_hit_passes',path='V2')")

# rodent benchmarks (S04_benchmarks) and Klotho noise (S07_klotho_noise)
for cid, clock, outcome in (('M0434', 'all_module_genes', 'Mortality'), ('M0435', 'all_module_genes', 'Chronological'),
                            ('M0437', 'all_genes', 'Mortality'), ('M0438', 'all_genes', 'Chronological')):
    E(cid, f"v('S04_benchmarks','r_oof_all',clock='{clock}',outcome='{outcome}')")
E('M0455', "v('S07_klotho_noise','median_abs_shift_6v6',tissue='Kidney')", 'wild-type halves, scaled to six versus six')
E('M0456', "v('S07_klotho_noise','median_abs_shift_6v6',tissue='Skeletal muscle')", 'wild-type halves, scaled to six versus six')

# second spike-in experiment: MDE80 (S16_mde2), pass rates (S16_pass2), within-hit labels (S16_within2),
# escalation of draws (S16_escalation2); first experiment, S1b nulls (S16_s1b); Table 2 criteria (S16_criteria)
MDE = "'S16_mde2','mde80',scenario='S1a',gate='{g}'"
for cid in ('M0003', 'M0399', 'M0776', 'M0790'):
    E(cid, "mn(" + MDE.format(g='program') + ")", 'smallest program-level MDE80')
for cid in ('M0004', 'M0400', 'M0777'):
    E(cid, "mx(" + MDE.format(g='program') + ")", 'largest program-level MDE80 that is reached at d <= 2')
for cid, ct in (('M0402', 'Macrophage_Alveolar'), ('M0403', 'Macrophage'), ('M0404', 'T'), ('M0405', 'cMonocyte'),
                ('M0406', 'T_Cytotoxic'), ('M0407', 'NK'), ('M0408', 'B'), ('M0409', 'Ciliated'), ('M0410', 'cDC2'),
                ('M0629', 'Macrophage'), ('M0630', 'Macrophage_Alveolar')):
    E(cid, "v(" + MDE.format(g='program') + f",celltype='{ct}')")
for cid, ct in (('M0411', 'Macrophage_Alveolar'), ('M0412', 'Macrophage'), ('M0413', 'cDC2')):
    E(cid, "v(" + MDE.format(g='clock') + f",celltype='{ct}')")
for cid, g in (('W0001', 'program'), ('W0057', 'program'), ('W0058', 'clock')):
    E(cid, f"n('S16_mde2',scenario='S1a',gate='{g}',query='mde80 <= 2')")
    KIND[cid] = 'computed'
WH = "v('S16_within2','{f}',definition='within_hit',population='false_hits')"
for cid in ('M0006', 'M0389', 'M0656'):
    E(cid, '100*' + WH.format(f='rate_weight_specific'))
for cid in ('M0007', 'M0388', 'M0655'):
    E(cid, '100*' + WH.format(f='rate_set_specific'))
E('M0370', "100*v('S16_pass2','rate',scenario='S1a',gate='program',d=2.0)")
E('M0371', "100*v('S16_pass2','rate',scenario='S1a',gate='clock',d=2.0)")
E('M0372', "100*v('S16_pass2','rate',scenario='S1b',gate='program',d=2.0)")
E('M0373', "100*v('S16_s1b','rate_random_only',d=2.0)")
E('M0375', "100*v('S16_s1b','rate_both',d=2.0)")
E('M0393', "100*v('S16_escalation2','rate_extended')")
E('M0396', "v('S16_escalation2','min_screen_fdr_not_extended') >= 0.10",
  'no pair with a first-stage p >= 0.05 reached a screen FDR below 0.10')
for cid, expr in (('W0101', "n('S16_criteria',met='Yes')"), ('W0102', "n('S16_criteria')"),
                  ('W0103', "n('S16_criteria',met='Yes')"), ('W0104', "n('S16_criteria',met='No')")):
    E(cid, expr)
    KIND[cid] = 'computed'

# resolution: HLCA expression-only reconstruction (S18_resolution), Angelidis grid rule (S18_positive)
E('M0562', "v('S18_resolution','expression_only_rho_adj',species='human')")
E('M0563', "v('S18_resolution','expression_only_p',species='human')")
E('W0168', "v('S18_positive','n_celltypes_min4')")
KIND['W0168'] = 'computed'

# all 253 planned tests (S09_planned)
E('M0632', "n('S09_planned')")
E('M0634', "n('S09_planned',hit_253=True)")
E('M0637', "v('S09_planned','fdr_IPF_253',celltype='NK',module='greenyellow')", 'FDR among all 253 planned tests')
E('S1230', "n('S09_planned')", 'all planned tests')
