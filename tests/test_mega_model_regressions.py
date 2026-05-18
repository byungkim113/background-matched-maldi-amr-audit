import ast
import pathlib
import re
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
MEGA_PATH = ROOT / "Mega_Model.py"
AUDIT_PATH = ROOT / "scripts" / "ecoli_drug_panel_audit.py"
CONTRASTIVE_PATH = ROOT / "scripts" / "background_matched_contrastive_kaggle.py"


def read_source():
    return MEGA_PATH.read_text(encoding="utf-8")


def parsed_tree():
    return ast.parse(read_source())


def get_constant_assignment(name):
    for node in parsed_tree().body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == name:
                    return ast.literal_eval(node.value)
        if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            if node.target.id == name:
                return ast.literal_eval(node.value)
    raise AssertionError(f"{name} assignment not found")


def get_function_source(name):
    source = read_source()
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return ast.get_source_segment(source, node)
    raise AssertionError(f"{name} function not found")


class MegaModelRegressionTests(unittest.TestCase):
    def test_mega_model_file_exists(self):
        self.assertTrue(MEGA_PATH.exists())

    def test_ecoli_drug_panel_audit_script_exists_and_is_no_training(self):
        if not AUDIT_PATH.exists():
            self.skipTest("ecoli_drug_panel_audit.py is not part of this repository snapshot")
        source = AUDIT_PATH.read_text()
        self.assertIn("ECOLI_DRUG_PANEL", source)
        self.assertIn("compute_pairwise_label_correlations", source)
        self.assertIn("ecoli_drug_panel_inventory.csv", source)
        self.assertIn("ecoli_drug_panel_correlations.csv", source)
        self.assertIn("ecoli_drug_panel_audit.md", source)
        self.assertIn("discover_antibiotic_columns", source)
        self.assertIn("collect_all_drug_labels", source)
        self.assertIn("compute_all_drug_eligibility", source)
        self.assertIn("all_drug_inventory.csv", source)
        self.assertIn("all_drug_eligibility.csv", source)
        self.assertIn("organism_drug_site_matrix.csv", source)
        self.assertIn("ecoli_candidate_panel.csv", source)
        self.assertIn("gram_negative_candidate_panel.csv", source)
        self.assertIn("mechanism_validation_candidate_panel.csv", source)
        self.assertIn("co_resistance_blocks.csv", source)
        self.assertIn("discovered_antibiotic_columns.csv", source)
        self.assertIn("parse_known_args", source)
        self.assertNotIn("torch", source)
        self.assertNotIn("DataLoader", source)
        self.assertNotIn("fit(", source)

    def test_background_matched_contrastive_script_is_eval_only(self):
        self.assertTrue(CONTRASTIVE_PATH.exists())
        source = CONTRASTIVE_PATH.read_text()
        self.assertIn("background_matched_predictions.csv", source)
        self.assertIn("background_matched_contrastive_summary.csv", source)
        self.assertIn("background_matched_contrastive_summary.md", source)
        self.assertIn("torch.load", source)
        self.assertIn("create_maldi_model", source)
        self.assertIn("parse_known_args", source)
        self.assertIn("num_workers=0", source)
        self.assertNotIn("train_all_seeds", source)
        self.assertNotIn("optimizer", source)
        self.assertNotIn("backward()", source)

    def test_background_matched_contrastive_script_matches_on_coresistance(self):
        source = CONTRASTIVE_PATH.read_text()
        self.assertIn("background_signature", source)
        self.assertIn("exclude_focal", source)
        self.assertIn("stratum_centered_auc", source)
        self.assertIn("pairwise_accuracy", source)
        self.assertIn("matched_retention", source)
        self.assertIn("min_pos_per_stratum", source)
        self.assertIn("min_neg_per_stratum", source)

    def test_background_matched_contrastive_script_autodetects_kaggle_input_runs(self):
        source = CONTRASTIVE_PATH.read_text()
        self.assertIn("resolve_run_dir", source)
        self.assertIn("discover_run_dirs", source)
        self.assertIn("rglob", source)
        self.assertIn("/kaggle/input", source)
        self.assertIn("bfdf121", source)
        self.assertIn("newruns/runs", source)
        self.assertIn("/kaggle/working/background_matched_contrastive", source)
        self.assertIn("config.json not found", source)
        self.assertIn("exp_ecoli_mechanism6_drugid_mae30", source)
        self.assertIn("newruns", source)
        self.assertIn("candidate run directories discovered", source)

    def test_background_matched_contrastive_script_autodetects_mega_model_file(self):
        source = CONTRASTIVE_PATH.read_text()
        self.assertIn("resolve_mega_model_path", source)
        self.assertIn("/kaggle/working/Mega_Model.py", source)
        self.assertIn("Mega_Model.py not found", source)
        self.assertIn("candidate Mega_Model.py paths", source)
        self.assertIn("Please upload or copy Mega_Model.py", source)

    def test_background_matched_contrastive_script_reports_empty_row_boundaries(self):
        source = CONTRASTIVE_PATH.read_text()
        self.assertIn("resolve_data_root", source)
        self.assertIn("candidate DRIAMS data roots", source)
        self.assertIn("diagnostics", source)
        self.assertIn("spectrum_files", source)
        self.assertIn("rows_with_spectrum", source)
        self.assertIn("drug_columns", source)

    def test_background_matched_contrastive_script_handles_single_stratum_column(self):
        source = CONTRASTIVE_PATH.read_text()
        self.assertIn("valid_indices", source)
        self.assertIn("group.index", source)
        self.assertNotIn("key_frame.loc[valid_keys]", source)

    def test_background_matched_contrastive_script_can_reuse_predictions_csv(self):
        source = CONTRASTIVE_PATH.read_text()
        self.assertIn("--predictions-csv", source)
        self.assertIn("--force-rescore", source)
        self.assertIn("Reusing scored predictions", source)
        self.assertIn("background_matched_predictions.csv", source)

    def test_background_matched_contrastive_script_adds_statistical_audit_outputs(self):
        source = CONTRASTIVE_PATH.read_text()
        self.assertIn("bootstrap_metric_ci", source)
        self.assertIn("permutation_null_within_strata", source)
        self.assertIn("assign_adequacy_label", source)
        self.assertIn("raw_auc_ci_low", source)
        self.assertIn("raw_auc_ci_high", source)
        self.assertIn("matched_auc_ci_low", source)
        self.assertIn("matched_auc_ci_high", source)
        self.assertIn("stratum_centered_auc_ci_low", source)
        self.assertIn("stratum_centered_auc_ci_high", source)
        self.assertIn("permutation_p", source)
        self.assertIn("adequacy_label", source)
        self.assertIn("--bootstrap-n", source)
        self.assertIn("--permutation-n", source)
        self.assertIn("--adequacy-min-n-matched", source)
        self.assertIn("--adequacy-min-retention", source)

    def test_background_matched_contrastive_script_writes_sensitivity_outputs(self):
        source = CONTRASTIVE_PATH.read_text()
        self.assertIn("parse_sensitivity_thresholds", source)
        self.assertIn("run_sensitivity_summaries", source)
        self.assertIn("background_matched_sensitivity.csv", source)
        self.assertIn("background_matched_sensitivity.md", source)
        self.assertIn("--sensitivity-thresholds", source)

    def test_background_matched_contrastive_script_can_score_lgbm(self):
        source = CONTRASTIVE_PATH.read_text()
        self.assertIn("--score-family", source)
        self.assertIn("--lgbm-model-dir", source)
        self.assertIn("--row-template-csv", source)
        self.assertIn("load_lgbm_models", source)
        self.assertIn("score_lgbm_rows", source)
        self.assertIn("background_matched_lgbm_single_predictions.csv", source)
        self.assertIn("background_matched_lgbm_multi_predictions.csv", source)
        self.assertIn("lgbm_single_org", source)
        self.assertIn("lgbm_multi.pkl", source)

    def test_pair_profiles_include_run14_clinical4_clinical5_and_gram_negative6(self):
        profiles = get_constant_assignment("PAIR_PROFILES")

        self.assertEqual(
            profiles["run14"],
            [
                ("Escherichia coli", "Ciprofloxacin"),
                ("Staphylococcus aureus", "Oxacillin"),
            ],
        )
        self.assertEqual(
            profiles["clinical5"],
            [
                ("Escherichia coli", "Ciprofloxacin"),
                ("Escherichia coli", "Amoxicillin-Clavulanic acid"),
                ("Staphylococcus aureus", "Oxacillin"),
                ("Staphylococcus aureus", "Penicillin"),
                ("Staphylococcus epidermidis", "Erythromycin"),
            ],
        )
        self.assertEqual(
            profiles["clinical4"],
            [
                ("Escherichia coli", "Ciprofloxacin"),
                ("Escherichia coli", "Amoxicillin-Clavulanic acid"),
                ("Staphylococcus aureus", "Oxacillin"),
                ("Staphylococcus epidermidis", "Erythromycin"),
            ],
        )
        self.assertEqual(
            profiles["gram_negative6"],
            [
                ("Escherichia coli", "Ciprofloxacin"),
                ("Escherichia coli", "Ceftriaxone"),
                ("Escherichia coli", "Ceftazidime"),
                ("Klebsiella pneumoniae", "Ciprofloxacin"),
                ("Klebsiella pneumoniae", "Ceftriaxone"),
                ("Staphylococcus aureus", "Oxacillin"),
            ],
        )

    def test_ecoli_drug_panel_profile_and_smiles_are_present(self):
        profiles = get_constant_assignment("PAIR_PROFILES")
        self.assertEqual(
            profiles["ecoli_drug_panel"],
            [
                ("Escherichia coli", "Ciprofloxacin"),
                ("Escherichia coli", "Amoxicillin-Clavulanic acid"),
                ("Escherichia coli", "Ceftriaxone"),
                ("Escherichia coli", "Ceftazidime"),
                ("Escherichia coli", "Piperacillin-Tazobactam"),
                ("Escherichia coli", "Gentamicin"),
                ("Escherichia coli", "Trimethoprim-Sulfamethoxazole"),
            ],
        )
        source = read_source()
        self.assertIn("DRUG_SMILES", source)
        self.assertIn("DRUG_FINGERPRINT_DIM", source)
        for drug in [drug for _, drug in profiles["ecoli_drug_panel"]]:
            self.assertIn(drug, source)

    def test_ecoli_mechanism6_profile_is_available(self):
        profiles = get_constant_assignment("PAIR_PROFILES")
        self.assertEqual(
            profiles["ecoli_mechanism6"],
            [
                ("Escherichia coli", "Ciprofloxacin"),
                ("Escherichia coli", "Norfloxacin"),
                ("Escherichia coli", "Amoxicillin-Clavulanic acid"),
                ("Escherichia coli", "Ceftriaxone"),
                ("Escherichia coli", "Ceftazidime"),
                ("Escherichia coli", "Cefepime"),
            ],
        )
        source = read_source()
        self.assertIn("Norfloxacin", source)
        self.assertIn("Cefepime", source)
        self.assertIn("ecoli_mechanism6", source)

    def test_notebook_defaults_saureus_panel_oxa_background(self):
        self.assertEqual(get_constant_assignment("NOTEBOOK_PAIR_PROFILE"), "saureus_panel")
        self.assertEqual(
            get_constant_assignment("NOTEBOOK_EXPERIMENT"),
            "exp_saureus_panel_oxa_background_mae30",
        )
        self.assertIs(get_constant_assignment("NOTEBOOK_NO_BN_ADAPT"), True)
        self.assertEqual(get_constant_assignment("NOTEBOOK_EARLY_STOP"), "macro")
        self.assertEqual(get_constant_assignment("NOTEBOOK_SEED_POLICY"), "all")
        self.assertEqual(get_constant_assignment("NOTEBOOK_PREVALENCE_SHIFT"), "none")
        self.assertEqual(get_constant_assignment("MAE_EPOCHS"), 30)
        self.assertEqual(get_constant_assignment("NOTEBOOK_MAE_EPOCHS"), 30)
        self.assertIs(get_constant_assignment("NOTEBOOK_NO_LGBM"), True)
        self.assertIs(get_constant_assignment("NOTEBOOK_NO_SALIENCY"), True)

    def test_cli_defaults_match_clinical4_diagnostic(self):
        source = read_source()
        self.assertIn('p.add_argument("--experiment", default=NOTEBOOK_EXPERIMENT', source)
        self.assertIn('p.add_argument("--no-bn-adapt", action="store_true", default=NOTEBOOK_NO_BN_ADAPT', source)
        self.assertIn('default=NOTEBOOK_SEED_POLICY', source)

    def test_dann_lambda_is_not_applied_twice(self):
        body = get_function_source("run_epoch_dann")
        self.assertIn("model.forward_dann(x_s, org_s, lam)", body)
        self.assertIn("model.forward_domain_only(x_d, lam)", body)
        self.assertIn("loss = res_loss + dom_loss_s + dom_loss_t", body)
        self.assertNotRegex(body, r"loss\s*=\s*res_loss\s*\+\s*lam\s*\*")

    def test_main_training_early_stops_on_macro_auc(self):
        body = get_function_source("train_all_seeds")
        self.assertRegex(body, r"stop_val\s*=\s*macro")
        self.assertNotIn("stop_val = primary if not math.isnan(primary) else macro", body)

    def test_prevalence_shift_uses_correct_threshold_direction_and_cap(self):
        shift_body = get_function_source("capped_logit_shift")
        prevalence_body = get_function_source("_prevalence_shift")
        self.assertIn("MAX_PREVALENCE_ODDS_SHIFT", read_source())
        self.assertIn("np.clip(shift", shift_body)
        self.assertIn("logit_base - logit_shift", prevalence_body)
        self.assertNotIn("logit_base + shift", prevalence_body)

    def test_active_pair_guard_reserves_validation_resistant_samples(self):
        source = read_source()
        self.assertIn("required_pretest_resistant", source)
        self.assertIn("MIN_R_TRAIN + MIN_R_VAL", source)
        self.assertIn("has_enough_pretest_resistance", source)

    def test_seed_quality_control_is_present(self):
        source = read_source()
        self.assertIn("SEED_POLICY", source)
        self.assertIn("TOP_K_SEEDS", source)
        self.assertIn("MIN_SEED_MACRO_AUC", source)
        self.assertIn("select_seed_ensemble", source)

    def test_run14_overlap_summary_is_present(self):
        source = read_source()
        self.assertIn("RUN14_OVERLAP_PAIRS", source)
        self.assertIn("compute_subset_metrics", source)
        self.assertIn("run14_overlap_summary.csv", source)

    def test_same_site_empty_samples_are_skipped(self):
        body = get_function_source("ensemble_predict")
        self.assertRegex(body, r"if\s+not\s+samples:")
        self.assertIn("return np.array([]), np.array([]), np.array([])", body)

    def test_mae_pretraining_excludes_same_site_test_samples(self):
        body = get_function_source("train_all_seeds")
        self.assertIn("mae_paths = list({p for p, _, _ in train_s + val_s})", body)
        self.assertNotIn("train_s + val_s + test_s", body)

    def test_missing_reporting_artifacts_are_present(self):
        source = read_source()
        self.assertIn("compute_per_organism_summary", source)
        self.assertIn("per_organism_summary.csv", source)
        self.assertIn("compute_threshold_sanity_report", source)
        self.assertIn("threshold_sanity_report.csv", source)
        self.assertIn("seed_summary.csv", source)

    def test_protocol_toggles_are_exposed_on_cli(self):
        source = read_source()
        for flag in (
            "--early-stop",
            "--seed-policy",
            "--top-k-seeds",
            "--min-seed-macro-auc",
            "--prevalence-shift",
            "--no-prevalence-shift",
            "--strict-mae-source-only",
        ):
            self.assertIn(flag, source)

    def test_drug_conditioning_toggles_are_exposed_and_wired(self):
        source = read_source()
        run_body = get_function_source("run_experiment")
        train_body = get_function_source("train_all_seeds")
        self.assertIn("DRUG_CONDITIONING_CHOICES", source)
        self.assertIn("NOTEBOOK_DRUG_CONDITIONING", source)
        self.assertIn("--drug-conditioning", source)
        self.assertIn("--with-leave-one-drug-out", source)
        self.assertIn("args.drug_conditioning", run_body)
        self.assertIn("drug_conditioning=", train_body)
        self.assertIn("leave_one_drug_out.csv", source)

    def test_drug_conditioned_model_uses_shared_head_and_fingerprint_film(self):
        source = read_source()
        self.assertIn("class DrugConditionedFiLMLayer", source)
        self.assertIn("class MALDICNNDrugConditioned", source)
        self.assertIn("self.res_head = nn.Linear(256, 1)", source)
        self.assertIn("drug_features", source)
        self.assertIn("create_maldi_model", source)
        body = get_function_source("create_maldi_model")
        self.assertIn("MALDICNNDrugConditioned", body)
        self.assertIn('"task_id"', body)

    def test_toggles_are_wired_into_behavior(self):
        train_body = get_function_source("train_all_seeds")
        run_body = get_function_source("run_experiment")
        site_thresh_body = get_function_source("site_thresholds_for")

        self.assertIn("early_stop=", train_body)
        self.assertIn('early_stop == "primary"', train_body)
        self.assertIn("seed_policy=", train_body)
        self.assertIn("strict_mae_source_only=", train_body)
        self.assertIn("prevalence_shift_mode", site_thresh_body)
        self.assertIn('prevalence_shift_mode == "none"', site_thresh_body)
        self.assertIn("args.strict_mae_source_only", run_body)
        self.assertIn("args.prevalence_shift", run_body)
        self.assertIn("args.seed_policy", run_body)

    def test_seed_warning_and_uda_labeling_are_explicit(self):
        source = read_source()
        self.assertIn("NEAR_CHANCE_MACRO_AUC", source)
        self.assertIn("near chance", source)
        self.assertIn("unsupervised domain adaptation", source)

    def test_review_patch_gpu_cache_lookup_avoids_eager_default(self):
        source = read_source()
        self.assertNotIn("cnn_cache.get(\n                site_label,\n                ensemble_predict", source)
        self.assertIn("cnn_cache[site_label] if site_label in cnn_cache", source)

    def test_review_patch_ablation_uses_screened_active_pairs(self):
        body = get_function_source("_train_ablation_variant")
        run_body = get_function_source("run_experiment")
        self.assertIn("active_pairs=None", body)
        self.assertIn("if active_pairs is None:", body)
        self.assertIn("active_pairs=active_pairs", run_body)

    def test_review_patch_mae_mask_is_vectorized(self):
        body = get_function_source("pretrain_mae")
        self.assertIn("mask.scatter_", body)
        self.assertNotIn("for i in range(B):", body)

    def test_review_patch_prevalence_zero_not_falsy(self):
        body = get_function_source("site_thresholds_for")
        self.assertIn("train_prev is not None", body)
        self.assertIn("target_prev is not None", body)
        self.assertNotIn("if train_prev and target_prev", body)

    def test_review_patch_debug_warnings(self):
        source = read_source()
        self.assertIn("MIN_TOTAL_SOURCE_SAMPLES_WARN", source)
        self.assertIn("Suspiciously small source sample count", source)
        self.assertIn("DANN disabled", source)

    def test_load_spectrum_has_clear_shape_and_length_guards(self):
        body = get_function_source("load_spectrum")
        self.assertIn("data.ndim != 2", body)
        self.assertIn("data.shape[1] < 2", body)
        self.assertIn("data.shape[0] < N_BINS", body)
        self.assertIn("expected at least {N_BINS}", body)
        self.assertIn("data[:N_BINS, 1]", body)
        self.assertIn("ValueError", body)

    def test_config_json_uses_runtime_seed_args_not_stale_globals(self):
        body = get_function_source("run_experiment")
        config_block = body[body.index("with open(out_dir / \"config.json\""):]
        self.assertIn("seed_policy=args.seed_policy", config_block)
        self.assertIn("top_k_seeds=args.top_k_seeds", config_block)
        self.assertIn("min_seed_macro_auc=args.min_seed_macro_auc", config_block)
        self.assertNotIn("seed_policy=SEED_POLICY", config_block)
        self.assertNotIn("top_k_seeds=TOP_K_SEEDS", config_block)
        self.assertNotIn("min_seed_macro_auc=MIN_SEED_MACRO_AUC", config_block)
        self.assertNotIn("requested_seed_policy", config_block)

    def test_ablation_training_respects_early_stop_config(self):
        body = get_function_source("_train_ablation_variant")
        run_body = get_function_source("run_experiment")
        self.assertIn('early_stop="macro"', body)
        self.assertIn('if early_stop == "primary":', body)
        self.assertIn("stop_val = primary", body)
        self.assertIn("stop_val = macro", body)
        self.assertIn("early_stop=args.early_stop", run_body)

    def test_scipy_import_is_optional_for_cli_smoke(self):
        source = read_source()
        temp_body = get_function_source("learn_ensemble_temperature")
        self.assertIn("scipy_optimize = None", source)
        self.assertIn("if scipy_optimize is None:", temp_body)
        self.assertIn("np.linspace", temp_body)

    def test_sklearn_metrics_import_is_optional_for_cli_smoke(self):
        source = read_source()
        self.assertIn("sklearn_metrics_available = False", source)
        self.assertIn("def roc_auc_score", source)
        self.assertIn("def average_precision_score", source)
        self.assertIn("def roc_curve", source)

    def test_paper_core_results_table_artifacts_are_present(self):
        source = read_source()
        self.assertIn("build_core_results_table", source)
        self.assertIn("core_results_table.csv", source)
        self.assertIn("core_results_table.md", source)
        self.assertIn("cnn_minus_best_lgbm_auc", source)
        self.assertIn("lgbm_single_auc", source)
        self.assertIn("lgbm_multi_auc", source)

    def test_saliency_interpretability_artifacts_are_present(self):
        source = read_source()
        self.assertIn("compute_occlusion_saliency", source)
        self.assertIn("saliency_top_bins.csv", source)
        self.assertIn("saliency_top_bins.md", source)
        self.assertIn("SALIENCY_WINDOW", source)
        self.assertIn("SALIENCY_MAX_SAMPLES_PER_PAIR", source)
        self.assertIn("--no-saliency", source)

    def test_saliency_uses_bounded_window_occlusion(self):
        body = get_function_source("compute_occlusion_saliency")
        self.assertIn("range(0, N_BINS, SALIENCY_STRIDE)", body)
        self.assertIn("SALIENCY_WINDOW", body)
        self.assertIn("SALIENCY_MAX_SAMPLES_PER_PAIR", body)
        self.assertIn("torch.no_grad()", body)
        self.assertEqual(get_constant_assignment("SALIENCY_WINDOW"), 200)
        self.assertEqual(get_constant_assignment("SALIENCY_STRIDE"), 100)
        self.assertEqual(get_constant_assignment("SALIENCY_MAX_SAMPLES_PER_PAIR"), 50)

    def test_saliency_mz_axis_maps_driams_bin_indices_to_da(self):
        body = get_function_source("load_mz_axis")
        self.assertIn("2000.0 + mz * 3.0", body)
        self.assertIn("2000.0 + np.arange(N_BINS", body)

    def test_mechanism_detectability_framing_is_present(self):
        source = read_source()
        self.assertIn("MECHANISM_DETECTABILITY", source)
        self.assertIn("expected_detectability", source)
        self.assertIn("chromosomal_structural", source)
        self.assertIn("mobile_heterogeneous", source)
        self.assertIn("mechanism_detectability_summary.csv", source)
        self.assertIn("mechanism_framing.md", source)

    def test_mechanism_summary_groups_more_and_less_detectable_pairs(self):
        source = read_source()
        body = get_function_source("compute_mechanism_detectability_summary")
        self.assertIn('("Escherichia coli", "Ciprofloxacin")', source)
        self.assertIn('("Staphylococcus aureus", "Oxacillin")', source)
        self.assertIn('("Escherichia coli", "Amoxicillin-Clavulanic acid")', source)
        self.assertIn('("Staphylococcus epidermidis", "Erythromycin")', source)
        self.assertIn("more_detectable", body)
        self.assertIn("less_detectable", body)
        self.assertIn("mean_auc", body)

    def test_optional_mae_epoch_performance_run_is_configurable(self):
        source = read_source()
        run_body = get_function_source("run_experiment")
        self.assertIn("NOTEBOOK_MAE_EPOCHS", source)
        self.assertIn("--mae-epochs", source)
        self.assertIn("args.mae_epochs", run_body)
        self.assertIn("global MAE_EPOCHS", run_body)

    def test_cross_site_saliency_stability_outputs_are_present(self):
        source = read_source()
        self.assertIn("compute_saliency_stability_summary", source)
        self.assertIn("saliency_top_bins_by_site.csv", source)
        self.assertIn("saliency_stability_summary.csv", source)
        self.assertIn("saliency_stability.md", source)
        self.assertIn("jaccard_top_bins", source)
        self.assertIn("mean_jaccard", source)

    def test_driams_d_shift_analysis_is_present(self):
        source = read_source()
        self.assertIn("SITE_CONTEXT", source)
        self.assertIn("compute_d_site_shift_analysis", source)
        self.assertIn("d_site_shift_analysis.csv", source)
        self.assertIn("d_site_shift_analysis.md", source)
        self.assertIn("community_acquired", source)
        self.assertIn("DRIAMS-D", source)

    def test_temporal_vs_random_cv_inflation_output_is_present(self):
        source = read_source()
        self.assertIn("NOTEBOOK_WITH_RANDOM_CV", source)
        self.assertIn("run_random_cv_inflation_analysis", source)
        self.assertIn("temporal_vs_random_cv.csv", source)
        self.assertIn("temporal_vs_random_cv.md", source)
        self.assertIn("auc_inflation", source)
        self.assertIn("--with-random-cv", source)
        self.assertIn("--no-random-cv", source)

    def test_requested_statistical_reporting_artifacts_are_present(self):
        source = read_source()
        self.assertIn("STAT_BOOTSTRAP_N", source)
        self.assertIn("compute_prediction_statistical_report", source)
        self.assertIn("auc_bootstrap_ci.csv", source)
        self.assertIn("macro_bootstrap_ci.csv", source)
        self.assertIn("mechanism_contrast_bootstrap_ci.csv", source)
        self.assertIn("mechanism_permutation_tests.csv", source)
        self.assertIn("bootstrap_ci_low", source)
        self.assertIn("permutation_p", source)
        self.assertIn("multiple-comparison caution", source)

    def test_requested_evaluation_and_confound_artifacts_are_present(self):
        source = read_source()
        self.assertIn("build_evaluation_critique_table", source)
        self.assertIn("evaluation_critique_table.csv", source)
        self.assertIn("evaluation_critique_table.md", source)
        self.assertIn("compute_mechanism_confound_checks", source)
        self.assertIn("mechanism_confound_checks.csv", source)
        self.assertIn("leave_one_pair_out", source)
        self.assertIn("excluding_DRIAMS-D", source)
        self.assertIn("prevalence", source)
        self.assertIn("sample_size", source)

    def test_requested_saliency_robustness_artifacts_are_present(self):
        source = read_source()
        self.assertIn("SALIENCY_STABILITY_TOP_KS", source)
        self.assertIn("SALIENCY_EXCLUDE_LOW_BINS_BELOW", source)
        self.assertIn("compute_saliency_robustness_summary", source)
        self.assertIn("saliency_stability_robustness.csv", source)
        self.assertIn("saliency_auc_alignment.csv", source)
        self.assertIn("saliency_null_control.csv", source)
        self.assertIn("saliency_stability_robustness_including_low_bins.csv", source)
        self.assertIn("saliency_stability_robustness_excluding_low_bins.csv", source)
        self.assertIn("top10_jaccard", source)
        self.assertIn("rank_correlation", source)
        self.assertIn("broad_window_jaccard", source)


if __name__ == "__main__":
    unittest.main()
