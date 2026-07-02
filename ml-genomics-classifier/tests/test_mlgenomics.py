import os
import sys

from sklearn.model_selection import train_test_split

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from mlgenomics import (
    evaluate_model,
    feature_importance,
    simulate_genotypes,
    train_logistic_regression,
    train_random_forest,
)


def test_simulate_genotypes_shape_and_range():
    genotypes, labels, causal_snps = simulate_genotypes(
        n_samples=200, n_snps=50, n_causal_snps=5, seed=1
    )
    assert genotypes.shape == (200, 50)
    assert len(labels) == 200
    assert set(genotypes.to_numpy().ravel()).issubset({0, 1, 2})
    assert set(labels.unique()).issubset({0, 1})
    assert len(causal_snps) == 5
    assert all(snp in genotypes.columns for snp in causal_snps)


def test_logistic_regression_recovers_causal_snps():
    genotypes, labels, causal_snps = simulate_genotypes(
        n_samples=1500, n_snps=200, n_causal_snps=10, effect_size_range=(1.0, 2.0), seed=3
    )
    X_train, X_test, y_train, y_test = train_test_split(
        genotypes, labels, test_size=0.25, stratify=labels, random_state=0
    )
    model = train_logistic_regression(X_train, y_train)
    results = evaluate_model(model, X_test, y_test)
    assert results["roc_auc"] > 0.7

    importance = feature_importance(model, list(genotypes.columns))
    top20 = set(importance.head(20).index)
    recovered = len(top20 & set(causal_snps))
    assert recovered >= 7  # most causal SNPs should surface near the top


def test_random_forest_beats_chance():
    genotypes, labels, causal_snps = simulate_genotypes(
        n_samples=1000, n_snps=150, n_causal_snps=8, effect_size_range=(1.0, 2.0), seed=5
    )
    X_train, X_test, y_train, y_test = train_test_split(
        genotypes, labels, test_size=0.25, stratify=labels, random_state=0
    )
    model = train_random_forest(X_train, y_train, n_estimators=100)
    results = evaluate_model(model, X_test, y_test)
    assert results["roc_auc"] > 0.6


def test_evaluate_model_keys():
    genotypes, labels, _ = simulate_genotypes(n_samples=300, n_snps=30, n_causal_snps=3, seed=2)
    X_train, X_test, y_train, y_test = train_test_split(
        genotypes, labels, test_size=0.3, stratify=labels, random_state=0
    )
    model = train_logistic_regression(X_train, y_train)
    results = evaluate_model(model, X_test, y_test)
    assert set(results.keys()) == {"accuracy", "roc_auc", "confusion_matrix", "fpr_tpr"}
    assert 0.0 <= results["accuracy"] <= 1.0
    assert 0.0 <= results["roc_auc"] <= 1.0
