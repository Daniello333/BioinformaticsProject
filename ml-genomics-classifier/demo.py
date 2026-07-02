"""Demo script: load the sample genotype dataset, train logistic regression and
random forest classifiers to predict disease status from SNP genotypes, evaluate
them, and check whether the models recover the truly causal SNPs."""
from __future__ import annotations

import os

import pandas as pd
from sklearn.model_selection import train_test_split

from mlgenomics import (
    evaluate_model,
    feature_importance,
    plot_feature_importance,
    plot_roc_curve,
    train_logistic_regression,
    train_random_forest,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "output")


def main() -> None:
    genotypes = pd.read_csv(os.path.join(DATA_DIR, "genotypes.csv"), index_col=0)
    labels = pd.read_csv(os.path.join(DATA_DIR, "labels.csv"), index_col=0).iloc[:, 0]
    causal_snps = set(pd.read_csv(os.path.join(DATA_DIR, "causal_snps.csv"))["causal_snp_id"])

    print(f"Loaded {genotypes.shape[0]} samples x {genotypes.shape[1]} SNPs")
    print(f"Cases: {labels.sum()}, Controls: {(labels == 0).sum()}")
    print(f"Truly causal SNPs ({len(causal_snps)}): {sorted(causal_snps)}\n")

    X_train, X_test, y_train, y_test = train_test_split(
        genotypes, labels, test_size=0.25, stratify=labels, random_state=0
    )

    models = {
        "Logistic Regression (L1)": train_logistic_regression(X_train, y_train),
        "Random Forest": train_random_forest(X_train, y_train),
    }

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    results_by_model = {}

    for name, model in models.items():
        results = evaluate_model(model, X_test, y_test)
        results_by_model[name] = results
        print(f"=== {name} ===")
        print(f"Accuracy: {results['accuracy']:.3f}")
        print(f"ROC AUC:  {results['roc_auc']:.3f}")
        print(f"Confusion matrix [[TN, FP], [FN, TP]]: {results['confusion_matrix']}")

        importance = feature_importance(model, list(genotypes.columns))
        top10 = importance.head(10)
        recovered = sum(1 for snp in top10.index if snp in causal_snps)
        print(f"Causal SNPs recovered in top 10 features: {recovered}/10")
        print(f"Top 10 SNPs: {list(top10.index)}\n")

        safe_name = name.lower().replace(" ", "_").replace("(", "").replace(")", "")
        plot_feature_importance(
            importance, os.path.join(OUTPUT_DIR, f"importance_{safe_name}.png"), causal_snps=causal_snps
        )

    plot_roc_curve(results_by_model, os.path.join(OUTPUT_DIR, "roc_curve.png"))
    print(f"Plots written to {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
