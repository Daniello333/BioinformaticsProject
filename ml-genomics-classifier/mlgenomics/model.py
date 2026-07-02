"""Model training and evaluation for genotype-based disease classification."""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

def train_logistic_regression(X_train: pd.DataFrame, y_train: pd.Series, l1: bool = True) -> Pipeline:
    """L1-penalized logistic regression, which drives noise-SNP coefficients toward zero
    and is a natural fit for sparse polygenic signal recovery. class_weight='balanced'
    compensates for the case/control imbalance typical of disease cohorts."""
    clf = LogisticRegression(
        penalty="l1" if l1 else "l2",
        solver="liblinear",
        C=0.5,
        max_iter=1000,
        class_weight="balanced",
        random_state=0,
    )
    pipeline = Pipeline([("scaler", StandardScaler()), ("clf", clf)])
    # scikit-learn >=1.8 emits a FutureWarning nudging callers toward l1_ratio instead
    # of penalty; liblinear's penalty= API is still fully supported for this demo.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*'penalty' was deprecated.*", category=FutureWarning)
        warnings.filterwarnings("ignore", message=".*Inconsistent values.*", category=UserWarning)
        pipeline.fit(X_train, y_train)
    return pipeline


def train_random_forest(X_train: pd.DataFrame, y_train: pd.Series, n_estimators: int = 300) -> RandomForestClassifier:
    clf = RandomForestClassifier(
        n_estimators=n_estimators, max_depth=6, class_weight="balanced", random_state=0, n_jobs=-1
    )
    clf.fit(X_train, y_train)
    return clf


def evaluate_model(model, X_test: pd.DataFrame, y_test: pd.Series) -> dict:
    y_pred = model.predict(X_test)
    y_proba = model.predict_proba(X_test)[:, 1]

    return {
        "accuracy": accuracy_score(y_test, y_pred),
        "roc_auc": roc_auc_score(y_test, y_proba),
        "confusion_matrix": confusion_matrix(y_test, y_pred).tolist(),
        "fpr_tpr": roc_curve(y_test, y_proba)[:2],
    }


def feature_importance(model, feature_names: list[str]) -> pd.Series:
    """Extract a ranked feature-importance Series from either a Pipeline wrapping
    a LogisticRegression, or a fitted tree ensemble."""
    if isinstance(model, Pipeline):
        coefs = model.named_steps["clf"].coef_.ravel()
        importance = np.abs(coefs)
    else:
        importance = model.feature_importances_
    return pd.Series(importance, index=feature_names).sort_values(ascending=False)
