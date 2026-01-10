#!/usr/bin/env python
"""
Train pathogenicity classifier on variant features.

Usage:
    python training/train_classifier.py --panel hereditary_cancer
    python training/train_classifier.py --input data/exports/custom_classifier_features.jsonl
"""

import json
import argparse
from pathlib import Path

import numpy as np


def load_data(filepath: Path):
    """Load and prepare training data."""
    X, y, ids = [], [], []

    feature_cols = [
        "sift_score", "polyphen2_hdiv", "polyphen2_hvar",
        "cadd_phred", "revel", "alphamissense", "clinpred",
        "dann", "phylop100", "phastcons100", "gerp", "gnomad_af"
    ]

    with open(filepath) as f:
        for line in f:
            row = json.loads(line)

            # Skip if no ClinVar label
            clinvar = row.get("clinvar_sig", "")
            if not clinvar:
                continue

            # Map ClinVar to binary label
            clinvar_lower = clinvar.lower()
            if "pathogenic" in clinvar_lower and "conflict" not in clinvar_lower:
                label = 1
            elif "benign" in clinvar_lower and "conflict" not in clinvar_lower:
                label = 0
            else:
                continue  # Skip uncertain/conflicting

            # Extract features
            features = []
            for col in feature_cols:
                val = row.get(col)
                if val is None or val == "" or val == ".":
                    features.append(np.nan)
                else:
                    try:
                        features.append(float(val))
                    except ValueError:
                        features.append(np.nan)

            X.append(features)
            y.append(label)
            ids.append(row["id"])

    return np.array(X), np.array(y), ids, feature_cols


def train_xgboost(X, y, feature_names):
    """Train XGBoost classifier."""
    try:
        import xgboost as xgb
        from sklearn.model_selection import train_test_split
        from sklearn.metrics import classification_report, roc_auc_score
    except ImportError:
        print("Install required packages: pip install xgboost scikit-learn")
        return None

    # Handle missing values
    from sklearn.impute import SimpleImputer
    imputer = SimpleImputer(strategy='median')
    X_imputed = imputer.fit_transform(X)

    # Split
    X_train, X_test, y_train, y_test = train_test_split(
        X_imputed, y, test_size=0.2, random_state=42, stratify=y
    )

    # Train
    model = xgb.XGBClassifier(
        n_estimators=100,
        max_depth=6,
        learning_rate=0.1,
        random_state=42,
        use_label_encoder=False,
        eval_metric='logloss'
    )

    model.fit(X_train, y_train)

    # Evaluate
    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)[:, 1]

    print("\n=== Classification Report ===")
    print(classification_report(y_test, y_pred, target_names=["Benign", "Pathogenic"]))

    print(f"ROC-AUC: {roc_auc_score(y_test, y_prob):.4f}")

    # Feature importance
    print("\n=== Feature Importance ===")
    importance = sorted(zip(feature_names, model.feature_importances_),
                       key=lambda x: x[1], reverse=True)
    for feat, imp in importance:
        print(f"  {feat}: {imp:.4f}")

    return model, imputer


def main():
    parser = argparse.ArgumentParser(description="Train pathogenicity classifier")
    parser.add_argument("--panel", "-p", help="Panel name to export and train on")
    parser.add_argument("--input", "-i", help="Pre-exported JSONL file")
    parser.add_argument("--output", "-o", default="models/classifier.json", help="Output model path")
    args = parser.parse_args()

    # Get training data
    if args.input:
        data_path = Path(args.input)
    elif args.panel:
        from src.export import export_for_classifier
        data_path = export_for_classifier(args.panel)
    else:
        print("Specify --panel or --input")
        return

    print(f"Loading data from {data_path}...")
    X, y, ids, feature_names = load_data(data_path)

    print(f"Loaded {len(y)} labeled variants")
    print(f"  Pathogenic: {sum(y)}")
    print(f"  Benign: {len(y) - sum(y)}")

    if len(y) < 100:
        print("Warning: Very small dataset. Results may not be reliable.")

    # Train
    result = train_xgboost(X, y, feature_names)

    if result:
        model, imputer = result
        # Save model
        output_dir = Path(args.output).parent
        output_dir.mkdir(parents=True, exist_ok=True)
        model.save_model(args.output.replace('.json', '.xgb'))
        print(f"\nModel saved to {args.output.replace('.json', '.xgb')}")


if __name__ == "__main__":
    main()
