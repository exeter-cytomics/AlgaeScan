# ==============================
# Pipeline: Supervised & Unsupervised Models
# ==============================

from __future__ import annotations

import itertools
import random
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import IsolationForest, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
)
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import OneClassSVM


RANDOM_STATE = 123
ISOLATION_FOREST_RANDOM_STATE = 42
ALGAE_LABEL = "algae"
NON_ALGAE_LABEL = "non-algae"

MODEL_NAMES = [
    "A_rf",
    "A_if",
    "A_ocsvm",
    "B_rf",
    "B_mlr",
]


def set_seed(seed=RANDOM_STATE):
    """Set Python and NumPy seeds.

    Parameters:
        seed (int): Reproducibility seed. Default is 123.

    Returns:
        None.
    """
    random.seed(seed)
    np.random.seed(seed)


set_seed()


# ==============================
# Data Helpers
# ==============================

def load_dataset(dataset_path):
    """Load a CSV or TSV dataset.

    Parameters:
        dataset_path (str or Path): Input file path.

    Returns:
        pandas.DataFrame: Loaded data.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file type is unsupported or the data is empty.
    """
    dataset_path = Path(dataset_path)
    if not dataset_path.exists():
        raise FileNotFoundError(f"Dataset not found: {dataset_path}")

    if dataset_path.suffix.lower() == ".csv":
        df = pd.read_csv(dataset_path)
    elif dataset_path.suffix.lower() in [".tsv", ".tab"]:
        df = pd.read_csv(dataset_path, sep="\t")
    else:
        raise ValueError("Only CSV and TSV files are supported.")

    if df.empty:
        raise ValueError("Dataset is empty.")
    return df


def save_dataset(df, output_path):
    """Save a DataFrame as CSV or TSV.

    Parameters:
        df (pandas.DataFrame): Data to save.
        output_path (str or Path): Output file path.

    Returns:
        Path: Saved file path.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.suffix.lower() == ".csv":
        df.to_csv(output_path, index=False)
    elif output_path.suffix.lower() in [".tsv", ".tab"]:
        df.to_csv(output_path, sep="\t", index=False)
    else:
        raise ValueError("Only CSV and TSV files are supported.")

    return output_path


def get_label_column(df, label_column=None):
    """Return the selected label column.

    Parameters:
        df (pandas.DataFrame): Dataset.
        label_column (str or None): Label column name. If None, the last column is used.

    Returns:
        str: Valid label column name.
    """
    if df.empty:
        raise ValueError("Dataset is empty.")

    label_column = df.columns[-1] if label_column is None else label_column
    if label_column not in df.columns:
        raise ValueError(f"Label column not found: {label_column}")
    return label_column


def put_label_last(df, label_column):
    """Move the selected label column to the end.

    Parameters:
        df (pandas.DataFrame): Dataset.
        label_column (str): Label column name.

    Returns:
        pandas.DataFrame: Dataset with label column last.
    """
    label_column = get_label_column(df, label_column)
    columns = [col for col in df.columns if col != label_column] + [label_column]
    return df.loc[:, columns]


def create_small_dataset(
    input_path,
    output_path,
    n_samples,
    label_column=None,
    sampling_strategy="stratified",
    random_state=RANDOM_STATE,
):
    """Create a smaller dataset from an existing dataset.

    Parameters:
        input_path (str or Path): Source dataset path.
        output_path (str or Path): Saved small dataset path.
        n_samples (int): Number of rows requested.
        label_column (str or None): Column used for stratified sampling.
        sampling_strategy (str): "stratified" or "random".
        random_state (int): Sampling seed. Default is 123.

    Returns:
        pandas.DataFrame: Sampled dataset. If n_samples is greater than the source
        size, a shuffled full copy is returned.
    """
    df = load_dataset(input_path)
    label_column = get_label_column(df, label_column)

    if not isinstance(n_samples, int) or n_samples <= 0:
        raise ValueError("n_samples must be a positive integer.")
    if sampling_strategy not in ["stratified", "random"]:
        raise ValueError("sampling_strategy must be 'stratified' or 'random'.")

    if n_samples >= len(df):
        sampled = df.sample(frac=1, random_state=random_state).reset_index(drop=True)
    elif sampling_strategy == "random":
        sampled = df.sample(n=n_samples, random_state=random_state).reset_index(drop=True)
    else:
        sampled = _stratified_sample(df, label_column, n_samples, random_state)

    sampled = put_label_last(sampled, label_column)
    save_dataset(sampled, output_path)
    return sampled


def _stratified_sample(df, label_column, n_samples, random_state):
    labels = df[label_column]
    class_counts = labels.value_counts()

    if n_samples < len(class_counts):
        warnings.warn(
            "n_samples is smaller than the number of classes; not every class can appear.",
            UserWarning,
        )

    sample_parts = []
    remaining = n_samples
    for i, (label, count) in enumerate(class_counts.items()):
        classes_left = len(class_counts) - i
        if classes_left == 1:
            class_n = min(count, remaining)
        else:
            class_n = round(n_samples * count / len(df))
            class_n = max(1, class_n)
            class_n = min(class_n, count, remaining - classes_left + 1)

        sample_parts.append(
            df[df[label_column] == label].sample(n=class_n, random_state=random_state)
        )
        remaining -= class_n

    sampled = pd.concat(sample_parts)
    return sampled.sample(frac=1, random_state=random_state).reset_index(drop=True)


def remove_characters_col(df, target_col=None, drop_columns=None):
    """Remove label and metadata columns before model fitting.

    Parameters:
        df (pandas.DataFrame): Input data.
        target_col (str or None): Label column to remove.
        drop_columns (list or None): Extra columns to remove.

    Returns:
        pandas.DataFrame: Numeric feature data.
    """
    drop_columns = [] if drop_columns is None else list(drop_columns)
    columns_to_drop = [col for col in drop_columns if col in df.columns]
    if target_col is not None and target_col in df.columns:
        columns_to_drop.append(target_col)

    X = df.drop(columns=list(dict.fromkeys(columns_to_drop)), errors="ignore")
    non_numeric = [
        col for col in X.columns if not pd.api.types.is_numeric_dtype(X[col])
    ]
    if non_numeric:
        raise ValueError(
            "Model features must be numeric. Remove these columns: "
            + ", ".join(non_numeric)
        )
    return X


def make_binary_label(y, non_algae_label=NON_ALGAE_LABEL):
    """Create algae/non-algae labels for Model A.

    Parameters:
        y (pandas.Series): Original label values.
        non_algae_label (str): Label treated as non-algae.

    Returns:
        pandas.Series: Binary algae/non-algae target.
    """
    y_text = y.astype(str)
    if non_algae_label not in set(y_text):
        raise ValueError(f"Non-algae label not found: {non_algae_label}")
    return pd.Series(
        np.where(y_text == non_algae_label, non_algae_label, ALGAE_LABEL),
        index=y.index,
        name=y.name,
    )


def filter_algae_data(df, species_col="Species", non_algae_label=NON_ALGAE_LABEL):
    """Keep only algae rows for Model B.

    Parameters:
        df (pandas.DataFrame): Dataset.
        species_col (str): Species label column.
        non_algae_label (str): Non-algae label value.

    Returns:
        pandas.DataFrame: Algae-only data.
    """
    if species_col not in df.columns:
        raise ValueError(f"Species column not found: {species_col}")

    algae_df = df[df[species_col].astype(str) != non_algae_label].copy()
    if algae_df[species_col].nunique() < 2:
        raise ValueError("Model B needs at least two algae species.")
    return algae_df


def create_preprocessor(scale=False):
    """Create simple numeric preprocessing.

    Parameters:
        scale (bool): Add StandardScaler when True.

    Returns:
        sklearn.pipeline.Pipeline: Imputer, and optionally scaler.
    """
    steps = [("imputer", SimpleImputer(strategy="median"))]
    if scale:
        steps.append(("scaler", StandardScaler()))
    return Pipeline(steps)


def _cv_object(y, cv):
    counts = pd.Series(y).value_counts()
    if counts.min() < 2:
        raise ValueError("Each class needs at least two rows for cross validation.")

    n_splits = min(cv, int(counts.min()))
    return StratifiedKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=RANDOM_STATE,
    )


def _grid_cv_model(pipe, param_grid, X, y, cv=5, scoring="accuracy"):
    cv_obj = _cv_object(y, cv)
    cv_model = GridSearchCV(
        estimator=pipe,
        param_grid=param_grid,
        cv=cv_obj,
        scoring=scoring,
    )
    cv_model.fit(X, y)
    model = cv_model.best_estimator_
    model.best_params_ = cv_model.best_params_
    model.best_cv_score_ = cv_model.best_score_
    return model


def _limit_mtry(mtry, n_features, model_name):
    if n_features < mtry:
        warnings.warn(
            f"{model_name}: mtry={mtry} was limited to {n_features} features.",
            UserWarning,
        )
        return n_features
    return mtry


def _save_feature_order(model, feature_columns, model_kind, non_algae_label=None):
    model.feature_columns_ = list(feature_columns)
    model.model_kind_ = model_kind
    model.non_algae_label_ = non_algae_label
    return model


def _align_features(test_data, model):
    missing = [col for col in model.feature_columns_ if col not in test_data.columns]
    if missing:
        raise ValueError("Missing feature columns: " + ", ".join(missing))
    return test_data.loc[:, model.feature_columns_]


def downsample(df, target_col, n=300, random_state=RANDOM_STATE):
    """Create a small balanced sample.

    Parameters:
        df (pandas.DataFrame): Source data.
        target_col (str): Label column.
        n (int): Requested row count.
        random_state (int): Sampling seed.

    Returns:
        tuple: X and y from the sampled data.
    """
    target_col = get_label_column(df, target_col)
    sampled = _stratified_sample(df, target_col, n, random_state)
    X = sampled.drop(columns=[target_col])
    y = sampled[target_col]
    return X, y


# ==============================
# Model A: algae vs non-algae
# ==============================

def Train_model_A_rf(
    train_data,
    downsampling=False,
    var_downs="type",
    drop_columns=None,
    cv=5,
    scoring="accuracy",
):
    """Train Model A Random Forest with cross validation.

    Parameters:
        train_data (pandas.DataFrame): Dataset with feature columns and labels.
        downsampling (bool): Apply balanced downsampling before training.
        var_downs (str): Binary label column, usually type.
        drop_columns (list or None): Metadata columns removed from features.
        cv (int): Number of cross-validation folds.
        scoring (str): GridSearchCV scoring metric.

    Returns:
        sklearn.pipeline.Pipeline: Best cross-validated model.
    """
    data = train_data.copy()
    if downsampling:
        X_sample, y_sample = downsample(data, var_downs)
        data = pd.concat([X_sample, y_sample], axis=1)

    X_train = remove_characters_col(data, target_col=var_downs, drop_columns=drop_columns)
    y_train = make_binary_label(data[var_downs])
    max_features = _limit_mtry(20, X_train.shape[1], "Model A Random Forest")

    pipe = Pipeline([
        ("preprocess", create_preprocessor(scale=False)),
        ("model", RandomForestClassifier(random_state=RANDOM_STATE)),
    ])
    param_grid = {
        "model__n_estimators": [15],
        "model__max_features": [max_features],
    }

    model = _grid_cv_model(pipe, param_grid, X_train, y_train, cv=cv, scoring=scoring)
    return _save_feature_order(model, X_train.columns, "supervised")


def Predict_model_A_rf(test_data, model):
    """Predict algae/non-algae with Model A Random Forest.

    Parameters:
        test_data (pandas.DataFrame): Data to predict.
        model (sklearn.pipeline.Pipeline): Trained model.

    Returns:
        pandas.DataFrame: Predicted labels and probabilities.
    """
    X_test = _align_features(test_data, model)
    result = pd.DataFrame({"predicted_label": model.predict(X_test)})
    probabilities = model.predict_proba(X_test)
    for label, values in zip(model.named_steps["model"].classes_, probabilities.T):
        result[f"probability_{label}"] = values
    return result


def Train_model_A_if(
    train_data,
    downsampling=False,
    var_downs="type",
    drop_columns=None,
    cv=5,
    scoring="balanced_accuracy",
):
    """Train Model A Isolation Forest with cross validation.

    Parameters:
        train_data (pandas.DataFrame): Dataset with algae and non-algae rows.
        downsampling (bool): Apply balanced downsampling before training.
        var_downs (str): Label column used to find algae rows.
        drop_columns (list or None): Metadata columns removed from features.
        cv (int): Number of cross-validation folds.
        scoring (str): Validation metric name.

    Returns:
        sklearn.pipeline.Pipeline: Best cross-validated model trained on algae rows.
    """
    data = train_data.copy()
    if downsampling:
        X_sample, y_sample = downsample(data, var_downs)
        data = pd.concat([X_sample, y_sample], axis=1)

    X_train = remove_characters_col(data, target_col=var_downs, drop_columns=drop_columns)
    y_train = make_binary_label(data[var_downs])
    param_grid = {
        "n_estimators": [100],
        "contamination": [0.1],
    }

    model = _cv_one_class_model(
        X_train,
        y_train,
        estimator_name="isolation_forest",
        param_grid=param_grid,
        cv=cv,
        scoring=scoring,
    )
    return _save_feature_order(model, X_train.columns, "one_class", NON_ALGAE_LABEL)


def Predict_model_A_if(test_data, model):
    """Predict algae/non-algae with Isolation Forest.

    Parameters:
        test_data (pandas.DataFrame): Data to predict.
        model (sklearn.pipeline.Pipeline): Trained model.

    Returns:
        pandas.DataFrame: Predicted labels, raw values, and scores.
    """
    return _predict_one_class(test_data, model)


def Train_model_A_ocsvm(
    train_data,
    downsampling=False,
    var_downs="type",
    drop_columns=None,
    cv=5,
    scoring="balanced_accuracy",
):
    """Train Model A One-Class SVM with cross validation.

    Parameters:
        train_data (pandas.DataFrame): Dataset with algae and non-algae rows.
        downsampling (bool): Apply balanced downsampling before training.
        var_downs (str): Label column used to find algae rows.
        drop_columns (list or None): Metadata columns removed from features.
        cv (int): Number of cross-validation folds.
        scoring (str): Validation metric name.

    Returns:
        sklearn.pipeline.Pipeline: Best cross-validated model trained on algae rows.
    """
    data = train_data.copy()
    if downsampling:
        X_sample, y_sample = downsample(data, var_downs)
        data = pd.concat([X_sample, y_sample], axis=1)

    X_train = remove_characters_col(data, target_col=var_downs, drop_columns=drop_columns)
    y_train = make_binary_label(data[var_downs])
    param_grid = {
        "kernel": ["rbf"],
        "nu": [0.1],
    }

    model = _cv_one_class_model(
        X_train,
        y_train,
        estimator_name="one_class_svm",
        param_grid=param_grid,
        cv=cv,
        scoring=scoring,
    )
    return _save_feature_order(model, X_train.columns, "one_class", NON_ALGAE_LABEL)


def Predict_model_A_ocsvm(test_data, model):
    """Predict algae/non-algae with One-Class SVM.

    Parameters:
        test_data (pandas.DataFrame): Data to predict.
        model (sklearn.pipeline.Pipeline): Trained model.

    Returns:
        pandas.DataFrame: Predicted labels, raw values, and scores.
    """
    return _predict_one_class(test_data, model)


def _cv_one_class_model(X, y_binary, estimator_name, param_grid, cv=5, scoring="balanced_accuracy"):
    cv_obj = _cv_object(y_binary, cv)
    keys = list(param_grid.keys())
    param_sets = [dict(zip(keys, values)) for values in itertools.product(*param_grid.values())]
    best_score = -np.inf
    best_params = None

    for params in param_sets:
        fold_scores = []
        for train_index, test_index in cv_obj.split(X, y_binary):
            X_fold_train = X.iloc[train_index]
            y_fold_train = y_binary.iloc[train_index]
            X_fold_test = X.iloc[test_index]
            y_fold_test = y_binary.iloc[test_index]

            algae_X = X_fold_train[y_fold_train == ALGAE_LABEL]
            if algae_X.empty:
                raise ValueError("No algae samples available for one-class training.")

            model = _build_one_class_pipe(estimator_name, params)
            model.fit(algae_X)
            predicted = _map_one_class_predictions(model.predict(X_fold_test))
            fold_scores.append(_score_labels(y_fold_test, predicted, scoring))

        mean_score = float(np.mean(fold_scores))
        if mean_score > best_score:
            best_score = mean_score
            best_params = params

    algae_X = X[y_binary == ALGAE_LABEL]
    final_model = _build_one_class_pipe(estimator_name, best_params)
    final_model.fit(algae_X)
    final_model.best_params_ = best_params
    final_model.best_cv_score_ = best_score
    final_model.training_rows_ = len(algae_X)
    return final_model


def _build_one_class_pipe(estimator_name, params):
    if estimator_name == "isolation_forest":
        estimator = IsolationForest(
            n_estimators=params["n_estimators"],
            contamination=params["contamination"],
            random_state=ISOLATION_FOREST_RANDOM_STATE,
        )
        scale = False
    elif estimator_name == "one_class_svm":
        estimator = OneClassSVM(kernel=params["kernel"], nu=params["nu"])
        scale = True
    else:
        raise ValueError(f"Unknown one-class model: {estimator_name}")

    return Pipeline([
        ("preprocess", create_preprocessor(scale=scale)),
        ("model", estimator),
    ])


def _map_one_class_predictions(raw_predictions):
    return pd.Series(
        np.where(raw_predictions == 1, ALGAE_LABEL, NON_ALGAE_LABEL)
    )


def _predict_one_class(test_data, model):
    X_test = _align_features(test_data, model)
    raw_predictions = model.predict(X_test)
    result = pd.DataFrame({
        "predicted_label": _map_one_class_predictions(raw_predictions),
        "raw_prediction": raw_predictions,
    })
    result["decision_score"] = model.decision_function(X_test)
    if hasattr(model.named_steps["model"], "score_samples"):
        result["anomaly_score"] = -model.score_samples(X_test)
    return result


def _score_labels(y_true, y_pred, scoring):
    if scoring == "balanced_accuracy":
        return balanced_accuracy_score(y_true, y_pred)
    if scoring == "accuracy":
        return accuracy_score(y_true, y_pred)
    if scoring == "macro_f1":
        return f1_score(y_true, y_pred, average="macro")
    raise ValueError("Unsupported scoring value.")


# ==============================
# Model B: algae species
# ==============================

def train_model_B_rf(
    train_data,
    downsampling=False,
    var_downs="Species",
    drop_columns=None,
    non_algae_label=NON_ALGAE_LABEL,
    cv=5,
    scoring="accuracy",
):
    """Train Model B Random Forest with cross validation.

    Parameters:
        train_data (pandas.DataFrame): Dataset with species labels.
        downsampling (bool): Apply balanced downsampling before training.
        var_downs (str): Species label column.
        drop_columns (list or None): Metadata columns removed from features.
        non_algae_label (str): Non-algae value removed before training.
        cv (int): Number of cross-validation folds.
        scoring (str): GridSearchCV scoring metric.

    Returns:
        sklearn.pipeline.Pipeline: Best cross-validated species model.
    """
    data = filter_algae_data(train_data, var_downs, non_algae_label)
    if downsampling:
        X_sample, y_sample = downsample(data, var_downs)
        data = pd.concat([X_sample, y_sample], axis=1)

    X_train = remove_characters_col(data, target_col=var_downs, drop_columns=drop_columns)
    y_train = data[var_downs]
    max_features = _limit_mtry(40, X_train.shape[1], "Model B Random Forest")

    pipe = Pipeline([
        ("preprocess", create_preprocessor(scale=False)),
        ("model", RandomForestClassifier(random_state=RANDOM_STATE)),
    ])
    param_grid = {
        "model__n_estimators": [15],
        "model__max_features": [max_features],
    }

    model = _grid_cv_model(pipe, param_grid, X_train, y_train, cv=cv, scoring=scoring)
    model.training_rows_ = len(data)
    return _save_feature_order(model, X_train.columns, "supervised")


def Predict_model_B_rf(test_data, model):
    """Predict algae species with Model B Random Forest.

    Parameters:
        test_data (pandas.DataFrame): Algae data to predict.
        model (sklearn.pipeline.Pipeline): Trained model.

    Returns:
        pandas.DataFrame: Predicted species and probabilities.
    """
    return _predict_supervised(test_data, model)


def train_model_B_mlr(
    train_data,
    downsampling=False,
    var_downs="Species",
    drop_columns=None,
    non_algae_label=NON_ALGAE_LABEL,
    decay=0.1,
    cv=5,
    scoring="accuracy",
):
    """Train Model B multinomial logistic regression with cross validation.

    Parameters:
        train_data (pandas.DataFrame): Dataset with species labels.
        downsampling (bool): Apply balanced downsampling before training.
        var_downs (str): Species label column.
        drop_columns (list or None): Metadata columns removed from features.
        non_algae_label (str): Non-algae value removed before training.
        decay (float): L2 regularization lambda. C is calculated as 1 / decay.
        cv (int): Number of cross-validation folds.
        scoring (str): GridSearchCV scoring metric.

    Returns:
        sklearn.pipeline.Pipeline: Best cross-validated multinomial model.
    """
    if decay <= 0:
        raise ValueError("decay must be greater than zero.")

    data = filter_algae_data(train_data, var_downs, non_algae_label)
    if downsampling:
        X_sample, y_sample = downsample(data, var_downs)
        data = pd.concat([X_sample, y_sample], axis=1)

    X_train = remove_characters_col(data, target_col=var_downs, drop_columns=drop_columns)
    y_train = data[var_downs]

    pipe = Pipeline([
        ("preprocess", create_preprocessor(scale=True)),
        ("model", LogisticRegression(solver="lbfgs", max_iter=1000, random_state=RANDOM_STATE)),
    ])
    param_grid = {
        "model__C": [1 / decay],
    }

    model = _grid_cv_model(pipe, param_grid, X_train, y_train, cv=cv, scoring=scoring)
    model.training_rows_ = len(data)
    model.decay_ = decay
    return _save_feature_order(model, X_train.columns, "supervised")


def Predict_model_B_mlr(test_data, model):
    """Predict algae species with multinomial logistic regression.

    Parameters:
        test_data (pandas.DataFrame): Algae data to predict.
        model (sklearn.pipeline.Pipeline): Trained model.

    Returns:
        pandas.DataFrame: Predicted species and probabilities.
    """
    return _predict_supervised(test_data, model)


def _predict_supervised(test_data, model):
    X_test = _align_features(test_data, model)
    result = pd.DataFrame({"predicted_label": model.predict(X_test)})
    if hasattr(model, "predict_proba"):
        probabilities = model.predict_proba(X_test)
        for label, values in zip(model.named_steps["model"].classes_, probabilities.T):
            result[f"probability_{label}"] = values
    return result


def prediction_function(test_data, model):
    """Predict with a trained pipeline.

    Parameters:
        test_data (pandas.DataFrame): Data to predict.
        model (sklearn.pipeline.Pipeline): Trained model.

    Returns:
        pandas.DataFrame: Prediction results.
    """
    if getattr(model, "model_kind_", "") == "one_class":
        return _predict_one_class(test_data, model)
    return _predict_supervised(test_data, model)


# ==============================
# Evaluation and Output
# ==============================

def evaluate_predictions(y_true, predictions):
    """Calculate classification metrics.

    Parameters:
        y_true (array-like): True labels.
        predictions (pandas.DataFrame): Prediction results with predicted_label.

    Returns:
        dict: Accuracy, balanced accuracy, macro F1, confusion matrix, and report.
    """
    y_pred = predictions["predicted_label"]
    labels = sorted(set(pd.Series(y_true).astype(str)) | set(pd.Series(y_pred).astype(str)))
    return {
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "macro_f1": f1_score(y_true, y_pred, average="macro"),
        "confusion_matrix": confusion_matrix(y_true, y_pred, labels=labels).tolist(),
        "classification_report": classification_report(
            y_true,
            y_pred,
            labels=labels,
            output_dict=True,
            zero_division=0,
        ),
    }


def save_prediction_results(predictions, y_true, model_name, output_dir):
    """Save predictions as a CSV file.

    Parameters:
        predictions (pandas.DataFrame): Prediction results.
        y_true (array-like): True labels.
        model_name (str): Model name for the file.
        output_dir (str or Path): Output folder.

    Returns:
        Path: Saved CSV path.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_df = predictions.copy()
    output_df.insert(0, "true_label", list(y_true))
    output_df.insert(0, "model_name", model_name)

    output_path = output_dir / f"{model_name}_predictions.csv"
    output_df.to_csv(output_path, index=False)
    return output_path


def save_model(model, output_path):
    """Save a trained model with joblib.

    Parameters:
        model (sklearn.pipeline.Pipeline): Trained model.
        output_path (str or Path): Model file path.

    Returns:
        Path: Saved model path.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(model, output_path)
    return output_path


# =========================================
# Small utility functions for this script
# =========================================

def save_metrics(metrics_rows, output_dir):
    """Save compact metric results."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    csv_path = output_dir / "metrics_summary.csv"
    json_path = output_dir / "metrics_summary.json"

    pd.DataFrame(metrics_rows).to_csv(csv_path, index=False)
    with json_path.open("w", encoding="utf-8") as file:
        json.dump(metrics_rows, file, indent=2)

    return csv_path, json_path


def add_metric_row(model_name, metrics, predictions):
    """Create one short metric row."""
    return {
        "model_name": model_name,
        "accuracy": metrics["accuracy"],
        "balanced_accuracy": metrics["balanced_accuracy"],
        "macro_f1": metrics["macro_f1"],
        "n_predictions": len(predictions),
    }


def print_metric_row(row):
    """Print one compact result line."""
    print(
        f"{row['model_name']}: "
        f"accuracy={row['accuracy']:.4f}, "
        f"balanced_accuracy={row['balanced_accuracy']:.4f}, "
        f"macro_f1={row['macro_f1']:.4f}, "
        f"predictions={row['n_predictions']}"
    )