"""Pure cohort-decision helpers for GSE143157 Step 2B."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd


EXPECTED_DONORS = 78
EXPECTED_CONDITIONS = {"Control": 37, "MSA": 41}
PRIMARY_EXCLUDED = {"MSA 30"}
NO_TWO_EXCLUDED = {"MSA 30", "Control 11"}


def as_bool(value: object) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if pd.isna(value):
        return False
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def _require_columns(frame: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def validate_inputs(donors: pd.DataFrame, outliers: pd.DataFrame) -> None:
    _require_columns(donors, ["donor_index", "donor_id", "condition"], "donor metadata")
    _require_columns(
        outliers,
        [
            "donor_id",
            "evidence_count",
            "outlier_evidence",
            "high_priority_review",
            "5mC_flag_pca",
            "5hmC_flag_pca",
            "flag_5hmC_missing_qc",
            "flag_5hmC_negative_qc",
            "flag_detection_qc",
        ],
        "Step 2 outlier review",
    )
    if len(donors) != EXPECTED_DONORS:
        raise ValueError(f"Expected {EXPECTED_DONORS} donors, found {len(donors)}.")
    if donors["donor_id"].duplicated().any():
        raise ValueError("Donor metadata contains duplicate donor IDs.")
    if outliers["donor_id"].duplicated().any():
        raise ValueError("Step 2 outlier review contains duplicate donor IDs.")
    if set(donors["donor_id"].astype(str)) != set(outliers["donor_id"].astype(str)):
        raise ValueError("Step 1 and Step 2 donor sets do not match.")
    condition_counts = donors["condition"].value_counts().to_dict()
    if condition_counts != EXPECTED_CONDITIONS:
        raise ValueError(
            f"Expected condition counts {EXPECTED_CONDITIONS}, found {condition_counts}."
        )
    ordered = pd.to_numeric(donors["donor_index"], errors="raise").astype(int).to_numpy()
    if not np.array_equal(ordered, np.arange(EXPECTED_DONORS)):
        raise ValueError("donor_index must be the aligned HDF5 order 0..77.")

    indexed = outliers.set_index("donor_id")
    msa30 = indexed.loc["MSA 30"]
    required_msa30 = [
        "5hmC_flag_pca",
        "flag_5hmC_missing_qc",
        "flag_5hmC_negative_qc",
        "flag_detection_qc",
    ]
    if not all(as_bool(msa30[column]) for column in required_msa30):
        raise ValueError(
            "MSA 30 does not have the four confirmed Step 2 technical/QC signals. "
            "The frozen cohort rule must be reviewed instead of silently changing it."
        )

    control11 = indexed.loc["Control 11"]
    if not (as_bool(control11["5mC_flag_pca"]) and as_bool(control11["5hmC_flag_pca"])):
        raise ValueError("Control 11 no longer has the confirmed cross-modality PCA evidence.")
    if any(
        as_bool(control11[column])
        for column in ["flag_5hmC_missing_qc", "flag_5hmC_negative_qc", "flag_detection_qc"]
    ):
        raise ValueError(
            "Control 11 now has an objective technical QC flag; the primary cohort rule needs review."
        )


def build_manifest(donors: pd.DataFrame, outliers: pd.DataFrame) -> pd.DataFrame:
    validate_inputs(donors, outliers)
    evidence = outliers[
        [
            "donor_id",
            "evidence_count",
            "outlier_evidence",
            "high_priority_review",
        ]
    ].copy()
    evidence = evidence.rename(
        columns={
            "evidence_count": "step2_evidence_count",
            "outlier_evidence": "step2_outlier_evidence",
            "high_priority_review": "step2_high_priority_review",
        }
    )
    manifest = donors[["donor_index", "donor_id", "condition"]].copy()
    manifest = manifest.rename(columns={"donor_index": "original_hdf5_column"})
    manifest = manifest.merge(evidence, on="donor_id", how="left", validate="one_to_one")

    manifest["include_sensitivity_full78"] = True
    manifest["sensitivity_full78_column"] = np.arange(len(manifest), dtype=int)

    manifest["include_primary_clean77"] = ~manifest["donor_id"].isin(PRIMARY_EXCLUDED)
    primary_index = pd.Series(pd.NA, index=manifest.index, dtype="Int64")
    primary_index.loc[manifest["include_primary_clean77"]] = np.arange(
        int(manifest["include_primary_clean77"].sum()), dtype=int
    )
    manifest["primary_clean77_column"] = primary_index
    manifest["primary_decision"] = np.where(
        manifest["donor_id"].eq("MSA 30"),
        "Excluded from paired primary cohort: BS/5hmC technical QC failure",
        "Included in paired primary cohort",
    )

    manifest["include_sensitivity_no2_76"] = ~manifest["donor_id"].isin(NO_TWO_EXCLUDED)
    no_two_index = pd.Series(pd.NA, index=manifest.index, dtype="Int64")
    no_two_index.loc[manifest["include_sensitivity_no2_76"]] = np.arange(
        int(manifest["include_sensitivity_no2_76"].sum()), dtype=int
    )
    manifest["sensitivity_no2_76_column"] = no_two_index
    manifest["sensitivity_no2_note"] = np.select(
        [
            manifest["donor_id"].eq("MSA 30"),
            manifest["donor_id"].eq("Control 11"),
        ],
        [
            "Removed for objective technical QC failure",
            "Removed only in stress-test because of cross-modality PCA evidence",
        ],
        default="Included",
    )
    return manifest


def subset_donor_metadata(
    donors: pd.DataFrame, manifest: pd.DataFrame, include_column: str, cohort_column: str
) -> pd.DataFrame:
    mapping = manifest.loc[
        manifest[include_column], ["donor_id", "original_hdf5_column", cohort_column]
    ].copy()
    result = mapping.merge(donors, on="donor_id", how="left", validate="one_to_one")
    result = result.drop(columns=["donor_index"]).rename(
        columns={cohort_column: "donor_index"}
    )
    order = ["donor_index", "original_hdf5_column", "donor_id"]
    return result[order + [column for column in result.columns if column not in order]].sort_values(
        "donor_index", kind="stable"
    )


def subset_array_metadata(
    arrays: pd.DataFrame, donor_metadata: pd.DataFrame
) -> pd.DataFrame:
    _require_columns(arrays, ["donor_id"], "array metadata")
    mapping = donor_metadata[["donor_id", "donor_index", "original_hdf5_column"]]
    result = arrays.loc[arrays["donor_id"].isin(mapping["donor_id"])].copy()
    result = result.merge(mapping, on="donor_id", how="left", validate="many_to_one")
    front = ["donor_index", "original_hdf5_column", "donor_id"]
    return result[front + [column for column in result.columns if column not in front]].sort_values(
        ["donor_index", "matrix_beta_index"], kind="stable"
    )


def cohort_counts(metadata: pd.DataFrame) -> dict[str, object]:
    age = pd.to_numeric(metadata.get("age"), errors="coerce")
    return {
        "donors": int(len(metadata)),
        "condition_counts": {
            str(key): int(value) for key, value in metadata["condition"].value_counts().items()
        },
        "sex_counts": {
            str(key): int(value) for key, value in metadata.get("sex", pd.Series(dtype=str)).value_counts().items()
        },
        "age_mean": None if age.notna().sum() == 0 else float(age.mean()),
        "age_min": None if age.notna().sum() == 0 else float(age.min()),
        "age_max": None if age.notna().sum() == 0 else float(age.max()),
    }
