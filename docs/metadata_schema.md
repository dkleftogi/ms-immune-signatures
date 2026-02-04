Required columns (minimum):

Filename — FCS filename (must exist inside input_dir)

batch — batch/run ID used for harmonisation

condition — biological condition if you want to preserve it (optional but recommended)

(optional) Patient_id, TimePoint, Cohort — useful for downstream stratification

Important: Filename values must match the actual .fcs file names exactly (including spaces). If you have trailing spaces in metadata rows, strip them.