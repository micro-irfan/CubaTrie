
## `ZA:Z` SAM Tag Format

`ZA:Z` is emitted when both `--sam` and `--anchors` are enabled:
- on mapped SAM records in anchor mode
- on unmapped (`FLAG 4`) SAM records if anchors were detected but short-reference mapping failed
- on unmapped (`FLAG 4`) SAM records in two-sided mode when full pairing fails but a single start anchor (`a5` or `a3rc`) is confidently detected (`partial=1`)

Expected format (full paired-anchor detection):

```text
ZA:Z:ori=<FWD|RC>;ins=<insert_start_1based>,<insert_len>;ast=5:<pass|fail|na>,3:<pass|fail|na>[;a5=<start_1based>-<end_1based>,ed=<int>,md=<MD_like>][;a3=<start_1based>-<end_1based>,ed=<int>,md=<MD_like>]
```

Expected format (two-sided mode partial start-anchor fallback on unmapped reads):

```text
ZA:Z:ori=<FWD|RC>;partial=1;ast=5:<pass|fail|na>,3:<pass|fail|na>;a5=<start_1based>-<end_1based>,ed=<int>,md=<MD_like>
ZA:Z:ori=<FWD|RC>;partial=1;ast=5:<pass|fail|na>,3:<pass|fail|na>;a3rc=<start_1based>-<end_1based>,ed=<int>,md=<MD_like>
ZA:Z:ori=<FWD|RC>;partial=1;ast=5:<pass|fail|na>,3:<pass|fail|na>;a5=<...>[;a3f=<start_1based>-<end_1based>,len=<int>,ed=<int>,md=<MD_like>]
ZA:Z:ori=<FWD|RC>;partial=1;ast=5:<pass|fail|na>,3:<pass|fail|na>;a3rc=<...>[;a5f=<start_1based>-<end_1based>,len=<int>,ed=<int>,md=<MD_like>]
```

Expected format (two-sided mode unmapped diagnostics):

```text
ZA:Z:ori=<FWD|RC>;ins=<insert_start_1based>,<observed_insert_len>;exp=<reference_len>;reason=<insert_len_lt|insert_len_gt|anchor_ambiguous|anchor_window_rejected>;ast=5:<pass|fail|na>,3:<pass|fail|na>
ZA:Z:ori=<FWD|RC>;ins=<insert_start_1based>,<observed_insert_len>;exp=<reference_len>;reason=<...>;ast=5:<pass|fail|na>,3:<pass|fail|na>;a5=<start_1based>-<end_1based>,ed=<int>,md=<MD_like>;a3=<start_1based>-<end_1based>,ed=<int>,md=<MD_like>
ZA:Z:reason=anchor_not_found;ast=5:<pass|fail|na>,3:<pass|fail|na>;a5f=<start_1based>-<end_1based>,len=<int>,ed=<int>,md=<MD_like>[;a3f=<start_1based>-<end_1based>,len=<int>,ed=<int>,md=<MD_like>]
```

Field meaning:

- `ori`: orientation of anchor detection in the read (`FWD` or `RC`).
- `ins`: extracted insert window used for mapping.
- `a5`: 5' anchor match window in read coordinates.
- `a3`: 3' anchor match window in read coordinates.
- `a3rc`: reverse-complement 3' anchor start hit (reported only in two-sided fallback mode with `ori=RC`).
- `partial=1`: indicates two-sided anchor mode fallback where only a start anchor (`a5` or `a3rc`) was confidently detected.
- `exp`: expected insert/reference length used by `count` mode.
- `reason`: diagnostic reason for unmapped anchor attempts.
- `ast`: anchor status summary where `5` and `3` are each `pass`, `fail`, or `na`.
- `na` means that anchor was not searched in that diagnostic path (or not configured).
- `a5f` / `a3f`: best-effort failed-anchor diagnostics with position range, matched segment length, edit distance, and MD-like string.
- `ed`: anchor edit distance used by anchor matching.
- `md`: MD-like string comparing read anchor segment against expected anchor sequence (`A/C/G/T`, digits for match runs, `^` for deletions from read relative to anchor).

Notes:

- Coordinates in `ins`, `a5`, and `a3` are 1-based.
- Coordinates in `a3rc` are also 1-based.
- `a5` and/or `a3` are present depending on anchor mode (both-sided or one-sided).
- In `count` mode, two-sided anchors (`a5...a3`) must bracket an insert whose length matches the reference length.
- In `count` mode, one-sided anchors keep fixed-length behavior (insert length follows reference length).
- In two-sided mode, if full pairing fails but a single start anchor is confidently found, unmapped SAM may include `partial=1` with either `a5` or `a3rc`.
- Unmapped SAM in anchor mode includes a diagnostic `reason` when extraction fails.
- `reason=anchor_not_found` indicates no confident anchor window was recoverable; when available, failed anchors include `a5f`/`a3f` metadata for debugging.

Examples:

```text
ZA:Z:ori=FWD;ins=31,20;ast=5:pass,3:pass;a5=13-30,ed=0,md=18;a3=51-68,ed=1,md=7A10
ZA:Z:ori=RC;ins=31,20;ast=5:pass,3:pass;a5=55-72,ed=1,md=7A10;a3=9-26,ed=0,md=18
ZA:Z:ori=FWD;ins=19,20;ast=5:pass,3:na;a5=1-18,ed=0,md=18
ZA:Z:ori=FWD;ins=31,22;exp=20;reason=insert_len_gt;ast=5:pass,3:pass;a5=13-30,ed=0,md=18;a3=53-70,ed=1,md=7A10
ZA:Z:ori=FWD;partial=1;ast=5:pass,3:fail;a5=1-18,ed=0,md=18;a3f=45-62,len=18,ed=3,md=4T5A7
ZA:Z:ori=RC;partial=1;ast=5:fail,3:pass;a3rc=3-20,ed=1,md=7A10;a5f=51-68,len=18,ed=2,md=5C12
ZA:Z:reason=anchor_not_found;ast=5:fail,3:na;a5f=2-5,len=4,ed=2,md=0GG0
```