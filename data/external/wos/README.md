# Local Web of Science input

Place the Web of Science Core Collection export used for the optional audit workflow here as:

```text
wos_raw_export.xls
```

The expected search is documented in:

```text
documentation/wos_search.md
```

This raw export is intentionally ignored by Git and is not redistributed in the public repository.

The core analytical pipeline does not require this file because the final source-verified supplementary records are retained separately under:

```text
data/curated/wos/wos_source_verified_records.xlsx
```
