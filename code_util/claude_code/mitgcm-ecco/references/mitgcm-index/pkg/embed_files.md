# pkg/embed_files

Embeds input files into the executable (e.g. for testing/portability).

**runtime switch:** `useEMBED_FILES`-style flag in `data.pkg` (check exact name in packages_boot.F)

## README
```

 =================
    EMBED_FILES
 =================


This package is quick and portable way to embed files (any general
binary and/or text data) within an executable for later extraction.

```

## Headers
- `EMBED_FILES_OPTIONS.h` — Place CPP define/undef flag here

## Routines (1)
`embed_files_init.F`

## Called from outside the package
- `EMBED_FILES_INIT` ← `model/src/packages_init_fixed.F:637`
