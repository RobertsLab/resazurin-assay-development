# Clam data conventions

This folder stores clam assay runs used by the dashboard.

## Required naming for plate reader files

Use this pattern for each timepoint file:

- `plate-<plateid>-T<hours>.txt`
- examples: `plate-A-T0.0.txt`, `plate-A-T2.5.txt`, `plate-B-T1.txt`

The dashboard parser reads:

- plate id from `<plateid>`
- timepoint (hours) from `T<hours>`
- well values from the `Results` table inside each file

Supported plate geometries are discovered from the `Results` rows and columns and work for single, 6-well, 12-well, and 24-well exports.

## Optional layout file

Each experiment directory can include `layout.txt`.

Accepted flexible patterns:

1. Explicit mapping lines (tab/comma/semicolon delimited), for example:

`A1\tcontrol\tambient`

2. Matrix format with numeric headers and row letters, for example:

```
	1	2	3	4
A	ctrl	ctrl	treat	treat
B	ctrl	ctrl	treat	treat
C	blank	blank	blank	blank
```

If `layout.txt` cannot be parsed, the dashboard still updates and uses well IDs only.

## Auto-update behavior

On push to `main` that changes `data/clam/**/*.txt`, GitHub Actions rebuilds dashboard data and redeploys the Quarto site.
