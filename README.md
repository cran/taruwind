# TaruWind 0.1.1

## Start here

**[Step-by-step getting-started guide](inst/doc/getting-started.md)** —
one connected Istanbul example, from coordinates and downloads to interpreting
results, extracting the four raw ERA5 nodes, saving/reusing responses, adding
a second site and mapping terrain inputs and a small local area.

The formatted guide is also installed locally:

```r
library(taruwind)
browseURL(system.file("doc", "getting-started.html", package = "taruwind"))
```

The HTML guide has a **Copy code** button on each example; long lines wrap to
fit the window. To regenerate it from the Markdown source while retaining
those controls, run this from the package source directory:

```sh
pandoc inst/doc/getting-started.md --standalone --toc \
  --include-in-header=tools/guide-ui.html \
  -o inst/doc/getting-started.html
```

Start with **`taruwind_downscale()`** when you have coordinates and a period.
You do not need to call all exported functions in sequence.

| Task | Entry point |
|---|---|
| Coordinates + period, automatic data preparation | `taruwind_downscale()` |
| Small local wind map for a period | `taruwind_downscale_area()` |
| Aligned elevation and roughness maps | `taruwind_site_layers()` |
| Draw an input or wind raster | `taruwind_plot_map()` |
| Previously prepared response + new period | `taruwind_downscale(responses = ...)` |
| Prepared response + your own wind components | `taruwind_predict()` |
| Georeferenced DEM and roughness rasters | `taruwind_prepare_site()` |
| Explicit metric-grid matrices | `taruwind_build()` |
| Independent metric-grid batches | `taruwind_build_many()` |
| Inputs only, without solving | `taruwind_download_era5()`, `taruwind_download_copernicus_dem()`, `taruwind_download_worldcover()` |
| Profile assumptions / resource policy | `taruwind_source_profile()`, `taruwind_neutral_depth()`, `taruwind_limits()`, `taruwind_workers()` |
| Inspect cached files and total size without deleting data | `taruwind_cache_info()` |

The geographic workflow reports six stages, beginning with an ERA5 availability
check before terrain preparation. It labels cached versus downloaded ERA5
chunks. Set `quiet = TRUE` to suppress progress. Stage numbers do not
estimate elapsed-time percentages. Explicit persistent caches are kept for reuse;
the default cache is session-temporary.

Requests must start in 1940 or later, end no later than the current UTC time,
and use complete hourly timestamps. The actual latest available ERA5 hour is
checked against the CDS ARCO time axis before terrain or source-profile work;
there is no hard-coded five- or seven-day availability lag. If cached metadata is behind
the source, use `refresh = "always"`. Each automatic downscaling request is
limited to ten calendar years, at most 100 sites and 10 million site-hours,
plus reflected terrain-grid limits. Local area maps have an additional
seven-day (168-hour) limit; reuse `attr(area$wind, "responses")` when batching
consecutive weeks.
Inspect `taruwind_limits()` for all thresholds. These are safety guards, not
an estimate or guarantee of network bytes, runtime or peak memory. Oversized
requests stop with an error and should be split into recoverable batches. The
standalone ERA5 download helper retains its deliberate `allow_large = TRUE`
override.

`wind_speed` in the geographic output is **vector-bilinear ERA5**;
`taruwind_speed` is the final downscaled speed, both in m/s. The `u`/`v`
columns remain ERA5 input components, not downscaled output components.
`u10`/`v10` (or `u100`/`v100`) are height-specific aliases of `u`/`v`.
To retrieve the four surrounding grid-node series, use
`taruwind_download_era5(include_nodes = TRUE)` and read `attr(era5, "nodes")`.
When reusing the main workflow's cache, pass `file.path(cache, "era5")` to this
helper. The guide includes a complete example and a separate node CSV export.

For a small area around the same example coordinate, request aligned static
inputs and hourly wind maps. A separate response is built for each output
cell, so begin with a coarse output grid:

```r
cache = file.path(tools::R_user_dir("taruwind", "cache"), "first-example")
input_maps = taruwind_site_layers(29.02, 41.02, roughness = "worldcover",
  domain_m = 990, cache_dir = file.path(cache, "dem"))
taruwind_plot_map(input_maps, layer = c("elevation_m", "roughness_m"))

area = taruwind_downscale_area(29.02, 41.02,
  start = "2020-01-01 00:00:00", end = "2020-01-01 02:00:00",
  roughness = "worldcover", area_m = 990, map_resolution = 500,
  domain_m = 990, workers = 1, cache_dir = cache)
taruwind_plot_map(area$maps, layer = "taruwind_20200101T000000Z")
```

The example creates a 3 × 3 output grid and uses a 990 m terrain domain as a
software check. The normal terrain-domain default is 30 km. The area result
contains the cell coordinates, site-hour table and time-labelled rasters.
Output-cell spacing and terrain-grid resolution are separate settings.

**Recommended order:** follow the one-coordinate guide → inspect output and
raw ERA5 nodes → save responses → map aligned inputs → try the local area
example → expand to a suitable domain and period. An initial terrain
solve is distinct from the much cheaper application of a prepared response.

## Technical overview and examples

Local first-release candidate by Cem Özen (ozenc@itu.edu.tr).
Prepared in the private GitHub repository `cemoze/taruwind` for cross-platform
release checks. Not yet published on CRAN or as a public release.

The current neutral reduced-mode terrain solver and finite-fetch surface
memory are exposed through an explicit metric-grid interface. Older research
algorithms are not silently substituted by geographic convenience functions.

```r
library(taruwind)
dem = matrix(0, 9, 11)
roughness = matrix(0.03, 9, 11)
response = taruwind_modal_prepare(
  dem, roughness, dx = 90, x0 = 0, y0 = 0, x = 450, y = 360,
  depth = 900, source_z0 = 0.03
)
taruwind_modal_predict(response, u = c(3, 0), v = c(4, 0))
```

The two predicted speeds are 5 and 0 m/s on flat, homogeneous,
matched-reference terrain.

Rows run south to north. Coordinates are projected metres; roughness values
are aerodynamic lengths in metres, not land-cover class codes. `source_z0`
is an explicit effective input, never inferred from target roughness. Wind
components must already be interpolated to the target in the grid frame.
Preparation uses one native thread. Observations are not used to fit responses.

`taruwind_build()` and `taruwind_predict()` use this same current core.
Prediction also accepts four-column source-node component matrices and four
interpolation weights. `taruwind_build_many()` supports portable PSOCK workers:
auto is half the logical cores, capped at 16 and the site count, with a padded-
cell memory guard. Each native solve remains single-threaded.

Separate ERA5 (10/100 m) and Copernicus DEM download utilities are included;
GLO-90 is the default. Network services and credentials remain provider-owned.
No download occurs when the package is loaded or during ordinary package tests.
`taruwind_prepare_site()` projects rasters to a target-centred metric grid;
`taruwind_downscale()` connects this current solver to ERA5 retrieval. DEM
acquisition is automatic unless a raster is supplied. Roughness may be a
georeferenced aerodynamic-length raster or explicitly `"worldcover"` to fetch
ESA WorldCover 2021 and apply the documented WindKit lookup. WorldCover is
not presented as a replacement for the study's CORINE-based inputs.

This geographic example requires CDS credentials and acceptance of the
relevant dataset terms. TaruWind obtains the documented reference profile
automatically:

```r
wind = taruwind_downscale(
  lon = 29.02, lat = 41.02,
  period = c("2020-01-01", "2020-01-02"),
  roughness = "worldcover"
)
later = taruwind_downscale(
  lon = 29.02, lat = 41.02,
  period = c("2020-02-01", "2020-02-02"),
  responses = attr(wind, "responses")
)
class(wind)
names(attributes(wind))
names(attr(wind, "responses"))
```

`wind` and `later` are tibbles, which also inherit from `data.frame`.
The `responses` attribute is a named list of prepared responses for reuse;
`forcing_provenance` records the ERA5 source. The HTML guide explains how to
inspect and save both the rows and their attributes. Complete DEM/roughness
rasters and four-node raw wind series are retrieved separately; they are not
embedded in the printed prediction tibble.

`sites` also accepts a data frame or matrix with `lon`/`lat` columns and an
optional `site_id`. Output is a tibble; UTC is default, `local_tz = TRUE` adds
coordinate-specific local timestamps with `lutz`. Geographic downloads and
raster preparation are serial to avoid cache races; response construction
uses automatic portable workers with a memory guard. The bounded local-area
helper computes independent responses at each output-cell centre. Missing
raster coverage is an error, never silently filled.

Source roughness is never guessed from target land cover. If omitted, it is
obtained from the ERA5 FSR reference climatology (default 2021, all months/hours)
through the documented diagnostic adapter. At 10 m this adapter bounds FSR to
0.0002–0.03 m; it does not assert that actual ERA5 roughness is globally 0.03 m.
The neutral depth follows the study's latitude-dependent formula. Both can
still be explicitly overridden with `source_z0` and `depth`.
Near the equator, automatic depth is rejected; this convention is restricted
to absolute latitudes 5–85 degrees. Source-profile downloads can enter a CDS
queue: repeat the same call after timeout to resume the cached job, without
submitting it again. Install optional `ncdf4` for this operation.
Automatic CORINE acquisition is not included in this release. The
geographic workflow is a new input-grid convenience layer, not a byte-for-byte
reproduction of the archived study grids. The current solver itself is unchanged.

The `.github/workflows/R-CMD-check.yaml` workflow checks Linux, Windows and
macOS. All three completed with `Status: OK` for the base source snapshot in
GitHub Actions run 35756209012 (22 September 2026). The input-validation
correction in this candidate has passed local Linux checks; its exact source
has not yet had hosted Windows/macOS checks. Package checks cover installation,
offline examples and tests, not live downloading on every platform.

The model is neutral and diagnostic; separation, thermal circulations and
canopy displacement are not explicitly resolved. The mapped source/target
profile convention must be chosen explicitly. Evaluation above 10 m is not
claimed by this interface; supported numeric heights extend up to 100 m.
