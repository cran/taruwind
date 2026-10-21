---
title: "TaruWind 0.1.0 — A first wind-downscaling workflow in Istanbul"
lang: en
---

## The story: one location, from coordinates to reusable predictions

You do not need to learn every function. **Start with `taruwind_downscale()`
when you have coordinates and a period.** It combines data acquisition,
raster preparation and the current solver. No target observations are required.

Suppose you want a wind-speed series at **29.02° E, 41.02° N in Istanbul**.
We will request the first three hours of 1 January 2020, compare the resulting
speeds, save our work and reuse the same terrain response for the following day.
This is an illustrative coordinate, not a claim of an instrumented station.

**Follow Sections 1–8 in order in the same R session.** Each step uses objects
created earlier. Installation commands are needed only once. The later reference
sections explain extensions; they are not extra steps required to complete the
example. No local NetCDF or raster files need to be prepared manually.

## 1. Prepare R and data access

Start with the software requirements, before defining a site or requesting
data. In RStudio, select **Session → Restart R** (`Ctrl+Shift+F10`). Check which
packages are already available; this downloads nothing:

```r
needed = c("taruwind", "terra", "ncdf4", "Rarr", "lutz")
setNames(vapply(needed, requireNamespace, logical(1), quietly = TRUE), needed)
```

Install missing dependencies as needed. `Rarr` is available through
Bioconductor; `BiocManager` is only its installer:

```r
cran_packages = c("terra", "ncdf4", "lutz")
missing_cran = cran_packages[!vapply(cran_packages, requireNamespace,
                                     logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran)
if (!requireNamespace("Rarr", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("Rarr", ask = FALSE, update = FALSE)
}
```

`terra`, `ncdf4` and `Rarr` support this automatic download workflow. `lutz`
is needed only if you later request automatic local time zones. TaruWind is
not yet on CRAN; use its locally installed package rather than
`install.packages("taruwind")`.

Now load the package and confirm which local version R is using:

```r
library(taruwind)
packageVersion("taruwind")
find.package("taruwind")
?taruwind
```

The package version shown here should be 0.1.0 for this guide.

### CDS access

Follow the [official CDS API setup instructions](https://cds.climate.copernicus.eu/how-to-api).
The package can read `~/.cdsapirc` or the `CDSAPIKEY` / `CDS_API_KEY` environment
variable. Accept the applicable terms for the ERA5 single-level and monthly-mean
datasets in your CDS account.

Existing credentials may already be sufficient. **Do not put your API key in
scripts, GitHub or chat messages.** DEM and WorldCover downloads do not require
a CDS key. The package does not require Python.

For the workflow below, these are the inputs and their sources:

| Input | What the package does | What you provide |
|---|---|---|
| Location and time | Selects the target and hours | Longitude/latitude in decimal degrees and a UTC period |
| ERA5 | Downloads wind from the four surrounding nodes | CDS account and API key |
| DEM | Downloads GLO-90 and projects it onto a metric grid | Optionally, your own DEM raster |
| Target roughness | Converts WorldCover classes to z₀ | `roughness = "worldcover"` or a prepared z₀ raster |
| Source profile | Retrieves reference FSR and applies the model convention | Optionally, explicit `source_z0` |
| Neutral depth | Calculates a latitude-dependent depth | Optionally, explicit `depth` |

DEM values represent **terrain elevation in metres**. Roughness values represent
**aerodynamic z₀ in metres**, not land-cover class identifiers.

## 2. Define our Istanbul site and period

With the R packages and CDS access ready, define the location and requested
hours once. All later steps use these same objects. Times are UTC and both
endpoints are included: this request contains **three**, not two, hourly records.

```r
site = data.frame(site_id = "istanbul_example", lon = 29.02, lat = 41.02)
start_utc = "2020-01-01 00:00:00"
end_utc = "2020-01-01 02:00:00"
cache = file.path(tools::R_user_dir("taruwind", "cache"), "first-example")
output_dir = file.path(getwd(), "taruwind-istanbul-example")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
```

`cache` holds reusable downloaded inputs; `output_dir` will hold our saved results.
Print either variable to see its full path. No files are deleted by this example.

## 3. Obtain the first three hours at our site

This example uses a terrain domain approximately 1 km wide. **The 990 m domain
is only an installation/integration check, not the normal study domain.**

By default, downloaded files are cached in the current R session's temporary
directory. The example below explicitly chooses a persistent user cache so
downloads and queued CDS requests can be reused after restarting R.

```r
wind = taruwind_downscale(
  sites = site,
  start = start_utc, end = end_utc,
  roughness = "worldcover",
  domain_m = 990, resolution = 90, workers = 1,
  cache_dir = cache
)
wind[, c("timestamp_utc", "wind_speed", "taruwind_speed")]
```

`wind` is a **tibble**, which is also an R data frame. Each printed row is
one target and one hour. It additionally carries a prepared terrain response
and forcing provenance as attributes; these do not appear among the printed
columns. Section 5 shows exactly how to inspect them.

The first call includes network access and potentially a CDS queue, so its
elapsed time is not just solver time. Requesting three hours does not necessarily
mean downloading a three-hour-sized file: provider chunks may be larger.

If you receive **“CDS source-profile job is still queued/running”**, repeat the
same call later with **the same cache directory**. The saved job identifier
allows the request to resume rather than submitting another job. Do not keep
clearing the cache and starting over.

For a normal site analysis, remove `domain_m = 990`: the default domain is
approximately **30 × 30 km**, with **90 m** metric cells. This requires more
computation; initial preparation at complex sites can take several minutes.
**Shortening the period does not shrink the terrain solve.** A 30 km domain
is not a universal guarantee of adequacy for every landscape.

## 4. Follow what happened during that call

### Progress and downloaded files

The console reports six stages: ERA5-hour availability, source profile,
terrain/roughness preparation, response construction, ERA5 acquisition and
response application. These are
stage counters, not percentages of elapsed time. Use `quiet = TRUE` in
`taruwind_downscale()` to suppress these messages. Warnings and errors remain.

ERA5 winds are read from Zarr, not downloaded as GRIB. The messages distinguish
coordinate arrays from wind chunks and report cache reuse. A stored wind chunk
can cover many more hours than requested; only the requested hours enter the
output. `Rarr` reads these compressed arrays; `BiocManager` is only its installer.
Neither is needed when applying the solver to your own wind components.

Inspect the example cache without deleting anything:

```r
cache_info = taruwind_cache_info(cache)
cache_info
attr(cache_info, "total_bytes") / 1e6
attr(cache_info, "cache_dir")
```

The divided total is in megabytes; the cache directory is the path used above.

An explicit persistent cache is retained after the call and after R restarts.
The default session-temporary cache is normally removed when R closes, not
when the function returns. The package does not automatically evict persistent
files. Reusing cached data does not itself reuse the terrain solve; save and
pass `responses` as shown below to avoid reconstructing that response.

DEM acquisition currently caches complete GLO-90 tiles. WorldCover acquisition
reads raster windows. A short wind period does not reduce these static inputs.

## 5. Inspect and plot our first result

The `wind` tibble now contains our three Istanbul records. The two speed
columns are the input baseline and the model output, not two observation series.

Inspect the object in layers, without printing its potentially large stored
response:

```r
class(wind)
nrow(wind)
names(wind)
names(attributes(wind))
names(attr(wind, "responses"))
prepared = attr(wind, "responses")[[site$site_id]]
class(prepared)
names(prepared)
names(prepared$provenance)
attr(wind, "responses")[[site$site_id]]$provenance$geographic$actual_domain_m
attr(wind, "forcing_provenance")
```

`class(wind)` includes `tbl_df` and `data.frame`. `names(wind)` lists the
ordinary data columns below. `names(attributes(wind))` also lists standard
R data-frame attributes such as `names`, `row.names` and `class`. The current
workflow carries these additional attributes:

| Attribute | What it contains | Why keep it |
|---|---|---|
| `responses` | A named list with one prepared terrain response per site | Reuse it for later hours without rebuilding that site's terrain response |
| `forcing_provenance` | Source, retrieval and interpolation metadata for ERA5 | Record which forcing produced the rows |
| `sites`, `source`, `license`, `citation`, `retrieved_utc`, `wind_height_m`, `component_variables`, `interpolation`, `cache_dir` | Metadata forwarded from the ERA5 download result | Inspect the target and downloaded forcing |

`prepared` is one `taruwind_modal_response` object. Its `matrix`, `surface`,
`residual`, `iterations` and `provenance` components describe the reusable
response; it is not another hourly data frame. For this quick check,
`actual_domain_m` is 990 m. To discover a result's contents in R, start with
`class()`, `names()` and `names(attributes())`; `?taruwind_downscale` documents
the return value. An RDS file preserves both the columns and these attributes;
CSV preserves only the columns. The result does not embed complete DEM or
roughness rasters or all four raw ERA5 node series. The source files remain in
the cache; Sections 5 and 8 show how to retrieve the nodes and aligned maps.

| Column | Meaning |
|---|---|
| `site_id`, `lon`, `lat` | Target identifier and coordinates |
| `timestamp_utc` | UTC timestamp |
| `wind_height_m` | Requested height: 10 or 100 m |
| `u`, `v` | ERA5 input components interpolated to the target using vector bilinear interpolation |
| `wind_speed` | ERA5 vector-bilinear speed, m/s |
| `taruwind_speed` | Final TaruWind speed estimate, m/s |

The `u10`/`v10` or `u100`/`v100` columns are also **input** components, not
downscaled output components. This version returns speed as its final product.
At 10 m, `u10` equals `u` and `v10` equals `v`: these are height-specific
aliases, not additional predictions or surrounding grid values. To inspect a
compact view without removing attributes from the original object:

```r
wind_view = wind[, c("site_id", "timestamp_utc", "wind_height_m",
                     "u", "v", "wind_speed", "taruwind_speed")]
wind_view
```

Specifically, `wind_speed = sqrt(u^2 + v^2)` **after** interpolating the
components. It is vector-bilinear speed, not scalar interpolation of the four
node speeds. `taruwind_speed` is the final model estimate.

```r
wind$increment = wind$taruwind_speed - wind$wind_speed
plot(wind$timestamp_utc, wind$wind_speed, type = "l", col = "#0072B2",
     xlab = "UTC", ylab = "Wind speed (m/s)",
     ylim = range(c(wind$wind_speed, wind$taruwind_speed), finite = TRUE))
lines(wind$timestamp_utc, wind$taruwind_speed, col = "#D55E00")
legend("topright", c("ERA5 Bilinear", "TaruWind"),
       col = c("#0072B2", "#D55E00"), lty = 1, bty = "n")
```

This uses ordinary R graphics, not a special TaruWind plotting function. In
RStudio the figure appears in the Plots pane. It is a short demonstration, not
a climatology. You can use the same tibble with ggplot2 for longer series.

A difference between models is not itself an accuracy gain. Independent
observations are needed to assess accuracy.

### Retrieve the four ERA5 input nodes for these same hours

The downloaded data are not limited to the target prediction. The ERA5 helper
can return each surrounding node's wind components, speed, coordinates and
interpolation weight. Reuse the **`era5` subdirectory** of our existing cache:

```r
era5 = taruwind_download_era5(
  sites = site, start = start_utc, end = end_utc, wind_height = 10,
  include_nodes = TRUE, cache_dir = file.path(cache, "era5")
)
nodes = attr(era5, "nodes")
nodes[, c("site_id", "timestamp_utc", "node", "node_lon", "node_lat",
          "weight", "u", "v", "wind_speed")]
saveRDS(era5, file.path(output_dir, "era5-with-nodes.rds"))
write.csv(as.data.frame(nodes), file.path(output_dir, "era5-four-nodes.csv"),
          row.names = FALSE)
```

For this example there are 12 rows: four nodes for each of three hours.
These are native ERA5 grid values exposed by the data store, before our
interpolation or terrain correction. Cached historical wind chunks are reused;
metadata may still be checked. `wind` itself does not request these raw series:
`include_nodes` currently belongs to `taruwind_download_era5()`, not
`taruwind_downscale()`. RDS preserves the node attribute; exporting only the
main `era5` tibble to CSV would omit it, hence the separate node CSV above.

## 6. Save our result and predict the following day

We now retain both the output and the terrain response built for this exact
location. The response can be applied to new hours without solving the same
terrain again. Download caching alone would not achieve this.

```r
responses = attr(wind, "responses")
saveRDS(responses, file.path(output_dir, "responses.rds"))
saveRDS(wind, file.path(output_dir, "first-three-hours.rds"))
write.csv(as.data.frame(wind), file.path(output_dir, "first-three-hours.csv"),
          row.names = FALSE)

later = taruwind_downscale(
  sites = site,
  period = c("2020-01-02", "2020-01-02"),
  responses = responses, cache_dir = cache
)
nrow(later)
saveRDS(later, file.path(output_dir, "following-day.rds"))
```

Identical start/end dates select all 24 hours of that UTC day; endpoints are
inclusive. This reused response still belongs to the **990 m test domain**.
Reusing it does not enlarge the domain. Build a new response for a normal 30 km
analysis.

In another session, recreate `site`, `cache` and `output_dir` from Step 2,
then load `responses = readRDS(file.path(output_dir, "responses.rds"))`
before running the `later` call.
CSV does not preserve model/provenance attributes; use RDS for that purpose.

Rebuild whenever you change the DEM, roughness, resolution, domain, coordinates,
height or source profile. Supplying new raster arguments alongside `responses`
does not update the saved response.

## 7. Inspect the inputs behind this same prediction

Our first workflow is complete. If you want to understand its source-profile
inputs, inspect the cached profile for the **same Istanbul coordinate**. This
does not require a file called `my-era5-fsr.nc` or a second manual setup.

```r
profile = taruwind_source_profile(
  lon = site$lon, lat = site$lat,
  cache_dir = file.path(cache, "source_profile"), quiet = FALSE
)
profile$source_z0
profile$depth
responses[[1]]$provenance
attr(wind, "forcing_provenance")
```

These describe the inputs and response that generated our result. They are
not instructions to fit a new model to observations. The source-reference
convention is described under Assumptions below.

## 8. Map the inputs and a small surrounding area

First, look at the **same elevation and roughness layers** used when preparing
our Istanbul response. `taruwind_site_layers()` uses the same alignment routine
as the solver. It reuses downloaded tiles in our cache and does not acquire
wind. The two layers are metres of elevation and metres of aerodynamic
roughness length, respectively.

```r
input_maps = taruwind_site_layers(
  lon = site$lon, lat = site$lat, roughness = "worldcover",
  domain_m = 990, resolution = 90,
  cache_dir = file.path(cache, "dem")
)
names(input_maps)
taruwind_plot_map(input_maps, layer = c("elevation_m", "roughness_m"))
```

Now extend the example from one target to a **3 × 3 output grid** around it.
`map_resolution = 500` is the output-cell spacing; `resolution = 90` remains
the terrain calculation spacing. The requested 990 m square rounds to an
actual 1,500 m square so the centre itself is an output cell. Every output
cell has its own local terrain response. The nine cells therefore require
more work than the original point call, even though the ERA5 input hours and
the source files are cached.

```r
area = taruwind_downscale_area(
  lon = site$lon, lat = site$lat,
  start = start_utc, end = end_utc,
  roughness = "worldcover",
  area_m = 990, map_resolution = 500,
  domain_m = 990, resolution = 90, workers = 1,
  cache_dir = cache
)
names(area$maps)
area$grid$actual_area_m
taruwind_plot_map(area$maps, layer = "taruwind_20200101T000000Z")
taruwind_plot_map(area$maps, layer = "bilinear_20200101T000000Z")
saveRDS(area, file.path(output_dir, "three-by-three-area.rds"))
```

`area$wind` contains the underlying site-hour tibble and `area$maps` contains
one TaruWind and one vector-bilinear ERA5 raster per hour. Inspect
`area$grid$sites` for the exact output coordinates. The map is a grid of
independently computed point estimates, not a single shared-domain flow solve.
Unlike `wind`, **`area` is a list**, so inspect its named components with
`names(area)`. Its `wind` component is a tibble, its `maps` component is a
`terra::SpatRaster`, and its `grid` component describes cell locations and
resolution:

```r
class(area)
names(area)
class(area$wind)
class(area$maps)
area$grid$sites
```

The 990 m terrain domains and nine coarse map cells demonstrate the software;
use an appropriate terrain domain and output spacing for an analysis. The
function defaults to at most 25 output cells and permits at most 100 by
explicitly raising `max_cells`; work scales with the number of cells. An area
call is also limited to 168 hourly timestamps (seven days). For the next
period, reuse the nine prepared responses instead of solving terrain again:

```r
area_next = taruwind_downscale_area(
  lon = site$lon, lat = site$lat,
  start = "2020-01-01 03:00:00", end = "2020-01-01 05:00:00",
  area_m = 990, map_resolution = 500,
  domain_m = 990, resolution = 90,
  responses = attr(area$wind, "responses"),
  cache_dir = cache
)
```

You can export a raster with `terra::writeRaster(area$maps[[1]], filename)`.
To save a PNG directly, supply a new filename to `taruwind_plot_map(...,
file = "my-map.png")`; it will not overwrite an existing image.

## Optional next step: use the normal terrain domain

Before interpreting this site scientifically, build a response with the normal
domain. This optional calculation uses the same site, hours and cache but
**constructs a new approximately 30 km terrain response**. It can take
substantially longer than the preceding steps.

```r
wind_30km = taruwind_downscale(
  sites = site, start = start_utc, end = end_utc,
  roughness = "worldcover", domain_m = 30000, resolution = 90,
  workers = 1, cache_dir = cache
)
class(wind_30km)
wind_30km[, c("timestamp_utc", "wind_speed", "taruwind_speed")]
attr(wind_30km, "responses")[[site$site_id]]$provenance$geographic$actual_domain_m
saveRDS(wind_30km, file.path(output_dir, "normal-domain-result.rds"))
```

`wind_30km` has the **same tibble/data-frame structure** as `wind`; the
different variable name does not create another result type. Their ERA5
input hours and target are the same. The stored `responses` attribute shows
which domain was used to construct each terrain response. Compare the two
model series directly:

```r
data.frame(timestamp_utc = wind$timestamp_utc,
           quick_check_speed = wind$taruwind_speed,
           normal_domain_speed = wind_30km$taruwind_speed)
```

The earlier `wind` remains the 990 m software check; `wind_30km` is the
separate normal-domain result. Saving either with `saveRDS()` keeps its
columns and attributes for later inspection.

### Multiple sites, local time and 100 m winds

Now suppose we also want a second location 0.02 degrees northeast of our first
site. Keep the original identifier and add the second coordinate. This optional
example keeps the same short period and small software-test domain:

```r
sites = rbind(site, data.frame(site_id = "istanbul_neighbour",
                              lon = site$lon + 0.02, lat = site$lat + 0.02))
network = taruwind_downscale(
  sites = sites, start = start_utc, end = end_utc,
  domain_m = 990, resolution = 90,
  roughness = "worldcover", workers = "auto", local_tz = TRUE,
  cache_dir = cache
)
network[, c("site_id", "timestamp_utc", "wind_speed", "taruwind_speed")]
saveRDS(network, file.path(output_dir, "two-site-result.rds"))
```

Expect six rows: two sites times three hours. This constructs responses for
both targets; it does not implicitly reuse the earlier `responses` object.
For normal analysis use 30 km domains as above. A matrix with named `lon` and
`lat` columns is also accepted. Keep `site_id` values unchanged when reusing
responses.

`local_tz = TRUE` adds `time_zone` and `timestamp_local` while retaining UTC.
**The requested period is still interpreted in UTC**, not local time.

For 100 m winds, make a new call with `wind_height = 100`. Native ERA5 u100/v100
are retrieved, rather than extrapolating from 10 m. Do not reuse a 10 m response
at 100 m. Numerical support does not imply observational validation at 100 m
in the accompanying study.

## Reference: supplying your own inputs

**The first-use workflow is already complete.** `taruwind_downscale()` handles
the source profile internally. Do not run a separate profile or raster setup
unless you want to inspect inputs or supply your own data.

### Supply existing files

If you already have a georeferenced DEM and an **aerodynamic z₀ raster in
metres**, pass those real raster objects or their real file paths as `dem`
and `roughness` to `taruwind_downscale()`. It can still obtain the reference
profile automatically. There are no placeholder files to copy or rename in
this guide. See `?taruwind_prepare_site` for the lower-level input interface
and `?taruwind_predict` for applying a prepared response to wind components.

Rasters must cover the complete local domain. Convert CORINE class codes using
an appropriate z₀ lookup before supplying them as roughness. Automatic CORINE
downloading is not provided. WorldCover inputs do not reproduce the study's
CORINE-based experiment exactly.

`taruwind_build()` is the advanced metric-matrix interface: x/y are in metres,
rows run south to north and columns west to east. Converting a terra raster
directly to a matrix may reverse the required row order. The geographic
interface handles this orientation.

## Reference: assumptions, runtime and limits

- The model is a neutral diagnostic approach. Turbulence, thermal circulations
  and canopy displacement are not explicitly resolved.
- GLO-90 is projected onto a metric 90 m grid. DEM resolution, 10/100 m wind
  height and 30 km domain width are separate quantities.
- WorldCover 2021 classes receive z₀ values from a documented lookup; these
  are not direct roughness measurements.
- The default source FSR reference year is 2021. It does not automatically
  follow the wind period; change it with `profile_year`. At 10 m, the effective
  source reference is bounded to 0.0002–0.03 m. This does not claim that actual
  ERA5 roughness is everywhere 0.03 m.
- Automatic neutral depth is not hourly boundary-layer height. The automatic
  profile path is restricted to absolute latitudes of 5–85 degrees.
- Auto workers use half the logical cores, with a minimum of 1 and maximum of
  16, further limited by the number of sites and memory guards. **Sixteen
  workers do not split one site's solve into sixteen parts.** Data preparation
  is sequential; independent numerical solves run in parallel.

```r
attr(wind, "responses")[[1]]$provenance
attr(wind, "forcing_provenance")
taruwind_limits()
taruwind_workers(n_sites = 10)
```

The automatic point workflow allows up to ten calendar years per call, 100
sites and 10 million site-hours.
Reflected solves are limited to 2 million computational cells per site and
60 million cells across a request. These cells include the solver's reflected
grid, not just the visible local square. The local-area workflow defaults to
25 output cells and permits at most 100 when `max_cells` is raised explicitly;
its period is limited to 168 hourly timestamps (seven days).
`allow_large` belongs to the separate ERA5 download helper, not the main
`taruwind_downscale()` function. Split larger requests into batches.

The site-hour limit still depends on the number of sites. Longer point series
and area-map periods can be processed in consecutive batches with reusable
responses. Requests before 1940, after the current
UTC time, or outside the hours currently present in the CDS ARCO record stop
with errors. The source commonly trails real time by several days, but no
fixed lag is assumed: TaruWind checks the actual time axis before downloading
terrain or requesting a source profile. By default, mutable metadata is
refreshed once per day; use `refresh = "always"` to check the latest source
state immediately. The package does not impose a network-byte quota or
promise that a request below its safety limits will be fast or fit in memory.

## Reference: troubleshooting

| Problem | First check |
|---|---|
| An older version loads | Restart R; inspect `find.package()` and `packageVersion()` |
| CDS key / 403 error | API configuration and acceptance of both datasets' terms |
| Queue timeout | Retry later with the same call and cache directory |
| Requested ERA5 hour unavailable | Inspect the available range in the error; use `refresh = "always"` if cached metadata may be stale |
| Request exceeds site-hour or reflected-cell limits | Reduce the period, site count or terrain domain; run recoverable batches |
| Missing R package | Run the dependency check in Section 1 |
| Raster gaps / unknown class | Check complete domain coverage; missing ocean cells are not automatically set to zero |
| Response coordinate/height mismatch | Match the target, `site_id` and `wind_height` |
| Slow initial preparation | Check terrain complexity and domain size; separate download time |
| NA output | Check input u/v for missing values |

This version produces point wind-speed time series and bounded local raster
maps, not national atlases, interactive maps, downscaled direction/u/v or
CERRA downloads. Linux, Windows
and macOS package checks passed in GitHub Actions run 35756209012 on
22 September 2026. Live data integration was additionally checked on Linux;
the platform checks do not download environmental datasets.

```r
?taruwind_downscale
?taruwind_prepare_site
?taruwind_predict
help(package = "taruwind")
browseURL(system.file("doc", "getting-started.html", package = "taruwind"))
```

In short: **start with one site, understand the output, save the response,
then scale up.**
