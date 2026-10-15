## Geoid heights

The geoid's height above the ellipsoid, against the gravity model each crate implements, at GeographicLib's 500 000 random test points (NGA's harmonic synthesis, to 1 µm). egm96's harmonic synthesis runs on every 50th point: it took 561 µs per point here. Longitudes are given from −180° to 180°; egm2008 refuses the same places written as 180° to 360°, and geo3d takes either.

### Geoid height against the model

| Library | Model and method | Points | Failures | Max error (m) | RMS error (m) |
|---|---|---|---|---|---|
| geo3d | EGM96, egm96-5 grid (19 MB), cubic | 500000 | 0 | 2.6e-3 | 6.6e-4 |
| geo3d | EGM96, egm96-5 grid, bilinear | 500000 | 0 | 1.2e-1 | 4.6e-3 |
| geo3d | EGM96, egm96-15 grid (2 MB), cubic | 500000 | 0 | 1.7e-1 | 7.0e-3 |
| geo3d | EGM2008, egm2008-1 grid (470 MB), cubic | 500000 | 0 | 2.2e-3 | 6.5e-4 |
| egm96 | EGM96, its 5′ raster (the same egm96-5 grid), bilinear: its default | 500000 | 0 | 1.9e0 | 1.4e-1 |
| egm96 | EGM96, spherical harmonics to degree 360 | 10000 | 0 | 1.3e-4 | 4.1e-5 |
| egm2008 | EGM2008, its built-in 3° grid | 500000 | 0 | 2.4e1 | 1.5e0 |

