# PolygonTorus

Is a tool for creating axially symetric CSG geometry. By revolving polygons around an axis we get toroids -- a torus-like surface, which is a composite of conical, planar and cylindrical surfaces.

Currently it generates CSG geometry for openMC and MCNP. When generating MCNP geometry the tool can be used by itself, whereas when generating openMC geometry, openMC is used for the construction of surfaces and regions.

The tool is in its early stages:
 * There are two seperate modules for openMC and MCNP inputs.
 * Code is not well documented.
 * For MCNP: The tool generates only regions. The generation of cells is made easier with the `split_into_lines()` function, but is still cumbersome.

For a theoretical explanation see `polygontorus_presentation.pdf`. For examples see `examples` folder.

## Features

**Offset():** Because of its robust definition of geometry this tool can be used for arbitrarily thin layers. When generating thick layers inside the shells (negative offset), caution should be used since vertices may get deleted.

## Known limitations

* Generation of a Polygon with many verticies in MCNP. This can result in a region definition that is too big for MCNP (MCNP 5.16 was tested).
* The `split_into_lines()` function is bugy at times.
* Suports only convex polygons. Region splitting should be implemented (also to ensure cell definitions are more simple - see 1st limitaiton).

## Ongoing activities

* Performance assesment vs. toroidal surfaces
