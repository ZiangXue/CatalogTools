# CatalogTools.py User Manual
2020.08.06

Ziang Xue\
Pomona College\
zxaa2018@mymail.pomona.edu

## Intro

The CatalogTools.py program provides a set of tools to load and query through catalogs with the structure of [Hubble Catalog of Variables (HCV)](https://arxiv.org/abs/1909.10757).

For the complete HCV catalog file, download the full catalog in .csv format at [the HCV explorer](http://hst.esac.esa.int/hcv-explorer/). [Bonanos et al.(2019)](https://arxiv.org/pdf/1909.10757.pdf) describe how the HCV is made and structured, and also what each entries means.

To run the program, you must have the following python packages on your device: `progress`, `astroquery`, `matplotlib`, `numpy`, `math` and `scipy`.

For Mac users, to install these packages, use the `pip` command, or `pip3` for Python3. For example, to install `progress`:

```
pip3 install progress
```

## General

Based on how the HCV is made and structured, the program sorts the data of the .csv file into 4 classes `Catalog`, `Target`, `Lightcurve` and `Entry`. The order of these four classes here shows their levels, with the former having a higher level than the latter.

Each line in the .csv file is an entry representing a HST observation, or one data point in a lightcurve. `Entry` objects from the same lightcurve (observations of the same target with the same filter) are gathered into a list belonging to a `Lightcurve` object. `Lightcurve` objects of the same target are gathered into a list belonging to a `Target` object. Each `Target` object and its matchID forms a dictionary belonging to the `Catalog` object, with matchID as key and the `Target` object as value. `Catalog` objects can represent the HCV catalog, a subcatalog of the HCV, or any catalog with the same columns.

The four classes have other attributes to store information as well.

- `Catalog`
   - `name`: Name of the catalog
   - `catalog`: The dictionary of matchID and `Target`

- `Target`
   - `matchID`: The matchID of the target, its only identifier.
   - `groupID`
   - `subgroup`
   - `ra`
   - `dec`
   - `pipeline_class`
   - `expert_class`
   - `lightcurves`: A list containing all of its `Lightcurve` objects.

- `Lightcurve`
   - `entries`: A list containing all of the `Entry` objects in this lightcurve.
   - `matchID`
   - `filter_type`: The filter and instrument used in this lightcurve.
   - `var_quality_flag`
   - `filter_detection_flag`
   - `num_in_lc`: The number of data points in this lightcurve.
   - `hsc_mean_m`
   - `hcv_mean_m`
   - `mad`
   - `chi2`

- `Entry`
   - `matchID`
   - `filter_type`
   - `lightcurve_d`
   - `lightcurve_m`
   - `lightcurve_cm`
   - `lightcurve_e`
   - `lightcurve_i`
   - `lightcurve_r`
   - `ci_d`
   - `ci_v`
   - `d_d`
   - `d_v`

   
For attributes that are not explained, refer to [Bonanos et al.(2018)](https://arxiv.org/abs/1909.10757).

## File I/O

### Loading a .csv file into a `Catalog`

The constructor of a `Catalog` objects takes a string of its name as a parameter.

The function `catalog_file` of the `Catalog` class takes a file path as a parameter. It read and load the .csv file in the given path into a `Catalog` object with all the `Target`, `Lightcurve` and `Entry`.

There is an example to construct and load the full HCV catalog.

```python
hcv_full=Catalog('HCVFullCatalog')
hcv_full.catalog_file('/Users/ziang/Documents/RAISE2020/Program/1564584578352RG-result.csv')
```

### Building a `Catalog`

While to begin any research, one must load the full HCV or a given subcatalog, one can also build catalogs as results of query or filtering.

Note that since `Targets` are what users are mostly looking for rather than individual lightcurves and entries, the program does not support adding only `Entry` or `Lightcurve` objects.

The `add_target` function of the `Catalog` file allow user to add a **complete**, **well-defined** `Target` object into the `Catalog` called upon.

For example, `some_target` is a `Target` object with well-documented lightcurves and entries that we have sorted out from the full catalog `hcv_full`. The following codes create a new `Catalog` and add this `Target` in.

```python
hcv_new=Catalog('hcv_new')
hcv_new.add_target(some_target)
```

### Exporting a `Catalog` into a .csv file

The function `write_catalog` takes a file path as a parameter (path does not have to exist). It writes the `Catalog` object called upon into a .csv file.

Here is an example to export the subcatalog we have just made.

```python
hcv_new.write_catalog('hcv_new.csv')
```

## Operating with `Catalog`s

### Getting informations from the data

CatalogTools.py comes with some built-in functions that allow the user to retrieve informations about the objects in question, for example the time baseline length of a given `Lightcurve` object, or calculating the galactic coordinates of a `Target`. It also has built-in functions to plot lightcurves and retrieve FITS images.

Some earlier tasks such as plotting the bar charts of the entry distribution with regards to filter type are wrote as built-in functions as well. However, it is recommended to write these as independent script that import CatalogTools.py.

Please refer to the source code and doc strings for information on operating these functions.

One can also define their own functions taking `Catalog`, `Target`, `Lightcurve` and/or `Entry` objects as parameter to return a desired information.

### Querying the `Catalog`

Under most circumstances the program does not support operations directly on lower-level classes: Any query or operation starts with iterating over the `Catalog`, often using a `for` loop or a `while` loop.

To perform a search across each and every entry of the catalog, refer to the following example, using the `hcv_full` aforementioned:

```python
for target in hcv_full.catalog.values():
    for lightcurve in target.lightcurves:
        for entry in lightcurve.entries:
            #operation
```

Note that we do not have to always do this three-layer iteration: if we are only looking for a `Target` with some characteristic we can perform only the first, and if we are looking for a certain kind of `Lightcurve` we can perform only the first two and omit the last.

### Searching by matchID

The dictionary structure of the `Catalog` class allows direct access by key. For example, to access the target with *matchID=352* in the full HCV catalog, one can simply use the following codes:

```python
hcv_full.catalog[352]
```

## Future development

While the program is now fully capable of catalog query and operations, there are some attributes that requires implements.

One problem that lies within this program is that there is no direct access from a lower-level object to its father-object. The key search feature of `Catalog` class will be implemented to the lower-level classes as well, allowing a direct search for lightcurves and entries, which solves this problem.

For any need and desired functions, please contact me directly.
