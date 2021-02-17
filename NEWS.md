# SimuNet 1.3.0
* Added a homemade class of packed matrices: `snPackMat` (SimuNet Packed Matrix) that stores only a vector of the relevant data, as well as a way to unpack the matrix to a regular one. Especially designed to avoid storing empry matrix triangles and performing useless operations on them. Now `simu_scan` has an optional use.snPackMat logical argument (`FALSE` by default). Early benchmarks are promising  
* Added two user-friendly functions to interact with the [Animal Social Network Repository](http://www.github.com/bansallab/asnr):  
    * `import_from_asnr` can mostly be used with:  
        * the Class and species folder (in the Networks folder of the github repository), and possibly graphml file name, as character strings  
        * the URL of a graphml object within the repository  
    an option allows the user to easily retrieve a igraph network object, or an adjacency matrix, from the graphml file
    * `asnr_network_df` retrieve the list (into a dataframe) of graphml file in the asnr repository. Internally used by `import_from_asnr`
        

# SimuNet 1.2.0.9000
* Switched internally to integer matrices instead of numeric matrices: significant improvement of -10% computation time and -30% memory allocation on a `n = 21` `total_scan = 9000` empirical network
* Added dependency/reliance on the `Matrix` package, notably for its printing of sparse matrices
* Improved handling of triangular matrices. `presenProb` and `obsProb` objects now rightfully store triangular matrices and supposedly only the required random draw are performed (TO CHECK BEFORE TRUE RELEASE)

# SimuNet 1.1.0
* implemented a plot method for `scan` and `empiScan` objects. Internal code probably cleanable, but working for now:  
    `plot()` can be used either on a `scan`/`empiScan` object (which will be passed through `summary()`) or on a `summary.scan`/`summary.empiScan` one. Notable `plot` arguments available: `method` to choose within `c("both","theoretical","group","focal")` for `empiScan` objects. Special case of the `layout` argument used to pass igraphs `layout` or `layout_` functions internally. This way, a layout is determined before calling `plot.igraph` (wrapped in `plot_emprical`) in the case of `method = "both"`, to ensure that all networks rely on the same layout to ease visual comparison.

# SimuNet 1.0.1
* fixed some ugly print outputs for long `focalList` objects

# SimuNet 1.0.0
* v1.0.0 Release of a fully functional and documented version of this network simulation framework, including a shift toward OOP internally

# SimuNet 0.6.0.9000
* Added track of `X.scaled` (X being `"theoretical"`, `"group"`, or `"focal"`) for `summary.scan` and `summary.empiScan` objects

# SimuNet 0.5.0.9000
* Renamed more explicitely the `X.scan` components of scan and empiScan object, to make way for `X.sum` and `X.sampled` new matrices (cf. below).
* The `summary`method used on `scan` and `empiScan` objects can now be used to create objects of new `summary.scan` and `summary.empiScan` objects, storing `X.sum` and `X.sampled` (X being `"theoretical"`, `"group"`, or `"focal"`), and with dedicated print methods
* `Adj` and `total_scan` arguments are now optional when `sampling.param` is passed when using `simu_scan`
* Added several required function for this purpose:
    * `sum_scan.list` is an equivalent of previous `sum_up.scan(s)`, but cleaner. It sums up any list of scans counting `NA`s as zeros.
    * added `resolve_NA`,  which is called before `sum_scan.list` in the case of empirical scans. Moved all `NA`s related functions to `empirical_NA_tools.R`
    * added `sum_scan.sampled` and `count_NA` functions to also count sampled (non-`NA`) edges

# SimuNet 0.4.0.9000
* Included `scans.to.do` variable (a scan index _or_ a vector of `1:total_scan`) in wrappers and nested functions to generate `scan` and `empiScan` objects. Now `scans.to.do` is stored in several internal and output objects
* Vectorized wrappers and nested functions to produce lists of scans (`raw.scan`, `theoretical.scan`, `group.scan` or `focal.scan` empirical scans)
* Updated several print methods to display vectorized or list components of `scan`, `empiScan` and `samplingParam` objects, truncate them when required to not overload output

# SimuNet 0.3.0.9000
* Added a wrapper to `generate_samplingParam`. User should interact with `simu_samplingParam` to create `samplingParam` objects for their simulations.
* Vectorized the `focal` and `scans.to.do` (which now replace `scan.number` for more transparency) components of `focal` objects. Next version should also vectorize some components of `samplingParam`, `obsProb`, `scan`, and `empiScan` objects
* Finalized a working version of `simu_scan` for single theoretical and empirical scan. Next version should allow the user to directly use `simu_scan` to generate a list of either types of scans (thus superseding former `iterate_scans`).
* Homogenize some function and variable names.


# SimuNet 0.2.0.9000
* First shift toward integration of OOP elements into the simulation framework:
    * Revamp from the ground up of the `do.scan` part of the non-OOP previous version
    * Attempt at homogenizing SimuNet's internal syntax:
        * the wrappers - to perform the whole network simulations (former `Boot_scan()`), iterate single scans into weighted adjacency matrix (former `iterate_scans()`), and draw a single random binary scan (former `do.scan()`) - will now follow the syntax `simu_*()` (starting with the `simu_scan()` wrapper to supersede `do.scan()`)
        * new simulation internal objects (cf. below) will have generator following the syntax `generate_className()`, and as needed `print.className()` and other relevant S3 class methods.
        * presently, their related code use variable names distinct from their `className` syntax, rather favoring a `variable = class.name` syntax internally
        * function names including internal ones are (hopefully) more explicit, and are now action verbs to follow coding "grammar" recommendations
* For the OOP-transition purpose, creation of several simulation objects (S3 class):
    * `presenceProb` objects, generator and related S3 methods: calculate and store infos on the presence probability `P` of a tie at each scan for each dyad, from inputted adjacency matrix (`Adj`), sampling effort (`total_scan`), and igraph network `mode`
    * `samplingParam` objects, generator and related S3 methods: store all usefull parameters, some in the form of new object classes (S3), related to the empirical aspect of the network simulation. Presently, stores what's needed for a single binary focal scan, i.e. infos on a single focal (out of a list of focals). Specifically, `samplingParam` objects store:
        * `method`: a character scalar between `"group"`,`"focal"` and `"both"`, indicating the chosen scan sampling method.
        * `mode`: a character scalar representing the igraph network `mode`.
        * `obsProb` objects (cf. below)
        * `focal` objects (cf. below)
    * `obsProb` objects, generator and related S3 methods: calculate and store a probability of observation `P` of an edge at each scan sampling for each dyad, from inputted user-defined function `obs.prob_fun` of `(i,j,Adj)` (that can also be used to pass a single [0,1] numeric to use as a constant, or the string `"random"` to have all dyad probability drawn from `runif(n*n,0,1)`)
    * `focalList` objects, generator and related S3 methods: draw and store a list `focals` of `total_scan` focals (used internally as their indices, but also keep track of their names), from inputted user-defined function `focal.prob_fun` of `(n,Adj)` (that should return a vector of `n` probability for each node to be drawn at each scan, but that can also be used to pass as the strings `"random"` or `"even"` to have them drawn from `sample()` or to maximize the evenness of drawn focals in the focal list)
    * `focal` objects, generator and related S3 methods (including `plot.scan` S3 method to plot the network relying on `plot.igraph`): determine and store which `focal` to sample, from inputted `focalList` object (cf. above) and `scan.number` (out of the sampling effort `total_scan`)
    * `scan` objects, generator and related S3 methods: draw and store, from inputted `presenceProb` object:
        * a `raw` scan: directed binary matrix drawn from the probability of presence of edges contained in a `presenceProb` object
        * a `theoretical` scan: the `raw` scan to which the selected mode has been applied
        * a `scan.type` scan string: for `scan` object not sampled from yet (i.e. not empirical scans of the `empiScan` class (cf.below)), `scan$scan.type = "theoretical"`, but later in the process can become `"empirical"`
        * the original adjacency matrix `Adj`, sampling effort `total_scan`, selected igraph network `mode`, the probability matrix `presence.prob` (from `presenceProb$P`), and other parameters mostly for internal use
    * `empiScan` objects, generator and related S3 methods: draw and store, from inputted `scan` and `samplingParam` objects:
        * inherits from the `scan` S3 class
        * but have their `scan.type = "empirical"`
        * a `method` string from `samplingParam$method`
        * and both/either a `group` and/or `focal` scan: matrix/matrices containing* `0`,`1`, or `NA`, which represents an edge that hasn't been sampled. Internally sample from `scan$theoretical`, to which the igraph network `mode` was already applied, and is set to keep `NA`s where they were drawn (relevant for `"group"` scan sampling), even if this result in a undirected adjacency matrix where `NA`s are non-symetrical.  
        A function to minimize solvable `NA`s will be soon introduced: for `mode = "max"`, this is e.g. when `scan$raw[i,j] = 1` and `scan$raw[j,i] = NA` or inversly, for which both values can be set to `1`; for `mode = "min"`, this is e.g. when `scan$raw[i,j] = 0` and `scan$raw[j,i] = NA` or inversly, for which both values can be set to `0`.  
        *: with `mode = "plus"`, possible values also include `2`
* Revamping from the ground up as the new object were formalized and introduced in the algorithm, the internal code also got cleaned and fixed in some places as code was recycled from its previous non-OOP state:
    * Some functions were updated to split some of their overall work into additional more explicit and more "single-step/purpose" functions (e.g. `generate_obsProb` now relies on `determine_obs.prob_type` and `calculate_obs.prob` instead of containing all the code in itself)
    * `empiScan$focal` now show not only the line of the `focal` sampled, but also its column
    * rethinking of the way the igraph network `mode` is applied: `apply_mode` now keep track of an always `"directed"` `raw` scan, makes it symmetrical into a `theoretical` scan if an undirected `mode` is selected (`"undirected"`, `"max"`, `"min"`, or `"plus"`), from which a `focal` and `group` scans are sampled through `sample_from_scan` (that internally relies on `group_sample` and `focal_sample`). Through this changes, internal algorithm changed significantly, especially in how `NA`s are handled in `group` scans. Also, `zero_NA` has been generalized into `replace_NA`, with subsequent changes in `compare_with_transposed`, but these are not used anymore within `apply_mode` and are kept for now just in case.
    * removed normally-superseded `*.R` files (`focal.list.R`, `obs.prob_._tools`), renamed `do.scan.R` into `do.scan.old.R`

# SimuNet 0.1.0
* Imported simulation-oriented functions from ConfiNet
* Structure package till it checks ok on its own
* Commented and removed mention of `decide_use.rare.opti()` for now till this routine is cleaned. Updated examples in `Boot_scans()` and bootstrap tools to not call for this function through setting `use.rare.opti = FALSE`
* Added a `NEWS.md` file to track changes to the package.

