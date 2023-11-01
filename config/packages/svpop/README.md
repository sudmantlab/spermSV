# SL's guide to the svpop cinematic universe
I parsed all three documentation markdown files so you don't have to. Thank me later.

# Getting your bearings
A **merge** combines variant files, where redundant variants are collapsed (along specific parameters that determine redundancy, of course). Merges can be performed across samples (yielding one variant file per caller), or across callers (yielding one variant file per sample).

An **intersection** is a subtype of a merge, but more flexible for sample-against-sample comparison. 
For example, you can compare...
* HG001 called with `pbsv` against HG001 called with `sniffles2`.
* HG001 called with `pbsv` against HG002 called with `pbsv`.
* HG001 called with `pbsv` against HG002 called with `sniffles2`.
And so on, and so forth.

# Constructing your merge 
Merges are described by providing a key-value pair to the configfile. The key is a string describing the types of variants to merge, and the value is a string describing how you want to conduct your merge. 

We'll start with constructing the value string, because that's the hard part. 

## Merge strategies
This can get VERY COMPLEX depending on how particular you want to be about how variants are merged. It's basically regex but for merging.

At the file level, all variant files are combined iteratively. An arbitrary first variant file is chosen, and variant entries are copied into a "reference" file. The next variant file is intersected with this "reference". Variant entries in this file that already exist in the reference are noted as annotations on the reference variant entries. Variants that are *not* present in the reference as added as new variant entries.

**Step 1**. You must pick one of two "prefixes" to use in the value to define the file-level combining strategy. This prefix only exists because of the vague idea that someday there could be non-iterative file-level combining strategies. 
    * `nr`: Short for "non-redundant merge". Does the iterative merging described above.
    * `nrsnv`: Short for "non-redundant SNV merge". Same as above, but specifically enforces REF/ALT matching. (Not sure why this has to be enforced specifically for SNVs, but I assume it's a codebase thing.)

**Step 2**: Now you get to define your merge strategy customizing five* different options. Remember that you can define different merge strategies for different types of SVs, so don't try to construct a mega-string right now.

Each option has its own parameters, which can either be defined positionally (in order that they appear in the parameters table) or by identifier (any order you want, as long as you provide the parameter ID).

---

`exact`: Size and position of the variant must match exactly.

Parameters:
None except `match`.

Examples:
`nr::exact`: Matches all variants by size and position.
---

`ro`: Short for "reciprocal overlap". Variants must overlap reciprocally by size and position. This is a classic merge strategy that works well for very large variants. 

For insertions, the length is added to the start position to obtain the end position for this strategy.

Note that for insertions, the length of the insertion is added to the single-point start position in the reference in order to represent the "span" of the SV in an equivalent form to deletions.
| Parameter |         Range        |  Default  |                Description                 |
| --------- | -------------------- | --------- | ------------------------------------------ |
| ro        |    0.0 < ro ≤ 1.0    |    0.5    | Proportion of overlap (defaults to 50% RO) |
| dist      | 0 ≤ dist ≤ unlimited | unlimited | Breakpoints must be within dist bp         |

Examples:
`nr::ro(0.8)`: The variants must be at least 80% overlapping in size (`ro` parameter).
---

`szro`: Short for "size-reciprocal overlap". It's like `ro`, except the size parameter is flexible around a set distance offset. This strategy works better for smaller variants where small changes is position create large differences by ro.

The offset can be restricted by either a set number of bases (dist) or as a multiplier of the variant size (szro), or even or both. For obvious reasons, you cannot set both dist and szdist to "unlimited".
| Parameter |         Range        |  Default  |                Description                         |
| --------- | -------------------- | --------- | -------------------------------------------------- |
| szro      |   0.0 < szro ≤ 1.0   |    0.5    | Proportion of overlap (defaults to 50% size match) |
| dist      | 0 ≤ dist ≤ unlimited |    200    | Breakpoints must be within dist bp                 |
| szdist    | 0 ≤ dist ≤ unlimited | unlimited | Breakpoints must be offset ≤ szdist * size         |

Examples:
`szro(0.8,,4)`: The variants must be at least 80% overlapping in size (`szro` parameter). `szdist` is used to further specify that the breakpoints must not deviate more than 4 times the size of their overlap. `dist` is left as the default value. 
`szro(0.2, szdist=2)`: The variants must be at least 20% overlapping in size, and the breakpoints must not deviate more than 2 times the size of their overlap.

---

`distance`: A strict distance match that ignores the variant size. This strategy is useful for comparing callsets with extreme variant size inaccuracies, typcially legacy short-read SV callsets (modern short-read callsets should have more accurate sizes).

There are technically three adjustable parameters here, but according to the docs, if you need to use any parameter other than `dist`, you're better off using `szro` as your strategy. (So why enable having three params? I can only imagine the hyper-specificity of that use case...)
| Parameter |         Range        |  Default  |                Description                         |
| --------- | -------------------- | --------- | -------------------------------------------------- |
| szro      |   0.0 < szro ≤ 1.0   |     NA    | Proportion of overlap (defaults to 50% size match) |
| dist      | 0 ≤ dist ≤ unlimited |    200    | Breakpoints must be within dist bp                 |
| szdist    | 0 ≤ dist ≤ unlimited | unlimited | Breakpoints must be offset ≤ szdist * size         |

Examples:
`distance(dist=500)`: Specifies that the breakpoints must be within a 500 bp range.

---

`match`: Use this to specify sequence similarity as its own standalone strategy. 

Variant sequences can be compared by alignment similarity (Smith-Waterman alignment, see `MERGE.md` for the full details) or Jaccard distance. Variants above a certain size limit are compared by Jaccard distance to prevent long CPU time. If the limit is set to unlimited, all variants are compared by alignment.

| Parameter |         Range        |  Default  |                Description                         |
| --------- | -------------------- | --------- | -------------------------------------------------- |
| score      |   0.0 < score ≤ 1.0   |    0.8    | Minimum match score by alignment/Jaccard |
| match      | 0.0 < match |    2.0    | S-W alignment base match score                 |
| mismatch    | mismatch < 0.0 | -1.0 | S-W alignment base mismatch score         |
| open    | open < 0.0 | -1.0 | S-W alignment gap opening score         |
| extend    | extend ≤ 0.0 | -0.25 | S-W alignment gap extension score         |
| limit    | 0 ≤ limit ≤ unlimited | 4000 | Use Jaccard when variant size > limit         |
| ksize    | 0 < ksize | 9 | The k-mer size for computing Jaccard distance.         |

---

**Step 3:** Create aliases for your merge strategies

If you can't tell, these value strings can get really long and really complicated. To avoid having to type the long value strings more than once, you want to "register" them in a `merge_def` entry to `config/config.json`. **This is MANDATORY for any merges you want to perform when intersecting files, because the alias is part of the intersect outfile path.**

Below, each entry is denoted with a key (alias) and value (definition for the merge strategy).

```
"merge_def": {
    "szro-80": "nr::szro(0.8,,4):match(0.8)",
    "snv-exact": "nrsnv::exact",
    "ro-80-nm": "nr::ro(0.8)"
}
```
**Step 4:** Assign the merge to a variant type 
Finally, we come back to the key. The key is a string describing the types of variants to merge. 
    * SVs: Make sure to prefix with `sv`, as in `"sv:ins, del"`.
    * Not SVs: Can be entered as keys directly, like `snv` and `indel`.

Now we can put it all together: 
* `"sv:ins,del": "nr::szro(0.8,,4):match(0.8)"` and `"sv:inv": "nr::ro(0.8)"` are two key-value pairs describing different merge parameters for insertions + deletions versus inversions. 
* `"snv": "nrsnv::exact"` specifies a merge specific to SNVs.

It's more concise to use aliases for the value strings:
* `"sv:ins,del": "szro-80"` and `"sv:inv": "ro-80-nm"`
* `"snv": "snv-exact"`

---

*(Optional)*: Make your merge really complicated (NOT RECOMMENDED IF AT ALL POSSIBLE)

You can "chain" the above options to add multi-level complexity to your merge. You can also provide an additional `match` parameter to each option, in case you only want sequence matching to play a role in some of your merge strategies and not others. 

`nr::ro(0.8,match(0.75)):szro:match`: Variants must be at least 80% overlapping in size, and must also be at least 75% similar in sequence. `ro`-merged variants are subject to a `szro` merge with default settings. The surviving variants are then subject to a final `match` merge with default settings, which overrides the original `0.75` match score from the `ro` merge with the default `0.8` match score.

This is a good example of "you can, but you should stop to think about whether or not you *should*".

---
---

# Deciding your merge "direction"
You can combine variant files by samples or caller, depending on your desired output.

* `sampleset`: Merge across samples. Pick this if you want to have one* file per caller that provides calls for all samples.
    * Example: You used `sniffles2` and `pbsv` on HG001, HG002, and HG003. You receive a `sniffles2` and `pbsv` outfile* showing (merged) calls across HG001, HG002, and HG003. Each variant entry is unique and provides annotations if the variant is observed in multiple samples.
* `callerset`: Merge across callers. Pick this if you want to have one* file per sample that provides calls from all callers.
    * Example: You used `sniffles2` and `pbsv` on HG001, HG002, and HG003. You receive HG001, HG002, and HG003 outfiles* showing (merged) calls from `sniffles2` and `pbsv`. Each variant entry is unique and provides annotations if the variant was identified by multiple callers.

Each combination method requires specific parameters, which are detailed below.

(*) - This is an oversimplification that will become clear once you read the Snakemake command(s) below.

## `sampleset`

### Sample group `config` entry
`samplelist` is an entry to `config/config.json` that is used to define groups of samples that you want to merge with `sampleset`. It's also just plain useful for listing groups of samples by provenance (HGSVC, 1kGP, PANPAN, etc...)

```
"samplelist": {
    "hgsvc2": [
        "HG00731", "HG00732",
        "HG00512", "HG00513",
        "NA19238", "NA19239",
        "NA24385", "HG03125",
        "NA12878", "HG03486", "HG02818",
        "HG03065", "HG03683", "HG02011",
        "HG03371", "NA12329", "HG00171",
        "NA18939", "HG03732", "HG00096",
        "NA20847", "HG03009", "NA20509",
        "HG00864", "HG01505", "NA18534",
        "NA19650", "HG02587", "HG01596",
        "HG01114", "NA19983", "HG02492",
        "HG00733", "HG00514", "NA19240"
    ],
    "PUR": ["HG00733", "HG00732", "HG00731"],
    "CHS": ["HG00514", "HG00513", "HG00512"],
    "YRI": ["NA19240", "NA19238", "NA19239"]
}
```

Note that samplesets are not mutually exclusive in membership: For example, `PUR` is a subset of `hgsvc2`.


### `sampleset` config entry

You must provide one configfile entry for each `sampleset` merge you want to perform.

| Parameter |                Description                         |
| --------- | -------------------------------------------------- |
| sourcetype | Set this to `"caller"`. |
| sourcename | Set this to the caller identifier used in the `NAME` field of `config/samples.tsv`. |
| merge    | Provide the merge key-value pair(s). |
| name    | Provide the name you want to appear in figures generated from the caller outfiles. |
| description    | Write your own relevant notes here. |

The entry below specifies a `sampleset` config entry named `pavmerge` that will generate files with merged calls from samples with `pav-hifi` in the `NAME` field of `config/samples.tsv`.

```
"sampleset": {
    "pavmerge": {
        "sourcetype": "caller",
        "sourcename": "pav-hifi",
        "merge": {
            "sv:ins,del": "nr::szro(0.8,,4):match(0.8)",
            "sv:inv": "nr::ro(0.8)",
            "indel": "nr::szro(0.8,,4):match(0.8)",
            "snv": "nrsnv::exact"
        },
        "name": "PAV HiFi",
        "description": "Pav HiFi"
    }
}
```

Note that it **does not** specify what samples are input to the merge, nor what types of variants you want to see in the outfile, because apparently you have to do that by defining an output for snakemake (see "Snakemake" section).

## `callerset`

| Parameter |                Description                         |
| --------- | -------------------------------------------------- |
| callsets | A three-item list (see below for more details). |
| merge    | Provide the merge key-value pair(s). |
| name    | Provide the name you want to appear in figures generated from the caller outfiles. |
| description    | Write your own relevant notes here. |

The `callsets` parameter: A list summarizing, in order, the following parameters:
| Parameter |                Description                         |
| --------- | -------------------------------------------------- |
| sourcetype | Set this to `"caller"`. |
| sourcename | Set this to the caller identifier used in the `NAME` field of `config/samples.tsv`. |
| sourcename** | This can be identical to the previous sourcename. See the snakemake example. |

### `callerset` config entry

The entry below specifies a `callerset` configuration named `"longreads"` that will generate sample outfiles containing information from each caller.

```
"callerset": {
    "longreads": {
        "callsets": [
            ["caller", "pav-hifi", "PAVHIFI"],
            ["caller", "pav-clr", "PAVCLR"],
            ["caller", "pbsv-hifi", "PBSVHIFI"],
            ["caller", "pbsv-clr", "PBSVCLR"]
        ],
        "merge": {
            "sv:ins,del": "nr::szro(0.8,,4):match(0.8)",
            "sv:inv": "nr::ro(0.8)",
            "indel": "nr::szro(0.8,,4):match(0.8)",
            "snv": "nrsnv::exact"
        },
        "name": "HGSVC2 LR",
        "description": "PAV/PBSV (HiFi & CLR)"
    }
}
```

Note that it **does not** specify what samples are input to the merge, nor what types of variants you want to see in the outfile, because apparently you have to do that by defining an output for snakemake (see "Snakemake" section).

# Snakemake specification

Finally, after all this, you want to tell `svpop` to create specific merges by passing a desired outfile name to `snakemake`. As if all the config parameters weren't enough, we have to further direct `svpop` using the outfile name.

Let's pretend that we want to receive a file that summarizes insertion calls made by PAV across all HGSVC2 samples, using the `sampleset` merge defined in the `pavmerge` config above.

```
# the desired output
output/svpop/sampleset/pavmerge/hgsvc2/all/all/bed/sv_ins.bed.gz

# the wildcard structure, showing how snakemake will pull wildcards from the desired outfile name
output/svpop/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz
```

Parsed wildcards:
* `{sourcetype}`: `sampleset`, because we want to perform a `sampleset` merge.
* `{sourcename}`: `pavmerge`, because we want to use the `pavmerge` config for merging. Disregard the fact that `sourcetype` is a param inside the `pavmerge` config, double-naming is apparently a common convention in the ******* lab's code.
* `{sample}`: `hgsvc2`, because we want to use all samples in the `hgsvc2` sample group under the `samplelist` entry to `config/config.json`. 
* `{filter}`: Set to `'all'` because we're not filtering out regions of the genome before merging variants. (There are "built-in" filters for hg38 described in the main svpop README).
* `{svset}`: Set to `'all'` because we're not subsetting variant types to call. (Again, the main README has more info).
* `{vartype}`: `sv`, because we're looking to merge insertion calls. Can be one of (snv, sv, indel, dup, rgn, or sub).
* `{svtype}`: `ins`, because we're looking to merge insertion calls. Can be one of (snv, sv, indel, dup, rgn, or sub), according to the README, but I think this must be incorrect given that snv for `{vartype}` would preclude indel.

Next, let's say that we want to create files for each sample in HGSVC2 that have merged insertion calls made by `pav` and `pbsv` for both types of PacBio data available (CLR and HiFi). We use the `callerset` merge according to the `longreads` config above:

```
# the desired outputs
output/svpop/callerset/PAVHIFI/hgsvc2/all/all/bed/sv_ins.bed.gz
output/svpop/callerset/PAVCLR/hgsvc2/all/all/bed/sv_ins.bed.gz
output/svpop/callerset/PBSVHIFI/hgsvc2/all/all/bed/sv_ins.bed.gz
output/svpop/callerset/PBSVCLR/hgsvc2/all/all/bed/sv_ins.bed.gz

# the wildcard structure, showing how snakemake will pull wildcards from the desired outfile name
output/svpop/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz
```

Parsed wildcards:
* `{sourcetype}`: `callerset`, because we want to perform a `sampleset` merge.
* `{sourcename}`: This is what is marked as `sourcename**` in the params. The redundant wildcard names come back to haunt us.
* `{sample}`: `hgsvc2`, because we want to use all samples in the `hgsvc2` sample group under the `samplelist` entry to `config/config.json`. 
* `{filter}`: Set to `'all'` because we're not filtering out regions of the genome before merging variants. (There are "built-in" filters for hg38 described in the main svpop README).
* `{svset}`: Set to `'all'` because we're not subsetting variant types to call. (Again, the main README has more info).
* `{vartype}`: `sv`, because we're looking to merge insertion calls. Can be one of (snv, sv, indel, dup, rgn, or sub).
* `{svtype}`: `ins`, because we're looking to merge insertion calls. Can be one of (snv, sv, indel, dup, rgn, or sub), according to the README, but I think this must be incorrect given that snv for `{vartype}` would preclude indel.

Finally, intersections. God save us.

```

# the wildcard structure
results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/{vartype}_{svtype}/intersect.tsv.gz
```