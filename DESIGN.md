This file currently holds an (incomplete) list of my thoughts regard design
notes I'd like to consider fleshing out. As is indicated above, it's not
complete list, nor is it etched in stone—items may be removed at any time as the
crate is being written.

I don't think there will be much value in browsing these quick-form thoughts
besides myself, but feel free if you're interested.

# General

- [ ] All non-trait functions use around `position::Number`s explicitly to make writing
      code easier. Briefly, using `impl Into<position::Number>>` everywhere would
      mean each invocation needs to be disambiguited by the compiler. This slows
      down code authoring as well as introduces churn into the writing process.
  - On the other hand, traits like `From` and `TryFrom` tend to use more generic bounds
    to make things maximally compatible for implicit behavior.
  - TLDR: it's intended for user's to generally use the non-trait functions like
    `new()`/`try_new()` instead of `from()`/`try_from()` , lest you suffer from
    constantly disambiguating your types. That being said, the traits may still
    be implemented for compatibility reasons.

## Strand

- [ ] The two strands in double-stranded molecules could have been named many
      things: `Positive`/`Negative`, `Sense`/`Antisense`, `Forward`/`Reverse`. Here,
      we chose the names `Positive` and `Negative` because they appeared to be the
      most widely applicable. In particular,
  - `Sense` and `Antisense` were considered to be too context specific,
    as they are typically referring to protein translation within a particular
    gene or locus.
  - `Forward` and `Reverse` implies a particular directionality or momentum.
    Since this library attempts to represent coordinates on both strands with
    equal weight (only really putting focus on the reference genome's selected
    strand), these words didn't seem to fit well.

## Interval

- [ ] A fair question to ask is "why create interval classes yourself"? Why not
      just use some existing facilities in the standard library for representing
      ranges?
  - First, at the time of writing, `std::ops::RangeBounds` assumes that that a
    range is monotonically increasing—something that is explicitly not true in
    this context on the negative strand.
  - Next, the overlap between (a) operations that are useful when examining
    the Rust standard library's concept of a range and (b) genomics ranges are
    quite small. This is particularly true when considering the different
    contexts of in-base and interbase coordinates. In our consideration,
    reimplementing the few methods that _did_ overlap seemed to be worth the
    cost of avoiding confusion that might arise by confounding
    in-base/interbase coordinates with the Rust standard library's ranges.
