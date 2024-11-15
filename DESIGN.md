- [ ] All impl-specific functions pass around `Number`s to disambiguate 
      things
      for the compiler. Traits like `From` and `TryFrom` tend to use generic bounds
      to make things maximally compatible.
  - TLDR: use the functions like `new()`/`try_new()` instead of
    `from()`/`try_from()` as a user, lest you suffer from constantly
    disambiguating your types.

## Strand

- [ ] `Positive` and `Negative` were chosen because
  - `Sense` and `Antisense` are in relation to protein translation exclusively
    and
  - `Forward` is wrapped up in "moving forward" and "moving backward".

## Interval

- Why not use std::ops::RangeBounds etc.
  - Because they assume monotonically increasing.
