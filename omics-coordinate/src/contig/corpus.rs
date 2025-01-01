//! The global corpus of contig names.

use std::sync::LazyLock;
use std::sync::RwLock;

use string_interner::StringInterner;
use string_interner::backend::StringBackend;

/// A global corpus for contig names.
pub static CORPUS: LazyLock<RwLock<StringInterner<StringBackend>>> =
    LazyLock::new(|| RwLock::new(StringInterner::<StringBackend>::new()));
