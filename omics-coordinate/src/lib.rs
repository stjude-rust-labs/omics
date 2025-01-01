//! Coordinates upon a molecule.
//!
//! A **coordinate** is the fundamental unit for describing a location within a
//! genome. Coordinates point to a single location within a contiguous molecule
//! (typically a nucleic acid molecule, such as DNA or RNA, or a protein) and
//! are specified at the _nucleotide_ level of abstraction.
//!
//! Coordinates are comprised of three components:
//!
//! * The name of the molecule upon which the coordinate sits is known as the
//!   [**contig**](crate::Contig).
//! * Each molecule is made of a contiguous series of elements. The offset of
//!   the selected element with respect to the starting element of the molecule
//!   is known as the [**position**](crate::Position).
//! * Optionally, if the molecule is stranded, the strand upon which the
//!   coordinate sits is known as the [**strand**](crate::Strand).
//!
//! Coordinates, via their positions, can fall within the _interbase_ coordinate
//! system (which is closely related to the 0-based, half-open coordinate
//! system) or the _in-base_ coordinate system (closely related to the 1-based,
//! full-closed coordinate system). If you want to learn more about the
//! supported coordinate systems, or if you want to learn why this crate uses
//! the terms that it does (e.g., "in-base" instead of "1-based"), please jump
//! to [this section](crate#positions) of the docs.
//!
//! ### Scope
//!
//! At present, `omics-coordinate` is focused almost exclusively on nucleic acid
//! molecules. In the future, however, we expect to expand this to cover
//! proteins as well.
//!
//! ### Quickstart
//!
//! To get started, you'll need to decide if you want to use 0-based or 1-based
//! coordinates. This decision largely depends on your use case, the consumers
//! of the data, and the context of both (a) where input data is coming from and
//! (b) where output data will be shared. Note that, if you're working with a
//! common bioinformatics file format, the coordinate system is often dictated
//! by the format itself. If you need help deciding which coordinate system to
//! use, you should start by reading [the positions section](#positions) of the
//! docs.
//!
//! Once you've decided on which coordinate system you'd like to use, you can
//! create coordinates like so:
//!
//! ```
//! use omics_coordinate::Coordinate;
//! use omics_coordinate::system::Base;
//! use omics_coordinate::system::Interbase;
//!
//! // An 0-based, interbase coordinate.
//! let coordinate = Coordinate::<Interbase>::try_new("seq0", "+", 0)?;
//! println!("{:#}", coordinate);
//!
//! // A 1-based, in-base coordinate.
//! let coordinate = Coordinate::<Base>::try_new("seq0", "+", 1)?;
//! println!("{:#}", coordinate);
//!
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! For convenience, the crate also provides type aliases for the 0-based and
//! 1-based variants of the relevant concepts. For example, you can use a
//! [`Position<Interbase>`] by instead simply importing a
//! [`zero::Position`](crate::position::zero::Position).
//!
//! ```
//! use omics_coordinate::interbase::Coordinate;
//!
//! let coordinate = Coordinate::try_new("seq0", "+", 0)?;
//! println!("{:#}", coordinate);
//!
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Background
//!
//! Coordinate systems can be surprisingly hard to find comprehensive,
//! authoritative material for and, thus, have a reputation for being confusing
//! to newcomers to the field. To address this lack of material and to describe
//! how terms are used within this crate, the authors lay out their
//! understanding of the history behind the terminology used in the community
//! and then cover their perspective on what terms are most appropriate to be
//! used within different contexts. Notably, this may not match the worldview of
//! other popular resources or papers out there. In these cases, departures from
//! convention are noted alongside carefully reasoned opinions on why the
//! departure was made.
//!
//! ## Biology Primer
//!
//! Before diving into the coordinate system-specific details, we must first lay
//! some groundwork for terms used within genomics in general. These definitions
//! serve as a quick overview to orient you to the discussion around coordinate
//! systems—if you're interested in more detailed information, you can learn
//! more at [https://learngenomics.dev](https://learngenomics.dev).
//!
//! * A **genome** is the complete set of genetic code stored within a cell
//!   ([learn more](https://www.genome.gov/genetics-glossary/Genome)).
//! * **Deoxyribose nucleic acid**, or **DNA**, is a molecule that warehouses
//!   the aforementioned genetic code. In eukaryotic cells, DNA resides in the
//!   nucleus of a cell.
//!     * DNA is stored as a sequence of **nucleotides** (i.e., `A`, `C`, `G`,
//!       and `T`).
//!     * DNA is double-stranded, meaning there are two, complementary sequences
//!       of nucleotides that run in antiparallel.
//! * **Ribonucleic acid**, or **RNA**, is a molecule that is _transcribed_ from
//!   a particular stretch of DNA.
//!     * RNA is _also_ stored as sequence of nucleotides (though, in this case,
//!       the nucleotides are `A`, `C`, `G`, and `U`).
//!     * RNA is single-stranded, meaning that it represents the transcription
//!       of only one of the strands of DNA.
//!      * RNA generally either (a) serves as a template for the production of a
//!        protein or (b) has some functional role in and of itself.
//! * **Proteins** are macromolecules that are assembled by _translating_ the
//!   nucleotide sequence stored with an RNA molecule into a chain of amino
//!   acids. Proteins play a wide variety of roles in the function of a cell.
//!
//! Though there are exceptions to this rule, the core idea is this: through a
//! series of steps described within [the central dogma of molecular
//! biology](https://en.wikipedia.org/wiki/Central_dogma_of_molecular_biology),
//! genetic code stored within DNA is commonly transcribed to RNA and either (a)
//! the RNA is used as a template to assemble a functional protein through the
//! process of translation [in the case of _coding_ RNA], or (b) that RNA plays
//! some functional role in and of itself [in the case of _non-coding_ RNA].
//!
//! This crate attempts to provide facilities to effectively describe
//! coordinates within the context of DNA molecules and RNA molecules in the
//! various notations used within the community. We'll start with the most
//! granular concepts (e.g., contigs, positions, and strands) and work our way
//! up to the most broad reaching concepts (e.g., intervals and coordinate
//! systems).
//!
//! ## Contigs
//!
//! Typically, genetic information that constitutes a genome is not stored as a
//! single, contiguous molecule. Instead, genomes are commonly broken up into
//! multiple, contiguous molecules of DNA known as **chromosomes**. Beyond the
//! chromosomes, other sequences, such as the [Epstein–Barr virus][chrEBV], the
//! [mitochondrial genome][chrMT], or decoy sequences are inserted as contigs
//! within a reference genome to serve various purposes. This broader category
//! of contiguous nucleotide sequences are colloquially referred to as
//! "contigs".
//!
//! As we learn more about the human genome, new versions, called **genome
//! builds** are released that describe the known genetic sequence therein. Each
//! contigs contained within a particular genome build is assigned a unique
//! identifier within that build (e.g., `chr1` within the `hg38` genome build).
//! Specifying the contiguous molecule upon which a coordinate is located is the
//! first step in anchoring the coordinate within a genome.
//!
//! For example, the [most recent release][t2t-genome] ([ref][t2t-publication])
//! of the human genome at the time of writing has _exactly_ 24 contigs—these
//! represent the 22 autosomes and the X/Y sex chromosomes present in the human
//! genome. Interestingly, earlier versions of the human genome, such as
//! [GRCh37][grch37-genome] and [GRCh38][grch38-genome], contain more contigs
//! that represent phenomenon such as unplaced sequences (i.e., sequences that
//! we know are located _somewhere_ in the human genome, but we didn't know
//! exactly where when the reference genome was released) and unlocalized
//! sequences (i.e., sequences where we know the chromosome upon which the
//! sequence was located but not the exact position).
//!
//! #### Design Considerations
//!
//! There are no current or planned restrictions on what a contig can be named,
//! as the crate needs to remain able to support all possible use cases. That
//! said, the authors may introduce (optional) convenience methods based on
//! common naming conventions in the future, such as the detection of `chr`
//! prefixes, which is a convention for the naming of chromosomes specifically.
//!
//! ## Positions
//!
//! This section lays out a detailed, conceptual model within which we can
//! compare and contrast the two kinds of positions used within genomic
//! coordinate systems: namely, _in-base_ positions and _interbase_ positions.
//! We then cover how these terms relate to commonly used terms in the community
//! (including a "0-based, half-open coordinate system" and a "1-based,
//! fully-closed coordinate system") and how you can use this crate to flexibly
//! represent a spectrum of locations within a genome.
//!
//! Before we begin, a word of caution—many materials attempt to make the
//! differences between in-base and interbase positions (or the closely related
//! 0-based, half-open and 1-based, fully closed coordinate systems) appear
//! small and unremarkable (e.g., by providing seemingly straightforward
//! formulas to convert between the two). In fact, after a quick scan of these
//! materials, you may even be tempted to view the two systems as simply a
//! difference in accounting and off-by-one hoopla!
//!
//! In the authors' opinion, not only is this not true, it also doesn't serve
//! you well to think of the coordinate systems as anything less than entirely
//! different universes that must be responsibly traversed between. To be clear,
//! we're not suggesting that the existing materials are _wrong_—often, you can
//! follow the conventions laid out, and, as long as the baked-in assumptions
//! are consistently true for your use case, everything will be well. That said,
//! we endeavour to go futher within this crate—to explore the very fabric of
//! these coordinate systems, point out the assumptions made in each coordinate
//! system, and enable you to understand and write code that works across the
//! spectrum of possible position representations.
//!
//! #### In-base and Interbase Positions
//!
//! Positions within a genomic coordinate system can be represented as either
//! _in-base_ positions or _interbase_ positions:
//!
//! * **In-base** positions point directly to and fully encapsulate a
//!   nucleotide. These types of positions are generally considered to be
//!   intuitive from a biological reasoning standpoint and are often used in
//!   contexts where data is reported back to a biological audience (e.g.,
//!   genome browsers and public variant databases). Though we use the term
//!   "in-base" exclusively in this document, these types of positions are also
//!   sometimes referred to as simply "base" positions in the broader community.
//! * **Interbase** positions point to the spaces _between_ nucleotides. These
//!   positions are generally considered to be easier to work with
//!   computationally for a variety of reasons that will become apparent in the
//!   text that follows. It is also possible to unambiguously represent certain
//!   types of variation, such as insertions and structural variant breakpoints,
//!   using interbase positions. As such, interbase positions are commonly used
//!   as the internal representation of positions within bioinformatics tools as
//!   well as in situations where the output is meant to be consumed
//!   computationally (e.g., APIs).
//!
//! For example, SAM files, which are intended to be human-readable, use in-base
//! positions to make themselves more easily interpretable and compatible with
//! genomic databases. Their non-human-readable, binary counterparts, known as
//! BAM files, use interbase positions for the reasons describe aboved. The
//! decision on which coordinate system to use was largely based on the
//! distinction on how the two file types were meant to be consumed (to learn
//! more about what the author of SAM/BAM said about the decision, read the end
//! of [this StackExchange
//! answer](https://bioinformatics.stackexchange.com/a/17757)).
//!
//! #### Conceptual Model
//!
//! Here, we introduce a conceptual model that is useful for comparing and
//! contrasting the two coordinate systems. Under this model, nucleotides and
//! the spaces between them are pulled apart and considered to coexist as
//! independent entities laid out along a discrete axis. Both nucleotides and
//! spaces represent a "slot", and the kind of slot may be distinguished by
//! designating it as a "nucleotide slot" and a "space slot" respectively.
//! Numbered positions are assigned equidistantly at every other slot within
//! either system, but the type of slot where positions are assigned is mutually
//! exclusive between the two systems:
//!
//! * Numbered positions are assigned to each of the nucleotide slots within the
//!   in-base coordinate system.
//! * Numbered positions are assigned to each of the space slots within the
//!   interbase coordinate system.
//!
//! Importantly, in both systems, **only slots with an assigned position can be
//! specified using a position**. This has incredibly important implications on
//! what locations can and cannot be expressed within the two coordinate
//! systems.
//!
//! The diagram below depicts the model applied over a short sequence of seven
//! nucleotides. Each slot has a series of double pipe characters (`║`) that
//! links a slot with its assigned, numbered position (if it exists) within the
//! in-base and interbase coordinate systems. Note that, though the two
//! positions systems are displayed in parallel in the diagram below, that is
//! only so that they can be compared/contrasted more easily. More specifically,
//! **they do not interact with each other in any way**.
//!
//! ```text
//! ========================== seq0 =========================
//! •   G   •   A   •   T   •   A   •   T   •   G   •   A   •
//! ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║
//! ║[--1--]║[--2--]║[--3--]║[--4--]║[--5--]║[--6--]║[--7--]║ In-base Positions
//! 0       1       2       3       4       5       6       7 Interbase Positions
//! ```
//!
//! As was alluded to above, reasoning about the in-base coordinate system under
//! this model is relatively straightforward—if one wants to create a position
//! representing the location of the first nucleotide (`G`), it can be done by
//! simply denoting the numbered position assigned to same slot as the `G`
//! nucleotide, which is position `1`.
//!
//! Creating a position that represents the same nucleotide using the interbase
//! coordinate system is more complicated. Recall that (a) no numbered positions
//! are assigned to nucleotide slots within the interbase coordinate system and
//! (b) only numbered slots may be referenced as a position. As such, referring
//! to the first nucleotide using a single, numbered position is impossible.
//! Indeed, in a strict sense, a _range_ of numbered positions must be used to
//! encapsulate even this single nucleotide (🤯)—namely, the range `[0-1]` (note
//! that the range of interbase positions is generally considered _exclusive_,
//! but that does not apply here when the space slots and nucleotide slots are
//! split).
//!
//! #### Starting Position
//!
//! By convention within the community, interbase positions almost always start
//! at position zero (`0`) and in-base positions almost always start at position
//! one (`1`). As far as the authors can tell, this is for three main reasons
//! (please contribute to the docs if you disagree with any of these assertions
//! or know of other reasons):
//!
//! * **History.** Biological coordinate systems and databases have historically
//!   used a starting position of `1`. Thus, in-base coordinates (which, again,
//!   are generally considered to be more suitable for a broader biological
//!   audience) tend to follow these same conventions.
//! * **Intention.** Interbase coordinates, on the other hand, depart from a
//!   biologically intuitive model in favor of a more computationally intuitive
//!   model. To that end, interbase positions typically mirror programming
//!   languages in that counting starts at `0`. This suggests that, many times,
//!   interbase coordinates are a more natural fit for existing data structures
//!   and algorithms.
//! * **Convention.** Beyond the reasons above (and, further, not strictly
//!   imposed by the definitions of interbase and in-base coordinate systems),
//!   the community has evolved to use the starting position of `0` or `1` to
//!   allude to the use of interbase and in-base positions, respectively.
//!
//! #### Design Considerations
//!
//! The previously decsribed inability to represent interbase positions as a
//! single number presents a number of practical problems.
//!
//! For example, to accurately model positions as described above, a crate would
//! need to support both _numerical_ positions and _interval-based_ positions at
//! the same time. All higher-order concepts that include positions, such as
//! coordinates and intervals, would need to somehow present an ergonomic
//! interface and mental model for working with these very different models of a
//! position. Among other drawbacks, modeling things in this way would introduce
//! an incredible duplication of effort and additional opportunities for bugs to
//! be introduced.
//!
//! Beyond these practical considerations, designing a range-based, singular
//! position is not trivial. For example, any range must have a more
//! fundamental, singular type that represents the start and the end of the
//! range:
//!
//! * Should the crate introduce an even lower level concept into the crate
//!   below positions (e.g., a "number"?) that enables this design? If so, this
//!   many levels of abstraction introduce significant additional mental load
//!   for would-be users of such a crate.
//! * How would these range-based positions interact with the aforementioned
//!   upstream facilities, such as intervals? Intervals start and end with a
//!   position—isn't it much more confusing for users of the crate if an
//!   interval starts and ends with an even _lower level_ concept of a
//!   range/interval?
//!
//! In pursuit of pragmatism, this crate codifies the heuristic included in many
//! that precede it: interbase positions are, instead, represented as single
//! number that includes the nucleotide following the numbered space slot. This
//! allows for a much simpler and interoperable representation of positions
//! between coordinate systems, as the interbase position representing the first
//! nucleotide `G` is now simply `0` while the in-base position for the first
//! nucleotide is still `1`. Further, this assumption works nicely with the
//! expected behavior of intervals, which is discussed further in [the intervals
//! section of the docs](#intervals).
//!
//! #### Final Thoughts
//!
//! Though the authors feel it is more intuitive to teach the positioning
//! systems using the "interbase" and "in-base" nomenclature (and, explicitly,
//! we wish these designations were used more pervasively in the community!),
//! these terms are not frequently used in the literature today. Indeed, it is
//! much more common to hear interbase positions referred to as "0-based"
//! positions and in-base positions referred to as "1-based" positions.
//!
//! As such, the following statements are true throughout the rest of this
//! document and within the crate itself:
//!
//! * The term **0-based** is used in place of and is interchangeable with the
//!   term "interbase" with the codified assumption that the coordinate system
//!   will always start at position zero.
//! * The term **1-based** is used in place of and is interchangeable with the
//!   terms "in-base" and "base" with the codified assumption that the
//!   coordinate system will always start at position one.
//!
//! ## Strand
//!
//! DNA is a double-stranded molecule that stores genetic code. This means that
//! two sequences of complementary nucleotides run in antiparallel. This is
//! often referred to as being read from [5' to
//! 3'](https://en.wikipedia.org/wiki/Directionality_%28molecular_biology%29),
//! referring to connections within the underlying chemical structure. For
//! example, below is a fictional double-stranded molecule with the name `seq0`.
//!
//! ```text
//! ---------------- Read this direction --------------->
//!
//! 5'                                                 3'
//! ===================== seq0 (+) ======================
//! G   A   T   A   T   G   A   A   T   A   T   G   A   G
//! |   |   |   |   |   |   |   |   |   |   |   |   |   |
//! C   T   A   T   A   C   T   T   A   T   A   C   T   C
//! ===================== seq0 (-) ======================
//! 3'                                                 5'
//!
//! <--------------- Read this direction ----------------
//! ```
//!
//! In a real-world, biological context, both strands contain genetic
//! information that is important to the function of the cell—though both
//! strands are biologically important, _some_ system of labelling must be
//! introduced to distinguish which of the two strands a genomic coordinate is
//! located on.
//!
//! To address this, a reference genome selects one of the strands as the
//! **positive** strand (also called the "sense" strand, the "reference" strand,
//! or the `+` strand) for each contiguous molecule. This implies that the
//! opposite, complementary strand is the **negative** strand (also called the
//! "antisense" strand, the "complementary" strand, or the `-` strand). Notably,
//! reference genomes only specify the nucleotide sequence for the _positive_
//! strand, as the negative strand's nucleotide sequence may be computed as the
//! reverse complement of the positive strand.
//!
//! The concept of strandedness is useful when describing the location of
//! coordinate on a molecule with two strands. Some nucleic acid molecules, such
//! as RNA are single-stranded molecules—RNA is _derived_ from a particular
//! strand of DNA, but the RNA molecule itself is not considered to be stranded.
//!
//! Within this crate, a [`Strand`] always refers to the strand of the
//! coordinate upon a molecule (if the molecule is stranded). If the molecule
//! upon which the nucleotide(s) sit is _not_ stranded, then no strand should be
//! specified.
//!
//! This means that,
//!
//! * Coordinates that lie upon a DNA molecule must always have a strand. The
//!   [`Strand::Positive`] and [`Strand::Negative`] variants are used to
//!   distinguish which strand a coordinate sits upon relative to the strand
//!   specified in the reference genome.
//! * Coordinates that lie upon an RNA molecule have no strand. In particular,
//!   the the original strand of DNA from which a position on RNA is derived is
//!   lost during any conversion from one to the other. If it is of interest,
//!   you may keep track of this kind of thing on your own at conversion time.
//!
//! ## Intervals
//!
//! Intervals describe a range of positions upon a contiguous molecule.
//! Generally speaking, you can think of an interval as simply a start
//! coordinate and end coordinate.
//!
//! As described above, positions can be either interbase (includes the
//! nucletide following the specified numbered space slot) or in-base (includes
//! the nucleotide at the specified numbered nucleotide slot). Given these
//! characteristics, intervals that are comprised of these two different types
//! of positions generally behave differently to accentuate their strong points:
//!
//! - Interbase intervals tend to be **half-open**, meaning that all nucleotides
//!   contained between the start and end positions (but not including the last
//!   position) are included within the range.
//! - In-base intervals tend to be **fully-closed**, meaning that both the
//!   nucleotides at the start and end positions of the interval are included in
//!   the range.
//!
//! The following figure illustrates this concept using the notation described
//! in [the position section of the docs](#positions).
//!
//! ```text
//! ========================== seq0 ===========================
//! •   G   •   A   •   T   •   A   •   T   •   G   •   A   •
//! ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║
//! ║   1   ║   2   ║   3   ║   4   ║   5   ║   6   ║   7   ║   In-base Positions
//! 0       1       2       3       4       5       6       7   Interbase Positions
//! ===========================================================
//! ┃   ┃                                               ┃   ░
//! ┃   ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛   ░   seq0:+:1-7 (1-based, fully-closed)
//! ┃                Both contain "GATATGA"                 ░
//! ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━   ░   seq0:+:0-7 (0-based, half-open)
//! ```
//!
//! By looking at this figure, the reason for not including the end position in
//! the interbase coordinate system should be relatively intuitive: inclusion of
//! position seven (`7`) in the interbase interval would mean that the
//! nucleotide following position seven would also be included in the range.
//!
//! Notably, this means that intervals in the two systems need to be treated
//! carefully internally. For example, the length of an interval in the
//! interbase coordinate system is found with the formula `end - start`, while
//! the length of an interval in the in-base coordinate system is `end - start +
//! 1`. That being said, this crate largely handles the differences in
//! implementation for these two coordiante systems, meaning that you can use
//! either with confidence via a common interface.
//!
//! # Crate Design
//!
//! Throughout the crate, you will see references to 0-based and 1-based
//! variants of the concepts above. For example, there is a core [`Position`]
//! struct that is defined like so:
//!
//! ```ignore
//! pub struct Position<S>
//! where
//!     S: System, {
//!     // private fields
//! }
//! ```

// TODO: this is a false positive missing doc link, remove this when it gets fixed.
#![allow(rustdoc::broken_intra_doc_links)]
//! The struct takes a single, generic parameter that is a [`System`]. In this
//! design, functionality that is fundamental to both 0-based and 1-based
//! position types are implemented in the core [`Position`] struct.
//! Functionality that is different between the two coordinate systems is
//! implemented through traits (in the case of positions,
//! [the `Position` trait](crate::position::r#trait::Position<S>)) and exposed
//! through trait-constrained methods (e.g., [`Position::try_new`]).

//! Note that some concepts, such as [`Contig`] and [`Strand`] are coordinate
//! system invariant. As such, they don't take a [`System`] generic type
//! parameter.
//!
//! ## Learning More
//!
//! In the original writing of these docs, it was difficult to find a single,
//! authoritative source regarding all of the conventions and assumptions that
//! go into coordinate systems. Here are a few links that the authors consulted
//! when writing this crate.
//!
//! * [This blog post](https://genome-blog.gi.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/)
//!   from the UCSC genome browser team does a pretty good job explaining the
//!   basics of 0-based versus 1-based coordinate systems and why they are used
//!   in different contexts.
//!     * Note that this crate does not follow the conventions UCSC uses for
//!       formatting the two coordinate systems differently (e.g. `seq0 0 1` for
//!       0-based coordinates and `seq1:1-1`). Instead, the two coordinate
//!       systems are distinguished by the Rust type system and are serialized
//!       similarly (e.g., `seq0:+:0-1` for 0-based coordinates and `seq0:+:1-1`
//!       for 1-based coordinates).
//! * [This blog post](https://tidyomics.com/blog/2018/12/09/2018-12-09-the-devil-0-and-1-coordinate-system-in-genomics/)
//!   also presents the two coordinate systems and gives some details about
//!   concrete file formats where each are used.
//! * [This cheat sheet](https://www.biostars.org/p/84686/) is a popular
//!   community resource (though, you should be sure to read the comments!).
//!
//! [chrEBV]: https://en.wikipedia.org/wiki/Epstein%E2%80%93Barr_virus
//! [grch37-genome]: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/
//! [grch38-genome]: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/
//! [chrMT]: https://en.wikipedia.org/wiki/Mitochondrial_DNA
//! [t2t-genome]: https://www.ncbi.nlm.nih.gov/assembly/GCF_009914755.1/
//! [t2t-publication]: https://www.science.org/doi/10.1126/science.abj6987

pub mod contig;
pub mod coordinate;
pub mod interval;
pub mod math;
pub mod position;
pub mod strand;
pub mod system;

pub use contig::Contig;
pub use coordinate::Coordinate;
pub use coordinate::base;
pub use coordinate::interbase;
pub use interval::Interval;
pub use position::Position;
pub use strand::Strand;
