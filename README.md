# Corpus Query Parser

A high-performance C++ tool for searching and analyzing linguistic corpora using a custom query language. This parser enables complex pattern matching across tokenized text documents with support for multiple linguistic attributes.

## Overview

This project implements an efficient corpus search engine that can:
- Parse and query large text corpora (tested with 5M+ tokens)
- Search using multiple linguistic attributes (word forms, lemmas, part-of-speech tags)
- Support both equality and inequality constraints
- Display search results with highlighted matches
- Optimize queries using indexed lookups and set operations

## Features

### Query Language

Queries are written using a bracket-based syntax where each clause represents a token:

```
[attribute="value"] [attribute="value"] ...
```

#### Supported Attributes

- `word` - The surface form of the word
- `lemma` - The base/dictionary form
- `pos` - Part-of-speech tag
- `c5` - C5 tagset (a simplified POS classification)

#### Operators

- `=` - Equality (token must match)
- `!=` - Inequality (token must not match)

### Example Queries

```
[word="the"] [pos="NOUN"]          # Find "the" followed by any noun
[lemma="be"] [word="running"]      # Find forms of "be" followed by "running"
[pos!="VERB"] [pos="VERB"]         # Find non-verb followed by verb
```

## Data Format

The corpus file should be a tab-separated CSV with four columns:

```
word    c5      lemma   pos
the     DT      the     DET
cat     NN1     cat     NOUN
sat     VVD     sit     VERB
```

- Empty lines separate sentences
- Lines starting with `#` are treated as comments

## Architecture

### Core Components

#### Corpus Structure
- Tokenized representation with efficient string interning
- Sentence boundary tracking
- Multiple indices for fast attribute lookups

#### Set Operations
Three set types optimize different scenarios:
- `DenseSet` - Contiguous ranges of tokens
- `IndexSet` - Sorted arrays from index lookups
- `ExplicitSet` - Arbitrary collections after operations

#### Query Execution
1. Parse query into clauses (literals)
2. Look up each literal in the appropriate index
3. Apply set intersection/difference operations
4. Filter results to sentence boundaries
5. Display matches with highlighting

### Performance Optimizations

- **Indexed lookups**: Pre-built sorted indices for each attribute
- **Set algebra**: Efficient intersection and difference operations
- **Adaptive algorithms**: Switches between merge-based and binary search depending on set sizes
- **Early filtering**: Processes smallest sets first to minimize work

## Building

Requires C++20 or later for `std::span` support.

```bash
g++ -std=c++20 -O3 parser.cpp -o corpus_parser
```

## Usage

```bash
./corpus_parser
```

The program will:
1. Load the corpus from `bnc-05M.csv`
2. Build indices (this may take a moment for large corpora)
3. Prompt for queries
4. Display up to 10 matches per query with highlighted matching tokens
5. Continue until you press Enter on an empty query

### Output Format

Matches are displayed with full sentences, showing all four attributes in columns:

```
WORD    C5      LEMMA   POS
```

Matching tokens appear in red (ANSI color codes) to highlight query matches.

## Implementation Details

### Index Building

Each attribute gets a sorted index mapping values to token positions, enabling O(log n) lookups with binary search.

### Query Optimization

- Clauses with multiple literals are intersected
- Results are sorted by set size (smallest first)
- Dense sets are collapsed before intersection with sparse sets
- Complement sets use DeMorgan's laws for efficient computation

### Memory Management

- String interning reduces memory footprint
- Index structures use native integer types
- Spans avoid unnecessary copying

## Limitations

- Queries must match consecutive tokens within a single sentence
- No support for wildcards or regular expressions
- Inequality constraints are less efficient than equality
- Maximum 10 results displayed per query (hardcoded)



