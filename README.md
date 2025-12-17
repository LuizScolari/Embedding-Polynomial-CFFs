# Embedding Polynomial CFFs

This project implements the generation and embedding (expansion) of **Cover-Free Families (CFFs)** using constructions based on polynomials over finite fields. The code is written in C and optimized for performance.

## 📚 What are CFFs?

**Cover-Free Families (CFFs)** are combinatorial structures represented as binary matrices. A matrix is considered a $d$-CFF if the union of any $d$ columns does not completely cover any other remaining column.

They have important applications in:

  * Group Testing
  * Cryptography and key distribution
  * Digital signatures

## 🧮 Implemented Constructions

The project focuses on algebraic construction using polynomials over finite fields ($\mathbb{F}_q$).

### Polynomial Construction (Standard)

In this construction, matrix rows are indexed by pairs of elements $(x, y) \in \mathbb{F}_q \times \mathbb{F}_q$, and columns are indexed by polynomials of degree up to $k$. A position in the matrix is set to $1$ if the corresponding polynomial evaluated at $x$ yields $y$ (i.e., $P(x) = y$), and $0$ otherwise.

### Monotone Construction

The Monotone construction is a variation designed for embedding operations. During expansion, rows are indexed by pairs $(x, y)$ where $x$ is restricted to a subset $B \subseteq \mathbb{F}_q$ corresponding to the smaller base field. The parameters $d$ and $k$ are always kept constant throughout the embedding process—only the field size $q$ changes.

## 🚀 Hash Table Optimizations

Naive generation of polynomial CFFs can be computationally expensive due to the need to evaluate many polynomials at many points.

To optimize this process, this project uses **Hash Tables (GHashTable from GLib)** to create an inverted index of evaluations:

1.  Instead of repeatedly recalculating $P(x)$, we precompute and store which polynomials yield a value $y$ for a given $x$.
2.  This transforms the search into an average constant-time $O(1)$ operation for filling the matrix, drastically speeding up the generation of large instances.

## 🛠️ Prerequisites

To compile this project, you will need the following libraries installed on your system:

  * **GCC** or **Clang**
  * **Make**
  * **FLINT 3.3.1** (Fast Library for Number Theory)
  * **GMP** (GNU Multiple Precision Arithmetic Library)
  * **GLib 2.0**
  * **OpenMP** (libomp)

### Supported Platforms

This project supports **Linux** and **macOS**. The Makefile automatically detects the operating system and configures the appropriate compiler flags.

| Platform | Compiler | Notes |
|----------|----------|-------|
| Linux | GCC | Standard configuration |
| macOS | Clang | Requires Homebrew for dependencies |

#### Installing Dependencies

**Linux (Ubuntu/Debian):**
```bash
sudo apt-get install build-essential libgmp-dev libglib2.0-dev libomp-dev
# FLINT 3.3.1 must be compiled from source (see flintlib.org)
```

**macOS (Homebrew):**
```bash
brew install gmp glib libomp flint
```

## 💻 How to Run

### 1\. Compilation

Use the provided `Makefile` to compile the project:

```bash
make
```

### 2\. Execution

The `generate_cff` executable supports different operation modes. Arguments vary according to the construction type (`p` for polynomial, `m` for monotone) and action (`f` to generate from scratch, `g` for embedding/expansion).

**General Syntax:**
`./generate_cff <type> <action> [parameters...]`

#### Examples:

  * **Generate a Polynomial CFF from scratch (`p f`):**

      * Parameters: `q` (field size), `k` (degree).

    <!-- end list -->

    ```bash
    ./generate_cff p f 3 1
    ```

  * **Polynomial Embedding (Expansion) (`p g`):**

      * Expands from a smaller field to a larger one.
      * Parameters: `initial_q` `final_q` `initial_k` `final_k`.

    <!-- end list -->

    ```bash
    ./generate_cff p g 3 9 1 1
    ```

  * **Monotone Embedding (Expansion) (`m g`):**

      * Parameters: `d` `initial_q` `final_q` `initial_k` `final_k`.

    <!-- end list -->

    ```bash
    ./generate_cff m g 2 3 9 1 1
    ```

Output files will be generated in the `CFFs/` folder.

## 📊 Benchmark

The project includes an automated benchmark system to measure the execution time of CFF generation. The benchmarks measure two metrics:

- **Time 1**: Total time for inverted index + CFF generation (+ concatenation for embedding)
- **Time 2**: Time for CFF matrix generation + concatenation only

### Quick Start

```bash
cd benchmark
make
./run_benchmarks.sh
```

For more details on how to run the benchmarks, configure iterations, and interpret results, see the [benchmark/README.md](benchmark/README.md).

-----

## 🎓 Supervision

This research is being supervised by **Professor Thaís Bardini Idalino** ([Website](https://thaisidalino.github.io), [Google Scholar](https://scholar.google.com/citations?user=hS_R7ZsAAAAJ&hl)).

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](https://www.google.com/search?q=LICENSE) file for details.