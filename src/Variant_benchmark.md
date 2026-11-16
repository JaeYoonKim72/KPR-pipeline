# 5.Variant_benchmark.sh

Benchmark SNVs/indels/SVs against suitable truth sets and confident reference intervals.

Part of the [KPR-pipeline](../README.md).

**Status: placeholder only.** This file is reserved for a future analysis.
Running `all` or a step currently exits with `NOT IMPLEMENTED`; it does
not calculate any statistic.

---

## 1. Set up

### 1.1 Software

To be specified when the analysis is implemented.

### 1.2 Input and reference files

To be specified; no expected input format has been finalized.

### 1.3 Environment

~~~bash
export REF_DIR=/path/to/reference
export OUT_DIR=/path/to/results
export THREADS=8
~~~

---

## 2. Run

~~~bash
bash 5.Variant_benchmark.sh --help
~~~

`bash 5.Variant_benchmark.sh all` is **not implemented**.

---

## 3. Steps

| Step | Description | Status |
|---|---|---|
| `prepare_input` | Input validation and preparation | Not implemented |
| `analysis` | Variant benchmarking | Not implemented |
| `summary` | Output summary | Not implemented |
| `all` | All steps in order | Not implemented |

---

## 4. Configuration

| Variable | Default | Description |
|---|---|---|
| `REF_DIR` | unset | Reference resources (future) |
| `OUT_DIR` | `./results` | Output root (future) |
| `THREADS` | `8` | Processing threads (future) |

---

## 5. Output

No outputs are created. File formats and locations will be documented
when the pipeline is implemented.

---

## Software

To be determined.