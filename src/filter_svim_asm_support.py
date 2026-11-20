#!/usr/bin/env python3
"""Retain SURVIVOR SV calls supported by SVIM-asm and >=2 read-based callers.

The position of each caller in INFO/SUPP_VEC is given by the SURVIVOR
input filelist, not by the lexical order of caller names.
"""

import argparse
import collections
import os
from pathlib import Path
import re
import sys
import tempfile

SUPP_VEC_HEADER = re.compile(r"^##INFO=<ID=SUPP_VEC(?:,|>)")


def read_filelist(path, required_vcf):
    with open(path, encoding="utf-8") as stream:
        filenames = [line.strip() for line in stream if line.strip()]
    if len(filenames) != 4:
        raise ValueError(
            f"Expected exactly four caller VCFs in {path}, found {len(filenames)}"
        )
    normalized = [os.path.realpath(name) for name in filenames]
    if len(set(normalized)) != len(normalized):
        raise ValueError("Caller VCF filelist contains duplicate paths")
    mandatory = os.path.realpath(required_vcf)
    if mandatory not in normalized:
        raise ValueError(f"Required SVIM-asm VCF absent from filelist: {required_vcf}")
    return normalized.index(mandatory), len(normalized)


def support_vector(info, n_callers, line_number):
    values = [token[9:] for token in info.split(";") if token.startswith("SUPP_VEC=")]
    if len(values) != 1:
        raise ValueError(f"VCF line {line_number}: missing/duplicate INFO/SUPP_VEC")
    value = values[0]
    if len(value) != n_callers or any(bit not in "01" for bit in value):
        raise ValueError(
            f"VCF line {line_number}: invalid SUPP_VEC={value!r}; "
            f"expected {n_callers} binary digits"
        )
    return value


def filter_vcf(source, dest, audit, asm_index, n_callers, min_other):
    output = Path(dest)
    output.parent.mkdir(parents=True, exist_ok=True)
    if os.path.realpath(source) == os.path.realpath(output):
        raise ValueError("Input and output VCF paths must differ")
    if not (1 <= min_other <= n_callers - 1):
        raise ValueError("min_other_callers is outside the valid range")

    counters = collections.Counter()
    vector_counts = collections.Counter()
    info_header = False
    column_header = False
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="", dir=output.parent,
            prefix=output.name + ".", suffix=".tmp", delete=False,
        ) as temp, open(source, encoding="utf-8") as handle:
            tmp_path = temp.name
            for line_num, line in enumerate(handle, 1):
                if line.startswith("##"):
                    info_header |= bool(SUPP_VEC_HEADER.match(line))
                    temp.write(line)
                    continue
                if line.startswith("#CHROM"):
                    if not info_header:
                        raise ValueError("VCF header does not declare INFO/SUPP_VEC")
                    column_header = True
                    temp.write(line)
                    continue
                if line.startswith("#"):
                    temp.write(line)
                    continue
                if not column_header:
                    raise ValueError(f"VCF line {line_num}: missing #CHROM header")
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 8:
                    raise ValueError(f"VCF line {line_num}: expected >=8 columns")
                vec = support_vector(fields[7], n_callers, line_num)
                vector_counts[vec] += 1
                counters["input_records"] += 1
                if vec[asm_index] != "1":
                    counters["excluded_no_svim_asm"] += 1
                    continue
                counters["svim_asm_supported_records"] += 1
                read_support = sum(
                    vec[i] == "1" for i in range(n_callers) if i != asm_index
                )
                if read_support < min_other:
                    counters["excluded_insufficient_other_callers"] += 1
                    continue
                counters["retained_records"] += 1
                temp.write(line)
            if not column_header:
                raise ValueError("Input VCF lacks a #CHROM column header")
        os.replace(tmp_path, output)
        tmp_path = None
        if audit:
            with open(audit, "w", encoding="utf-8") as handle:
                handle.write("metric\tcount\n")
                for name in (
                    "input_records",
                    "svim_asm_supported_records",
                    "excluded_no_svim_asm",
                    "excluded_insufficient_other_callers",
                    "retained_records",
                ):
                    handle.write(f"{name}\t{counters[name]}\n")
                for vector in sorted(vector_counts):
                    handle.write(f"support_vector_{vector}\t{vector_counts[vector]}\n")
        return counters["retained_records"]
    finally:
        if tmp_path is not None and os.path.exists(tmp_path):
            os.unlink(tmp_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="SURVIVOR merged VCF")
    parser.add_argument("--output", required=True, help="filtered per-sample VCF")
    parser.add_argument("--filelist", required=True, help="SURVIVOR input filelist")
    parser.add_argument("--required-vcf", required=True, help="SVIM-asm caller VCF")
    parser.add_argument("--min-other-callers", type=int, default=2)
    parser.add_argument("--audit", help="write selection counts in TSV format")
    args = parser.parse_args()
    try:
        index, n_callers = read_filelist(args.filelist, args.required_vcf)
        selected = filter_vcf(
            args.input, args.output, args.audit, index, n_callers,
            args.min_other_callers,
        )
    except (ValueError, OSError) as err:
        parser.exit(1, f"ERROR: {err}\n")
    print(
        f"Selected {selected} SVs supported by SVIM-asm and "
        f">={args.min_other_callers} of the other {n_callers - 1} callers",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
