#!/usr/bin/env bash
set -euo pipefail

INPUT=${1:?usage: bash decon_split.sh <input.vcf.gz> [sv_length_threshold]}
MINSV=${2:-50}

NAME=$(basename "${INPUT}")
DIR=$(dirname "${INPUT}")
OUT="${DIR}/${NAME%.vcf.gz}"

# Include all target records, or emit a header-only VCF if the selection is empty.
select_targets() {
  local SOURCE=$1 DEST=$2 TARGETS
  TARGETS=$(mktemp "${OUT}.targets.XXXXXX")
  cat > "${TARGETS}"
  if [[ -s "${TARGETS}" ]]; then
    bcftools view -T <(awk 'BEGIN{OFS="\t"} NF>=2 {print $1,$2}' "${TARGETS}") \
      "${SOURCE}" -Oz -o "${DEST}"
  else
    bcftools view -h "${SOURCE}" | bgzip -c > "${DEST}"
  fi
  rm -f "${TARGETS}"
}

# Allele lengths
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${INPUT}" \
  | awk '{
      split($4, alts, ",")
      for (i in alts) {
        diff = length(alts[i]) - length($3)
        if (diff < 0) diff = -diff
        print $1, $2, $3, alts[i], diff
      }
    }' > "${OUT}.var_sizes.txt"

# Size classification
awk -v m="${MINSV}" '$5 > m {print $1, $2}' "${OUT}.var_sizes.txt" \
  | sort -u -k1,1V -k2,2n > "${OUT}.var_sizes.SVs_list.txt"

awk -v m="${MINSV}" 'NR==FNR {exclude[$1,$2]; next}
     $5 <= m && !(($1,$2) in exclude) {print $1, $2}' \
  "${OUT}.var_sizes.SVs_list.txt" "${OUT}.var_sizes.txt" \
  | sort -u -k1,1V -k2,2n > "${OUT}.var_sizes.Smalls_list.txt"

select_targets "${INPUT}" "${OUT}.var_sizes.SVs_list.vcf.gz" \
  < "${OUT}.var_sizes.SVs_list.txt"

select_targets "${INPUT}" "${OUT}.var_sizes.Smalls_list.vcf.gz" \
  < "${OUT}.var_sizes.Smalls_list.txt"

# Structural variants
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${OUT}.var_sizes.SVs_list.vcf.gz" \
  | awk -v m="${MINSV}" '{
      split($4, alts, ","); all_large = 1
      for (i in alts)
        if ((length(alts[i]) - length($3)) <= m && (length($3) - length(alts[i])) <= m) {
          all_large = 0; break
        }
      if (all_large) print $1 "\t" $2
    }' | sort -u -k1,1V -k2,2n > "${OUT}.var_sizes.SVs_list.all_50_alts.txt"

select_targets "${OUT}.var_sizes.SVs_list.vcf.gz" \
  "${OUT}.var_sizes.SVs_list.all_50_alts.vcf.gz" \
  < "${OUT}.var_sizes.SVs_list.all_50_alts.txt"

# sites that do not mix substitutions with length-changing alleles
bcftools view -H "${OUT}.var_sizes.SVs_list.all_50_alts.vcf.gz" \
  | awk -F'\t' '{
      split($5, alts, ","); has_snp = 0; has_indel = 0
      for (i in alts) {
        if (length(alts[i]) == length($4)) has_snp = 1
        if (length(alts[i]) != length($4)) has_indel = 1
      }
      if (!(has_snp && has_indel)) print
    }' \
  | select_targets "${OUT}.var_sizes.SVs_list.all_50_alts.vcf.gz" "${OUT}.var_sizes.SVs_list.all_50_alts.no_mixed.vcf.gz"

bcftools view -H "${OUT}.var_sizes.SVs_list.all_50_alts.vcf.gz" \
  | awk -F'\t' '{
      split($5, alts, ","); has_snp = 0; has_indel = 0
      for (i in alts) {
        if (length(alts[i]) == length($4)) has_snp = 1
        if (length(alts[i]) != length($4)) has_indel = 1
      }
      if (has_snp && has_indel) print
    }' \
  | select_targets "${OUT}.var_sizes.SVs_list.all_50_alts.vcf.gz" "${OUT}.var_sizes.SVs_list.all_50_alts.only_mixed.vcf.gz"

# structural deletions
bcftools view -H "${OUT}.var_sizes.SVs_list.all_50_alts.no_mixed.vcf.gz" \
  | awk -F'\t' '{
      split($5, alts, ","); all_deletion = 1
      for (i in alts) if (length(alts[i]) >= length($4)) { all_deletion = 0; break }
      if (all_deletion) print
    }' \
  | select_targets "${OUT}.var_sizes.SVs_list.all_50_alts.no_mixed.vcf.gz" "${OUT}.var_sizes.SVs_list.all_50_alts.no_mixed.DEL.vcf.gz"

bcftools filter -i 'strlen(REF) > strlen(ALT)' \
  "${OUT}.var_sizes.SVs_list.all_50_alts.no_mixed.DEL.vcf.gz" -Oz \
  | bcftools norm --rm-dup all -Oz \
    -o "${OUT}.var_sizes.SVs_list.all_50_alts.no_mixed.DEL.DEL.vcf.gz"

# structural insertions
bcftools view -H "${OUT}.var_sizes.SVs_list.all_50_alts.no_mixed.vcf.gz" \
  | awk -F'\t' '{
      split($5, alts, ","); all_insertion = 1
      for (i in alts) if (length(alts[i]) <= length($4)) { all_insertion = 0; break }
      if (all_insertion) print
    }' \
  | select_targets "${OUT}.var_sizes.SVs_list.all_50_alts.no_mixed.vcf.gz" "${OUT}.var_sizes.SVs_list.all_50_alts.no_mixed.INS.vcf.gz"

bcftools filter -i 'strlen(ALT) > strlen(REF)' \
  "${OUT}.var_sizes.SVs_list.all_50_alts.no_mixed.INS.vcf.gz" -Oz \
  | bcftools norm --rm-dup all -Oz \
    -o "${OUT}.var_sizes.SVs_list.all_50_alts.no_mixed.INS.INS.vcf.gz"

# 4. Small variants: deletions, insertions and substitutions
bcftools view -H "${OUT}.var_sizes.Smalls_list.vcf.gz" \
  | awk -F'\t' '{
      split($5, alts, ","); has_snp = 0; has_indel = 0
      for (i in alts) {
        if (length(alts[i]) == length($4)) has_snp = 1
        if (length(alts[i]) != length($4)) has_indel = 1
      }
      if (!(has_snp && has_indel)) print
    }' \
  | select_targets "${OUT}.var_sizes.Smalls_list.vcf.gz" "${OUT}.var_sizes.Smalls_list.no_mixed.vcf.gz"

bcftools view -H "${OUT}.var_sizes.Smalls_list.vcf.gz" \
  | awk -F'\t' '{
      split($5, alts, ","); has_snp = 0; has_indel = 0
      for (i in alts) {
        if (length(alts[i]) == length($4)) has_snp = 1
        if (length(alts[i]) != length($4)) has_indel = 1
      }
      if (has_snp && has_indel) print
    }' \
  | select_targets "${OUT}.var_sizes.Smalls_list.vcf.gz" "${OUT}.var_sizes.Smalls_list.only_mixed.vcf.gz"

bcftools view -H "${OUT}.var_sizes.Smalls_list.no_mixed.vcf.gz" \
  | awk -F'\t' '{
      split($5, alts, ","); all_deletion = 1
      for (i in alts) if (length(alts[i]) >= length($4)) { all_deletion = 0; break }
      if (all_deletion) print
    }' \
  | select_targets "${OUT}.var_sizes.Smalls_list.no_mixed.vcf.gz" "${OUT}.var_sizes.Smalls_list.no_mixed.DEL.vcf.gz"

bcftools filter -i 'strlen(REF) > strlen(ALT)' \
  "${OUT}.var_sizes.Smalls_list.no_mixed.DEL.vcf.gz" -Oz \
  | bcftools norm --rm-dup all -Oz \
    -o "${OUT}.var_sizes.Smalls_list.no_mixed.DEL.DEL.vcf.gz"

bcftools view -H "${OUT}.var_sizes.Smalls_list.no_mixed.vcf.gz" \
  | awk -F'\t' '{
      split($5, alts, ","); all_insertion = 1
      for (i in alts) if (length(alts[i]) <= length($4)) { all_insertion = 0; break }
      if (all_insertion) print
    }' \
  | select_targets "${OUT}.var_sizes.Smalls_list.no_mixed.vcf.gz" "${OUT}.var_sizes.Smalls_list.no_mixed.INS.vcf.gz"

bcftools filter -i 'strlen(ALT) > strlen(REF)' \
  "${OUT}.var_sizes.Smalls_list.no_mixed.INS.vcf.gz" -Oz \
  | bcftools norm --rm-dup all -Oz \
    -o "${OUT}.var_sizes.Smalls_list.no_mixed.INS.INS.vcf.gz"

bcftools view -v snps "${OUT}.var_sizes.Smalls_list.no_mixed.vcf.gz" \
  -Oz -o "${OUT}.var_sizes.Smalls_list.no_mixed.SNP.vcf.gz"

bcftools filter -i 'strlen(REF)==1 && strlen(ALT)==1' \
  "${OUT}.var_sizes.Smalls_list.no_mixed.SNP.vcf.gz" -Oz \
  | bcftools norm --rm-dup all -Oz \
    -o "${OUT}.var_sizes.Smalls_list.no_mixed.SNP.SNP.vcf.gz"
