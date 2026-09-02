#!/usr/bin/env bash
#
# decon_split.sh
#
# Split a VCF produced by vg deconstruct by variant size and by variant type.
#
#   small variants   every alternative allele differs from the reference by
#                    less than 50 bp, and the site carries no allele of 50 bp
#                    or more
#   structural       every alternative allele differs from the reference by
#   variants         50 bp or more
#
# Sites whose alternative alleles mix substitutions and length-changing
# alleles are separated into "only_mixed" files and excluded from the
# insertion and deletion sets, so that each final file holds one unambiguous
# category.
#
# Usage: bash decon_split.sh <input.vcf.gz> [min_sv_length]

set -euo pipefail

INPUT=${1:?usage: bash decon_split.sh <input.vcf.gz> [min_sv_length]}
MINSV=${2:-50}

NAME=$(basename "${INPUT}")
DIR=$(dirname "${INPUT}")
OUT="${DIR}/${NAME%.vcf.gz}"

###############################################################################
# 1. Length difference between the reference allele and every alternative
#    allele of each site
###############################################################################

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${INPUT}" \
  | awk '{
      split($4, alts, ",")
      for (i in alts)
        print $1, $2, $3, alts[i],
              (length(alts[i]) > length($3)) ? length(alts[i]) - length($3)
                                             : length($3) - length(alts[i])
    }' > "${OUT}.var_sizes.txt"

###############################################################################
# 2. Sites holding at least one structural allele, and sites holding only
#    small alleles
###############################################################################

awk -v m="${MINSV}" '$5 >= m {print $1, $2}' "${OUT}.var_sizes.txt" \
  | sort -u -k1,1V -k2,2n > "${OUT}.var_sizes.SVs_list.txt"

awk -v m="${MINSV}" 'NR==FNR {exclude[$1,$2]; next}
     $5 < m && !(($1,$2) in exclude) {print $1, $2}' \
  "${OUT}.var_sizes.SVs_list.txt" "${OUT}.var_sizes.txt" \
  | sort -u -k1,1V -k2,2n > "${OUT}.var_sizes.Smalls_list.txt"

bcftools view -T <(tr ' ' '\t' < "${OUT}.var_sizes.SVs_list.txt") "${INPUT}" \
  -Oz -o "${OUT}.var_sizes.SVs_list.vcf.gz"

bcftools view -T <(tr ' ' '\t' < "${OUT}.var_sizes.Smalls_list.txt") "${INPUT}" \
  -Oz -o "${OUT}.var_sizes.Smalls_list.vcf.gz"

###############################################################################
# 3. Structural variants: keep only the sites where every alternative allele
#    is structural
###############################################################################

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${OUT}.var_sizes.SVs_list.vcf.gz" \
  | awk -v m="${MINSV}" '{
      split($4, alts, ","); all_large = 1
      for (i in alts)
        if ((length(alts[i]) - length($3)) < m && (length($3) - length(alts[i])) < m) {
          all_large = 0; break
        }
      if (all_large) print $1 "\t" $2
    }' | sort -u -k1,1V -k2,2n > "${OUT}.var_sizes.SVs_list.all_50_alts.txt"

bcftools view -T "${OUT}.var_sizes.SVs_list.all_50_alts.txt" \
  "${OUT}.var_sizes.SVs_list.vcf.gz" \
  -Oz -o "${OUT}.var_sizes.SVs_list.all_50_alts.vcf.gz"

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
  | bcftools view -T - "${OUT}.var_sizes.SVs_list.all_50_alts.vcf.gz" \
    -Oz -o "${OUT}.var_sizes.SVs_list.all_50_alts.no_mixed.vcf.gz"

bcftools view -H "${OUT}.var_sizes.SVs_list.all_50_alts.vcf.gz" \
  | awk -F'\t' '{
      split($5, alts, ","); has_snp = 0; has_indel = 0
      for (i in alts) {
        if (length(alts[i]) == length($4)) has_snp = 1
        if (length(alts[i]) != length($4)) has_indel = 1
      }
      if (has_snp && has_indel) print
    }' \
  | bcftools view -T - "${OUT}.var_sizes.SVs_list.all_50_alts.vcf.gz" \
    -Oz -o "${OUT}.var_sizes.SVs_list.all_50_alts.only_mixed.vcf.gz"

# structural deletions
bcftools view -H "${OUT}.var_sizes.SVs_list.all_50_alts.no_mixed.vcf.gz" \
  | awk -F'\t' '{
      split($5, alts, ","); all_deletion = 1
      for (i in alts) if (length(alts[i]) >= length($4)) { all_deletion = 0; break }
      if (all_deletion) print
    }' \
  | bcftools view -T - "${OUT}.var_sizes.SVs_list.all_50_alts.no_mixed.vcf.gz" \
    -Oz -o "${OUT}.var_sizes.SVs_list.all_50_alts.no_mixed.DEL.vcf.gz"

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
  | bcftools view -T - "${OUT}.var_sizes.SVs_list.all_50_alts.no_mixed.vcf.gz" \
    -Oz -o "${OUT}.var_sizes.SVs_list.all_50_alts.no_mixed.INS.vcf.gz"

bcftools filter -i 'strlen(ALT) > strlen(REF)' \
  "${OUT}.var_sizes.SVs_list.all_50_alts.no_mixed.INS.vcf.gz" -Oz \
  | bcftools norm --rm-dup all -Oz \
    -o "${OUT}.var_sizes.SVs_list.all_50_alts.no_mixed.INS.INS.vcf.gz"

###############################################################################
# 4. Small variants: deletions, insertions and substitutions
###############################################################################

bcftools view -H "${OUT}.var_sizes.Smalls_list.vcf.gz" \
  | awk -F'\t' '{
      split($5, alts, ","); has_snp = 0; has_indel = 0
      for (i in alts) {
        if (length(alts[i]) == length($4)) has_snp = 1
        if (length(alts[i]) != length($4)) has_indel = 1
      }
      if (!(has_snp && has_indel)) print
    }' \
  | bcftools view -T - "${OUT}.var_sizes.Smalls_list.vcf.gz" \
    -Oz -o "${OUT}.var_sizes.Smalls_list.no_mixed.vcf.gz"

bcftools view -H "${OUT}.var_sizes.Smalls_list.vcf.gz" \
  | awk -F'\t' '{
      split($5, alts, ","); has_snp = 0; has_indel = 0
      for (i in alts) {
        if (length(alts[i]) == length($4)) has_snp = 1
        if (length(alts[i]) != length($4)) has_indel = 1
      }
      if (has_snp && has_indel) print
    }' \
  | bcftools view -T - "${OUT}.var_sizes.Smalls_list.vcf.gz" \
    -Oz -o "${OUT}.var_sizes.Smalls_list.only_mixed.vcf.gz"

bcftools view -H "${OUT}.var_sizes.Smalls_list.no_mixed.vcf.gz" \
  | awk -F'\t' '{
      split($5, alts, ","); all_deletion = 1
      for (i in alts) if (length(alts[i]) >= length($4)) { all_deletion = 0; break }
      if (all_deletion) print
    }' \
  | bcftools view -T - "${OUT}.var_sizes.Smalls_list.no_mixed.vcf.gz" \
    -Oz -o "${OUT}.var_sizes.Smalls_list.no_mixed.DEL.vcf.gz"

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
  | bcftools view -T - "${OUT}.var_sizes.Smalls_list.no_mixed.vcf.gz" \
    -Oz -o "${OUT}.var_sizes.Smalls_list.no_mixed.INS.vcf.gz"

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
