# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#
__author__ = 'zerodel'

import sys
import random

try:
    import gffutils
except ImportError:
    sys.exit(-1)

read_length = 100
def main(quant_input, gtf_db, fasta_file, fasta_circ, fasta_linear):
    with open(quant_input) as getAssignment:
        all_lines = [line for line in getAssignment.readlines() if not line.startswith("#")]



    db = gffutils.FeatureDB(gtf_db)
    isoforms_raw = [line.strip().split('\t')[0].strip('"') for line in all_lines]
    # check whether db has such a key
    isoforms = []
    for isoform in isoforms_raw:
        try:
            some_iso = db[isoform]
            isoforms.append(isoform)
        except gffutils.FeatureNotFoundError:
            print isoform
            continue

    lines = []
    linear_lines = []
    for isoform in isoforms:
        raw_seq = db[isoform].sequence(fasta_file).strip()

        if isoform.startswith("chr"):
            seq_circ = ">%s\n%s%s\n" % (isoform, raw_seq[-read_length:], raw_seq)
            lines.append(seq_circ)
        else:
            seq_linear = ">%s\n%s\n" % (isoform, raw_seq)
            lines.append(seq_linear)
            linear_lines.append(seq_linear)

    with open(fasta_circ, "w") as export_fasta:
        export_fasta.writelines(lines)

    with open(fasta_linear, "w") as export_linear:
        export_linear.writelines(linear_lines)


def find_your_gene(quant_input, gtf_db, gene_map):
    with open(quant_input) as getAssignment:
        all_lines = [line for line in getAssignment.readlines() if not line.startswith("#")]



    db = gffutils.FeatureDB(gtf_db)
    isoforms_raw = [line.strip().split('\t')[0].strip('"') for line in all_lines]
    # check whether db has such a key
    isoforms = []
    for isoform in isoforms_raw:
        try:
            some_iso = db[isoform]
            isoforms.append(isoform)
        except gffutils.FeatureNotFoundError:
            print isoform
            continue

    iso_gene_pairs = []
    for isoform in isoforms:
        try:
            gene_id = [gene for gene in db.parents(isoform, featuretype="gene")][0].id
            iso_gene_pairs.append((isoform,gene_id))
            if len(iso_gene_pairs)%100 == 0:
                print len(iso_gene_pairs)
        except Exception:
            print isoform

    with open(gene_map, "w") as export_mapping:
        export_mapping.write("Transcript\tgene\n")
        export_mapping.writelines(["%s\t%s\n" % paire for paire in iso_gene_pairs])


if __name__ == "__main__":
    fasta_file = "/Users/zerodel/workspace/ref/genome.fa"
    gtf_db = "/Users/zerodel/workspace/ref/465c.db"

    target_path = "/Users/zerodel/workspace/quants/assignment/all1000/"
    quant = target_path + "iso.txt"
    out_fasta = target_path + "iso.fa"
    out_linear_fasta = target_path + "linear.fa"
    out_gene_map = target_path + "genemap.txt"


    main(quant, gtf_db, fasta_file, out_fasta, out_linear_fasta)
    #find_your_gene(quant, gtf_db, out_gene_map)