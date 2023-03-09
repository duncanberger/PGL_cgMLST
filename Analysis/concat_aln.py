#! /usr/bin/env python

import sys
import argparse
import re
from collections import OrderedDict
from Bio import AlignIO

def main(argv, out):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('--input', nargs="*", help="FASTA alignment files", required=True)
	parser.add_argument('--gap', help="Length of gaps between concatenated alignments [0]")
	parser.add_argument('--delimiter', help="Sequence name delimeter [None]")
	parser.add_argument('--output', help="Supermatrix output", required=True)
	args = parser.parse_args(argv)

	cumulativelength = 0
	aln_proc = 0
	aln_out = OrderedDict() 
	alignments = args.input

	for alignfile in alignments:
		aln_proc += 1
		alignment = AlignIO.read(alignfile, "fasta")
		alignmentlength = alignment.get_alignment_length()
		targetlength = cumulativelength + alignmentlength
		existingkeys = set(aln_out.keys())

		for seqrec in alignment:
			if args.delimiter:
				taxon_id = str(seqrec.id).split(args.delimiter)[0]
			else:
				taxon_id = str(seqrec.id)
			seq = aln_out.get(taxon_id)
			if seq:
				existingkeys.remove(taxon_id)
				aln_out[taxon_id] = seq + str(seqrec.seq)
			else:
				aln_out[taxon_id] = "-" * cumulativelength + str(seqrec.seq)
			if args.gap:
				aln_out[taxon_id] += str(int(args.gap) * "N")

		for notfoundkey in existingkeys:
			aln_out[notfoundkey] += "-" * alignmentlength
			if args.gap:
				aln_out[notfoundkey] += str(int(args.gap) * "N")
		cumulativelength += alignmentlength
		
	print("Alignments: {}\nTotal taxa: {}\n\nFinished".format(aln_proc, len(aln_out)), file=sys.stderr)

	with open(args.output,'w') as sm:
		for taxon,sequence in aln_out.items():
			print(">{}\n{}".format(taxon, sequence) , file=sm )

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
