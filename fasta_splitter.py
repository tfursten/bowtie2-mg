import argparse
import os
import sys


def readfq(fp):  # this is a generator function
    """
    Generator function to yield iterable of all sequence/quality pairs as
    2-tuples. Adapted from https://github.com/lh3/readfq/blob/master/readfq.py,
    original C implementation is MIT licensed also.

    Should be able to handle multi-line FASTQ and FASTA, but this is untested.

    :param handle: An open file handle (in text mode) to iterate over.
    :return: Yields 3-tuples of headers, sequences and quality strings
    """

    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Chunk a FASTA file.')
    parser.add_argument('--to-split', '-s', required=True,
        help='Absolute path to Fasta file which needs to be split.')

    parser.add_argument('--output_dir', '-o', required=True,
        help='Absolute path to output directory of split files.')

    parser.add_argument('--num-gigs', '-n', required=True, type=int,
        help='Maximum number of GB per chunk.')

    args = parser.parse_args()

    db_filename = args.to_split
    max_chunk_bytes = args.num_gigs * 1024 * 1024 * 1024
    output_dir = args.output_dir

    db_size = os.stat(db_filename).st_size

    print('Reading from: {} ({} bytes), writing {}GB increments to {}'.format(
        db_filename, db_size, args.num_gigs, output_dir))

    current_split_number = 0
    base_output_path = os.path.join(output_dir, 'ntdbout{}.fasta')
    current_outfilename = base_output_path.format(current_split_number)
    current_out_file = open(current_outfilename, 'w')

    print('Writing to {}...'.format(current_outfilename))

    with open(db_filename, 'r') as db_file:
        for header, sequence, _ in readfq(db_file):
            current_out_file.write('>' + header + '\n' + sequence)

            output_size = os.stat(current_outfilename).st_size
            if output_size >= max_chunk_bytes:
                current_split_number += 1
                current_outfilename = base_output_path.format(
                    current_split_number)

                current_out_file.close()
                current_out_file = open(current_outfilename, 'w')
                print('Writing to {}...'.format(current_outfilename))

        current_out_file.close()
    print('All done!')
