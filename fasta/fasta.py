#!/usr/bin/env python3
import re, requests

bases = ['A', 'G', 'C', 'T']
base_pairs = {'A':'U', 'T':'A', 'G':'C', 'C':'G'}
codons = {'UUU':'F', 'UUC':'F', 'UUA':'L', 'UUG':'L', 'CUU':'L', 'CUC':'L', 'CUA':'L', 'CUG':'L',
          'AUU':'I', 'AUC':'I', 'AUA':'I', 'AUG':'M', 'GUU':'V', 'GUC':'V', 'GUA':'V', 'GUG':'V',
          'UCU':'S', 'UCC':'S', 'UCA':'S', 'UCG':'S', 'CCU':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
          'ACU':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCU':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
          'UAU':'Y', 'UAC':'Y', 'UAA':'_', 'UAG':'_', 'CAU':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
          'AAU':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'GAU':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
          'UGU':'C', 'UGC':'C', 'UGA':'_', 'UGG':'W', 'CGU':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
          'AGU':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', 'GGU':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'
         }
codes = {'A', 'C', 'G', 'T', 'U', 'R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 
         'V', 'N', '-', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'}
aa_only = {'E', 'F', 'I', 'L', 'P', 'Q', 'X', 'Z', '*'}

def parse_fasta(fasta):
    desc = ''
    seq = ''
    is_nucleic = True    # Changed to false if an aa only base is encountered
    is_valid = True      # Changed to false if an invalid character is encountered
    first = True         # To avoid yielding (desc, seq) when encountering the first description
    
    for line in fasta.splitlines():
        if line.find('>') != -1:
            seq += line[:line.find('>')].strip().upper()
            current = line[line.find('>') + 1:]
            
            if first:
                first = False
            else:
                if seq:
                    for c in seq:
                        if c not in codes and c not in aa_only:
                            yield(desc, 'Invalid character: ' + c)
                            is_valid = False
                        elif c in aa_only:
                            is_nucleic = False
                    # Only valid sequences need to have numerical digits replaced
                    if is_valid:
                        # Default assumption is that all sequences are nucleic
                        if is_nucleic:
                            seq = re.sub(r'\d', 'N', seq)
                        # If an amino acid-only base is found, the aa substitution for numerals is used
                        else:
                            seq = re.sub(r'\d', 'X', seq)

                    yield (desc, seq)

                    seq = ''
                    is_nucleic = True
                    is_valid = True

                else:
                    yield (desc, 'No sequence')
            
            # Handle the case when multiple descriptions are on a single line
            while current.find('>') != -1:
                desc = current[:current.find('>')]
                yield (desc, 'No sequence')
                current = current[current.find('>') + 1:]    
                
            desc = current
        else:
            seq += line.strip().upper()
    
    # yield the last sequence (there are no more '>' to trigger a yield)
    if seq:
        for c in seq:
            if c not in codes and c not in aa_only:
                yield(desc, 'Invalid character: ' + c)
                is_valid = False
            elif c in aa_only:
                is_nucleic = False
        # Only valid sequences need to have numerical digits replaced
        if is_valid:
            # Default assumption is that all sequences are nucleic
            if is_nucleic:
                seq = re.sub(r'\d', 'N', seq)
            # If an amino acid-only base is found, the aa substitution for numerals is used
            else:
                seq = re.sub(r'\d', 'X', seq)

        yield (desc, seq)

        seq = ''
        is_nucleic = True
        is_valid = True

    else:
        yield (desc, 'No sequence')

def ORFs_in_fasta(filename):
    # print('Determining ORFs in ' + filename)
    for pair in parse_fasta(filename):
        full_seq = {}
        ORFs = {}
        ORF_start = {}
        ORF_len = {}
        ORF_end = {}
        
        # Complement is read in reverse
        complement = pair[1][::-1]
        complement = [base_pairs[c] for c in complement]
        complement = ''.join(complement)

        for frame in [1, 2, 3]:
            index = frame - 1
            current_seq = []
            while index < len(pair[1])-2:
                current_codon = pair[1][index:index + 3].replace('T', 'U')
                current_seq.append(codons[current_codon])
                index += 3
            full_seq[frame] = ''.join(current_seq)
        
            # Find all ORFs, then the index at which they begin and end
            ORFs[frame] = re.findall('M[^_]*_', full_seq[frame])
            ORF_start[frame] = [(full_seq[frame].index(orf)) * 3 + frame for orf in ORFs[frame]]
            ORF_len[frame] = [len(orf) * 3 for orf in ORFs[frame]]
            ORF_end[frame] = [ORF_start[frame][i] + ORF_len[frame][i] - 1 for i in range(len(ORFs[frame]))]

        # Frames -1, -2, and -3 are done with respect to the reversed sequence complement
        for frame in [-1, -2, -3]:
            index = -1 * frame - 1
            current_seq = []
            while index < len(complement)-2:
                current_codon = complement[index:index + 3].replace('T', 'U')
                current_seq.append(codons[current_codon])
                index += 3
            full_seq[frame] = ''.join(current_seq)
        
            # Find all ORFs, then the index at which they begin and end
            ORFs[frame] = re.findall('M[^_]*_', full_seq[frame])
            ORF_start[frame] = [len(complement) - 1 - ((full_seq[frame].index(orf)) * 3 + frame) for orf in ORFs[frame]]
            ORF_len[frame] = [len(orf) * 3 for orf in ORFs[frame]]
            ORF_end[frame] = [ORF_start[frame][i] - ORF_len[frame][i] + 1 for i in range(len(ORFs[frame]))]
        
        ORF_text = []
        for frame in ORFs:
            for i in range(len(ORFs[frame])):
                ORF_text.append(f'* {frame} | {ORF_start[frame][i]} | {ORF_end[frame][i]} | {ORF_len[frame][i]} | {ORFs[frame][i]}')
        
        return '\n'.join(ORF_text)

def acc_to_fasta(acc):
        #assemble the esearch URL
        base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
        url = base + f"efetch.fcgi?db=nuccore&id={acc}" #"&WebEnv={web}"
        url += "&rettype=fasta&retmode=text"

        #post the efetch URL
        fasta = requests.get(url)
        return fasta.text
        # f = open(f'{acc}.fa', 'w')
        # f.write(fasta.text)
        # f.close()
        # print(f'The sequence has been written to {acc}.fa')

def acc_to_ORFs(acc):
    fasta = acc_to_fasta(acc)
    return ORFs_in_fasta(fasta)
