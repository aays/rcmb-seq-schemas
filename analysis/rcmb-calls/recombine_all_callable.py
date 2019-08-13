'''
recombine_all_callable.py - 'simulate' recombination events
with all callable sites

takes in a space-separated file that contains desired locations of recombination events:

type position length_bp
CO 200 0
CO-GC 800 40
NCO 500 100

corresponds to simple CO at 200 bp, CO-GC at 800 bp with 40 bp GC,
and NCO at 500 bp of length 100 bp

parent 1 will be all 'A' and parent 2 all 'T'
'''

import argparse
from tqdm import tqdm
import random
import numpy as np
from copy import deepcopy

def args():
    parser = argparse.ArgumentParser(
        description='', 
        usage='python3.5 script.py [options]')

    parser.add_argument('-f', '--fname', required=True,
                        type=str, help='File containing recombination events')
    parser.add_argument('-l', '--length', required=True,
                        type=int, help='Length of parental haps')
    parser.add_argument('-o', '--out', required=True,
                        type=str, help='File to write to')

    args = parser.parse_args()

    return args.fname, args.length, args.out


def parse_events(fname: str) -> dict:
    '''
    parses input file containing recombination events
    should be a space separated file with columns event, position, and length
    '''
    events = {}
    with open(fname, 'r') as f:
        for line in f:
            if not line.startswith('event'):
                event, position, length_bp = line.split(' ')
                position, length_bp = int(position), int(length_bp)
                events[position] = [event, length_bp]
    return events

def draw_from_hap(hap_in: np.ndarray, event: str):
    '''
    helper function for recombine() that takes in an input
    hap and event and alters genotypes
    if CO, will return single hap (since switch point is instant)
    but for NCO/CO-GC, will return two (at start of event + at the
    end of the tract)
    '''
    hap = deepcopy(hap_in)
    if event == 'CO':
        # draw haps to switch
        parent1, parent2 = np.where(hap == 'A'), np.where(hap == 'T')
        switch1, switch2 = random.choice(*parent1), random.choice(*parent2)
        # switch phases
        hap[switch1] = 'T'
        hap[switch2] = 'A'
        return hap
    elif event == 'NCO':
        hap2 = deepcopy(hap)
        to_switch = random.choice([0,1,2,3])
        if hap[to_switch] == 'A':
            hap[to_switch] = 'T'
        elif hap[to_switch] == 'T':
            hap[to_switch] = 'A'
        return hap, hap2
    elif event == 'CO-GC':
        # determine phase switches
        parent1, parent2 = np.where(hap == 'A'), np.where(hap == 'T')
        switch1, switch2 = random.choice(*parent1), random.choice(*parent2)
        hap2 = deepcopy(hap)
        # create GC tract
        to_switch = random.choice([switch1, switch2])
        if hap[to_switch] == 'A':
            hap[to_switch] = 'T'
        elif hap[to_switch] == 'T':
            hap[to_switch] = 'A'
        # switch phases after GC tract
        hap2[switch1] = 'T'
        hap2[switch2] = 'A'
        return hap, hap2

        
def recombine(input_events: dict) -> dict:
    '''
    uses input event file to record recombination events
    assuming parents are of hap A or hap T
    '''
    initial_hap = ['A', 'A', 'T', 'T']
    phase_tracker = {0: np.array(initial_hap)}
    for position in tqdm(sorted(input_events.keys())):
        # convert existing positions to np array
        current_positions = np.array(sorted(list(phase_tracker.keys())))
        # get hap at highest position < position of event
        max_pos = current_positions[np.where(current_positions < position)[0].max()]
        current_hap = phase_tracker[max_pos]
        event, length_bp = input_events[position]
        if event == 'CO':
            out_hap = draw_from_hap(current_hap, 'CO')
            phase_tracker[position] = out_hap
        elif event in ['NCO', 'CO-GC']:
            out_hap1, out_hap2 = draw_from_hap(current_hap, event)
            phase_tracker[position] = out_hap1
            phase_tracker[position + length_bp] = out_hap2

    return phase_tracker
            
def write_recombinants(phase_changes: dict, length: int, out: str):
    '''
    uses event dict output by recombine() to create output
    sequences in fasta format
    '''
    change_positions = sorted(phase_changes.keys())
    with open(out, 'w') as f:
        for hap in tqdm(range(4)):
            f.write('>{hap}'.format(hap=hap) + '\n')
            sequence = ''
            for i in range(len(change_positions) - 1):
                position = change_positions[i]
                next_position = change_positions[i + 1]
                current_base = phase_changes[position][hap]
                sequence += ''.join([current_base for i in range(next_position - position)])
            remaining = length - len(sequence)
            last_base = sequence[-1]
            sequence += ''.join([last_base for i in range(remaining)])
            f.write(sequence + '\n')
            
def main():
    fname, length, out = args()
    print('Parsing input file...')
    input_events = parse_events(fname)
    print('Determining phase changes...')
    phase_changes = recombine(input_events)
    print('Creating sequences and writing to file...')
    write_recombinants(phase_changes, length, out)
    print('Done.')
    print('Good job!')

if __name__ == '__main__':
    main()

        
