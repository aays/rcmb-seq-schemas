
porting from `tetrad-draws` log:

## 7/8/2019

creating 'all callable' sites:
- create two 10 kb fasta records
    - parent 1 - all As
    - parent 2 - all Ts
- create simulated 'tetrad' from cross b/w parent 1 and 2
    - walk along chromosome and determine crossover/NCO locations using Poisson process
    - update four tetrad strings accordingly
    - could be faster if done in 'long' format before transposing?
- do 'read draws' from single individuals in tetrad
    - what SNP combinations indicate one form of recombination over another?


creating two 10 kb 'parents':

```python
>>> outname = 'data/rcmb-calls/parents_all_callable.fasta'
>>> d = dict.fromkeys(['parent1', 'parent2'], '')
>>> d['parent1'] = ''.join(['A' for i in range(10000)])
>>> d['parent2'] = ''.join(['T' for i in range(10000)])
>>> with open(outname, 'w') as f:
  2     for id, seq in d.items():
  3         f.write('>' + id + '\n')
  4         f.write(seq + '\n')
```

## 9/8/2019

'recombining' the parents

some fixed rcmb events -
- COs
    - 200 kb
    - 800 kb - this could be CO-GC
- NCOs
    - 500 kb - 100 bp long
    - 700 kb - 50 bp long

script should encode 'instructions' for where phase changes happen instead
of storing large strings in memory

at start - create four 'offspring' - assign 2 to parent 1 (A) and
2 to parent 2 (T)

- if CO
    - draw one 'A' offspring and one 'T' offspring and switch to other phase
- if CO-GC
    - encode phase change _n_ bases downstream of initial 'switch point' _s_
    - in distance b/w _s_ and _n_ draw one offspring and switch phase for that length
- if NCO-GC
    - draw one offspring and temporarily switch phase for _n_ bases



