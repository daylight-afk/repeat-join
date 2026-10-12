"""Synthetic algorithm smoke tests; no experimental reads or external aligners."""
import contextlib
import importlib.util
import io
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
import pysam
from Bio import SeqIO

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('repeat_join',ROOT/'repeat_join/repeat_join.py')
rj=importlib.util.module_from_spec(spec); spec.loader.exec_module(rj)

class JoinTests(unittest.TestCase):
    def run_join(self,path,marker=None):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'joined.fasta'
            a=SimpleNamespace(target=str(ROOT/'tests/fixtures/segments.fasta'),path=path,
                out=str(out),name='synthetic',rotate_to=marker,marker_bp=4)
            with contextlib.redirect_stdout(io.StringIO()): rj.cmd_join(a)
            return str(SeqIO.read(out,'fasta').seq)
    def test_repeat_path_and_reverse_complement(self):
        self.assertEqual(self.run_join('U1,R,-U2,R'),'ACGTTGGGTAGGATGGTA')
    def test_marker_rotation_preserves_sequence(self):
        self.assertEqual(self.run_join('U1,R,-U2,R',str(ROOT/'tests/fixtures/marker.fasta')),
                         'TTGGGTAGGATGGTAACG')
    def test_unknown_segment_rejected(self):
        with self.assertRaises(SystemExit): self.run_join('U1,missing')

class JunctionTests(unittest.TestCase):
    def count(self,cases):
        p=SimpleNamespace(min_aln=1000,min_mapq=20,end_window=500,adj_gap=2000)
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'synthetic.bam'
            with pysam.AlignmentFile(str(path),'wb',header={'HD':{'VN':'1.6'},'SQ':[{'SN':x,'LN':10000} for x in ['A','B']]}) as out:
                for n,c in enumerate(cases):
                    span=c.get('span',1200); gap=c.get('gap',0); length=2*span+gap
                    for side in [0,1]:
                        a=pysam.AlignedSegment(); a.query_name=f'synthetic{n}'
                        a.query_sequence='A'*length; a.reference_id=side
                        a.reference_start=c.get('start',8500) if side==0 else c.get('second_start',100)
                        a.mapping_quality=c.get('mapq',60)
                        a.flag=0 if side==0 else c.get('second_flag',2048)
                        a.cigartuples=[(0,span),(4,span+gap)] if side==0 else [(4,span+gap),(0,span)]
                        out.write(a)
            return dict(rj.count_junctions(str(path),{'A':10000,'B':10000},p))
    def test_supplementary_pairs_and_bin_proxy(self):
        hits=self.count([{'start':8500},{'start':8510},{'start':8800}])
        self.assertEqual(hits,{('A:3','B:5'):[8500,8510,8800]})
        self.assertEqual(len({s//200 for s in hits[('A:3','B:5')]}),2)
    def test_secondary_and_low_mapq_excluded(self):
        self.assertEqual(self.count([{'second_flag':256},{'mapq':19}]),{})
    def test_short_alignment_large_gap_and_internal_end_excluded(self):
        self.assertEqual(self.count([{'span':999},{'gap':2001},{'second_start':4000}]),{})

if __name__=='__main__': unittest.main()
