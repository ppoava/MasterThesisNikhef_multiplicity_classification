#!/bin/bash
# root -b -q 'status_analysis_bb.C("output_hardbbbar_1e5.root","BplusBminus_hardbbbar_1e5.root")'
# root -b -q 'status_analysis_bb.C("output_hardbbbar_JUNCTIONS_1e4.root","BplusBminus_hardbbbar_JUNCTIONS_1e4.root")'
root -b -q 'status_analysis_bb.C("output_minbias_bb_1e6.root","BplusBminus_minbias_bb_1e6.root")'
# root -b -q 'status_analysis_bb.C("output_minbias_bb_JUNCTIONS_1e4.root","BplusBminus_minbias_bb_JUNCTIONS_1e4.root")'

# root -b -q 'status_analysis_cc.C("output_hardccbar_1e5.root","DplusDminus_hardccbar_1e5.root")'
# root -b -q 'status_analysis_cc.C("output_hardccbar_JUNCTIONS_1e4.root","DplusDminus_hardccbar_JUNCTIONS_1e4.root")'
# COULD BE CC MINIBIAS?
root -b -q 'status_analysis_cc.C("output_minbias_bb_1e6.root","DplusDminus_minbias_cc_1e6.root")'
# ---------------------
# root -b -q 'status_analysis_cc.C("output_minbias_cc_JUNCTIONS_1e4.root","DplusDminus_minbias_cc_JUNCTIONS_1e4.root")'

root -b -q 'status_analysis_bb.C("output_hardQCD_1e7.root","BplusBminus_hardQCD_1e7.root")'
# root -b -q 'status_analysis_bb.C("output_hardQCD_JUNCTIONS_1e4.root","BplusBminus_hardQCD_JUNCTIONS_1e4.root")'
root -b -q 'status_analysis_cc.C("output_hardQCD_1e7.root","DplusDminus_hardQCD_1e7.root")'
# root -b -q 'status_analysis_cc.C("output_hardQCD_JUNCTIONS_1e4.root","DplusDminus_hardQCD_JUNCTIONS_1e4.root")'