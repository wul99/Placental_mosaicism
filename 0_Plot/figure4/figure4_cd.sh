#!/bin/bash

## figure 4c
plotProfile -m data/PhyloP30way.1k.50bp.matrix.gz -out PhyloP30way.1k.50bp.pdf --refPointLabel mut --plotType se --samplesLabel PhyloP30way --colors "#D74B34" "#F3AAA4" "#4BB0C8" "#E5B7A1" "#F6DCD8" "#C2D8E0" --regionsLabel "BP_shared" "PP_shared" "P_unique" "BP_random" "PP_random" "P_random" --dpi 600 --legendLocation lower-left --plotHeight 8 --plotWidth 8 --yMin 0


## figure 4d
plotProfile -m data/DR_score.1k.50bp.matrix.gz -out DR_score.1k.50bp.pdf --refPointLabel mut --plotType se --samplesLabel "Depletion rank" --colors "#D74B34" "#F3AAA4" "#4BB0C8" "#E5B7A1" "#F6DCD8" "#C2D8E0" --regionsLabel "BP_shared" "PP_shared" "P_unique" "BP_random" "PP_random" "P_random" --dpi 600 --legendLocation lower-left --plotHeight 8 --plotWidth 8 --yMin 0.42


## supplementary figure 5a
plotProfile -m data/GC_content.1k.50bp.matrix.gz -out GC_content/all.1k.50bp.pdf --refPointLabel mut --plotType se --samplesLabel GC_5base --colors "#D74B34" "#F3AAA4" "#4BB0C8" "#E5B7A1" "#F6DCD8" "#C2D8E0" --regionsLabel "BP_shared" "PP_shared" "P_unique" "BP_random" "PP_random" "P_random" --dpi 600 --legendLocation lower-left --plotHeight 8 --plotWidth 8 --yMin 37
