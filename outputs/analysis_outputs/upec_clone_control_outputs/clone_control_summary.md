# UPEC Clone-Controlled MALDI/WGS Analysis

Bruker MALDI-linked isolates with nonzero peak features: 407
ST131 isolates: 55

## Subset AUCs

- is_ST131_binary / overall: n=407, class_1=55, AUC=0.932, status=ok
- ciprofloxacin_R / overall: n=350, class_1=53, AUC=0.755, status=ok
- ciprofloxacin_R / ST131_only: n=47, class_1=30, AUC=0.725, status=ok
- ciprofloxacin_R / non_ST131_only: n=303, class_1=23, AUC=0.644, status=ok
- ceftriaxone_R / overall: n=360, class_1=30, AUC=0.689, status=ok
- ceftriaxone_R / ST131_only: n=48, class_1=21, AUC=0.515, status=ok
- ceftriaxone_R / non_ST131_only: n=312, class_1=9, AUC=0.547, status=ok

## Lineage-Controlled AUCs

- ciprofloxacin_R controlled by is_ST131_binary: raw=0.755, background_only=0.702, centered=0.625, retention=1.00, valid_groups=2
- ciprofloxacin_R controlled by ST: raw=0.755, background_only=0.898, centered=0.564, retention=0.36, valid_groups=7
- ciprofloxacin_R controlled by phylogroup_rough: raw=0.755, background_only=0.515, centered=0.762, retention=0.93, valid_groups=4
- ceftriaxone_R controlled by is_ST131_binary: raw=0.689, background_only=0.733, centered=0.488, retention=1.00, valid_groups=2
- ceftriaxone_R controlled by ST: raw=0.689, background_only=0.825, centered=0.503, retention=0.41, valid_groups=5
- ceftriaxone_R controlled by phylogroup_rough: raw=0.689, background_only=0.612, centered=0.687, retention=0.89, valid_groups=4

## Interpretation

- A high ST131 AUC means MALDI peaks encode lineage strongly.
- If ciprofloxacin AUC drops after ST/ST131/phylogroup centering, part of the resistance signal is lineage/background-associated.
- If ciprofloxacin AUC remains above chance within ST131 or outside ST131, that is evidence for residual focal resistance-associated signal beyond the dominant clone.
- Exact-ST centering is the strictest control and may retain fewer isolates because many STs contain only susceptible or only resistant isolates.
