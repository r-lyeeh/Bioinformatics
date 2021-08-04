#!/bin/bash
query=$1
awk 'FNR==NR{a[$1]=$1;next}($1 in a){print $0}' UKBB_CAD_QC.fam ukb${query}.tab >${query}.t1
head -1 ukb${query}.tab >${query}.header
cat ${query}.header ${query}.t1 >${query}.t2
sed -i -e "s/\t/\n/g" ${query}.header
# AGE 21022: 21022.0.0 702 
# SEX known
# Ethic background 21000: 21000.0.0 21000.1.0 21000.2.0 690 691 692 
# BMI 21001: 21001.0.0 21001.1.0 21001.2.0 693 694 695
# Systolic blood pressure 4080: 4080.0.0 4080.0.1 4080.1.0 4080.1.1 4080.2.0 4080.2.1 189,190,191,192,193,194
awk '{print $1,$189,$190,$191,$192,$193,$194,$690,$691,$692,$693,$694,$702}' ${query}.t2 >${query}.t3
# Townsend 189: 158
# Cholesterol 30690 30690.0.0 30690.1.0 :20 21
# HDL cholesterol 30760 30760.0.0 30760.1.0 :34 35
# LDL cholesterol 30780 30780.0.0 30780.1.0 :38 39
# Triglycerides 30870 30870.0.0 30870.1.0 :56 57
# Diabetes 2443 2443.0.0 2443.1.0 2443.2.0 :536 537 538
# Smoking 20116 20116.0.0 20116.1.0 20116.2.0 :4462 4463 4464
# Lipoprotein A 30790 30790.0.0 30790.1.0 :20
# Ciagrettes 3456 3456.0.0 3456.1.0 3456.2.0 :616 617 618
# MedicationS 6153 6153.0.0 6153.0.1 6153.0.2 6153.0.3 6153.1.0 6153.1.1 6153.1.2 6153.1.3 6153.2.0 6153.2.1 6153.2.2 6153.2.3 :2988: 2999
# Father illness 20107 4247:4276 20107.0.0:20107.2.9
# Mother illness 20110 4332:4364 20110.0.0:20110.2.10
# Siblings illness 20111 4365:4400 20111.0.0:20111.2.11
awk -v b=2 -v e=7 'BEGIN{FS=OFS=" "} {for (i=b;i<=e;i++) printf "%s%s", $i, (i<e ? OFS : ORS)}' ${query}.t3
# difference between tab and csv file
# e.g medS 4117890
# e.g smoke 3091962
