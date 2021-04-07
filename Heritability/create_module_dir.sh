#!/bin/bash
while read p;do
	mkdir ${p}
	cp /home/pang/Desktop/MP_Johan_20200221/snp_module/snplist_module_${p} ${p}/.
done <module
