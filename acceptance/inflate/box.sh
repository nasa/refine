
rm -f box.egads
serveCSM -batch box.csm | tee box-csm.txt

rm -f box-vol.meshb
ref bootstrap box.egads | tee box-boot.txt

rm -f box01.meshb
ref adapt box-vol.meshb \
	   -g box.egads \
	   -x box.meshb | tee box-crv.txt

