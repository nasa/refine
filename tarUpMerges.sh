#!/usr/bin/env sh
 
echo '#!/usr/bin/env sh' > send.sh
echo 'tar czvf refineMerge`date +%Y%m%d%H%M`.tgz \' >> send.sh
cvs -q up | grep -v -e '^?' | cut -d ' ' -f 2 | sed 's/.$/& \\/' >> send.sh
echo '' >> send.sh
chmod u+x send.sh

cat send.sh

echo "running tar:"
echo "------------"
./send.sh

rm send.sh