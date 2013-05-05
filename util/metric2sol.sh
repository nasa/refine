#!/usr/bin/env bash

if [ ! "$#" -eq 2 ];then
  echo usage: $0 file.metric file.sol
fi

cat > $2 <<EOF
MeshVersionFormatted 2

Dimension 3

SolAtVertices
`cat $1 | wc -l`
1 3

EOF

awk '{print $1, $2, $4, $3, $5, $6}' $1 >> $2

