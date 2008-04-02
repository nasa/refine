#!/bin/bash -x

./autogen.sh

( cd `uname` && make distcheck )

