#!/bin/bash -x

./autogen.sh --without-SDK 

( cd `uname` && make distcheck )

