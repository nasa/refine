#!/bin/bash -x

./autogen.sh --without-SDK

( cd `uname` && make check GEO_SDK=none )

